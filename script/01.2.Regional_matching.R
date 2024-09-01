setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
library(tidyverse)
library(data.table)
library(terra)
library(exactextractr)
library(Cairo)
library(fixest)
library(plm)
library(spdep)

dat_dir <- "D:/shaoyanchuan/data/"
source("script/Inequality_cal_function.R")

## DOSE inflation factors
Deflator <- fread("data/input/Deflator.csv")%>%
  dplyr::select(-cpi_2015,-PPP,-cpi_usd)

Deflator_mod_China <- Deflator%>%
  filter(country == "China")%>%
  dplyr::select(-country,-deflator_2015)%>%
  rename(Year = year)

County_sf_China <- readRDS(paste0(dat_dir,"Shp\\county_China\\county_shp.Rds"))%>%
  dplyr::select(-pop_2000, -pop_2010, -pop_2020)%>%
  rename(CountyID = PAC)

## Spatial information
Grid_info_China <- fread(paste0(dat_dir, "Shp/GS(2019)1822/Grid_Additional_Information_0.1g.csv"))%>%
  dplyr::select(X_Lon,Y_Lat,GridID,ProvID,CountyID,Region)

Pop_China <- fread(paste0(dat_dir,"Population/RU_Pop/Urban_Pop_0.1g_2000-2020_China.csv"))%>%
  rename(X_Lon = x, Y_Lat = y, Pop = pop_all)%>%
  dplyr::select(X_Lon,Y_Lat,Year,Pop)

## Absolute source contributions
PM_source_extend_China <- fread("result/Interpretation/SHAP_China_TAP/PM_Source_2001-2019_TAP_China.csv")%>%
  mutate_at(paste0("PM25_",sector), ~ifelse(.x < 0,0,.x))

## Function to label the data by population characteristic (now by GDP and urban population ratio)
rast_rename <- function(rast, name){
  names(rast) <- name
  return(rast)
}


County_data_label <- function(data, data_sf = NULL, yearseq = c(2001,2019)){
  ## population change during the period
  pop_change <- data%>%
    filter(Year%in%c(yearseq))%>%
    group_by(CountyID)%>%
    summarise(Popc = Pop[Year == yearseq[2]] - Pop[Year ==  yearseq[1]])%>%
    ungroup()

  ## multiyear average ranking: GDP per capita and urbanization
  multiyear_ave_rank <- data%>%
    group_by(CountyID)%>%
    summarise(GDP_USD = mean(GDP_USD),
              Pop_Uratio = mean(Pop_Uratio))%>%
    ungroup()%>%
    mutate(GDPQ = cut(GDP_USD, quantile(GDP_USD, seq(0,1,0.2)),
                      include.lowest = TRUE, labels = 1:5),
           Pop_UratioQ = cut(Pop_Uratio, quantile(Pop_Uratio, seq(0,1,0.2)),
                             include.lowest = TRUE, labels = 1:5), .keep = "unused")
  
  start_year_rank <- data%>%
    filter(Year ==  yearseq[1])%>%
    dplyr::select(CountyID, GDP_USD, Pop_Uratio)%>%
    mutate(GDPQ_SY = cut(GDP_USD, quantile(GDP_USD, seq(0,1,0.2)),
                         include.lowest = TRUE, labels = 1:5),
           Pop_UratioQ_SY = cut(Pop_Uratio, quantile(Pop_Uratio, seq(0,1,0.2)),
                                include.lowest = TRUE, labels = 1:5), .keep = "unused")
  ## convert data to sf object and join the label columns
  data_sf_join <- data_sf%>%
    dplyr::select(CountyID)%>%
    left_join(data)%>%
    na.omit()%>%
    st_make_valid()%>%
    arrange(CountyID, Year)%>%
    list(pop_change, multiyear_ave_rank, start_year_rank)%>%
    reduce(left_join)

  ## Remove county without neighboring counties in China
  nb_lst <- data_sf_join%>%
    filter(Year ==  yearseq[1])%>%
    poly2nb(., queen = TRUE, row.names = .$CountyID)
  ## Remove ID
  rm_id <- data_sf_join%>%
    filter(Year ==  yearseq[1])%>%
    .[card(nb_lst) == 0, "CountyID"]%>%
    st_drop_geometry()
  
  data_join <- mutate(data_sf_join, 
                      nb_ind = ifelse(CountyID%in%rm_id$CountyID,1,0))%>%
    st_drop_geometry()
  return(data_join)
}

Grid_data_label <- function(data, yearseq = c(2001,2019)){
  data_change <- data%>%
    dplyr::select(GridID,Year,Pop)%>%
    filter(Year%in%yearseq)%>%
    group_by(GridID)%>%
    summarise(Popc = Pop[Year == yearseq[2]] - Pop[Year == yearseq[1]])%>%
    ungroup()

  Ave_year_rank <- data%>%
    group_by(GridID)%>%
    summarise(GDP_USD = mean(GDP_USD))%>%
    ungroup()%>%
    mutate(GDPQ = cut(GDP_USD, quantile(GDP_USD, seq(0,1,0.2)),
                      include.lowest = TRUE, labels = 1:5), .keep = "unused")
  
  Start_year_rank <- data%>%
    filter(Year == yearseq[1])%>%
    dplyr::select(GridID, GDP_USD)%>%
    mutate(GDPQ_SY = cut(GDP_USD, quantile(GDP_USD, seq(0,1,0.2)),
                         include.lowest = TRUE, labels = 1:5), .keep = "unused")

  data_final <- data%>%
    list(data_change, Ave_year_rank, Start_year_rank)%>%
    reduce(left_join)

  return(data_final)
}

##  Extract data from administrative areas by population-weighted values
Admin_weight <- function(data, data_sf, decompvar, weightvar = "Pop",
                         proj = "WGS84", is_rast = FALSE){
  if(is_rast){
    rast <- data
  }else{
    rast <- data%>%
      # mutate_at(decompvar, ~.x*Pop)%>%
      dplyr::select(X_Lon,Y_Lat,all_of(decompvar),all_of(weightvar))%>%
      rast(crs = proj)
  }

  rast%>%
    exact_extract(data_sf, c("weighted_mean"), weights = rast[[weightvar]], stack_apply = TRUE)%>%
    rename_at(vars(everything()), ~ sub("weighted_mean.", "", .x))%>%
    dplyr::select(-all_of(weightvar))
}

EMI_data_read <- function(sector_sel = "residential", year = 2000,
                          species = c("BC","OC","NOx","SO2","CO","NH3","VOC","PM25"),
                          # species = c("BC","OC","NOx","SO2","CO","NH3","NMVOC"),
                          Inv_name = "CEDS"){
  dat_dir <- "D:/shaoyanchuan/data/"
  if(Inv_name == "MEIC"){
    ## MEIC emissions. Unit: kg/m2/s
    EMI_file <- paste0(dat_dir, "Air_pollution/Inventory/MEIC/GC_Monthly/",year,"/",sector_sel,".nc")
  
    species%>%
      purrr::map(~ rast(EMI_file)%>%
                   .[[paste0(.x,"_",1:12)]]%>%
                   mean()%>%
                   rast_rename(paste0(sector_sel,"_",.x)))%>%
      reduce(c)
  }else if(Inv_name == "CEDS"){
    ## Aggregated monthly CEDS emissions. Unit: kg/m2/s
    EMI_file <- paste0(dat_dir, "Air_pollution/Inventory/CEDS/Monthly/",year,"/",year,"_emission_g0.1_China.Rds")
    read_rds(EMI_file)%>%
      filter(Sector == sector_sel)%>%
      group_by(Lon, Lat)%>%
      summarise_at(all_of(species), mean)%>%
      ungroup()%>%
      rename_at(all_of(species), ~ paste0(sector_sel,"_",.x))%>%
      rast(crs = "WGS84")
  }
}

##======County level data in China (pcGDP)======
lu_v <- c("Total_grid","Urban_grid","Agriculture_grid","Forest_grid","Grassland_grid","Wetland_grid")
sector_v <- c("agriculture","industry","power","residential","transportation")
met_v <- c("u10","v10","t2m","blh","tcc","sp","e","lsp")

socio_v_china <- c("GDP_pc_from_CSMAR","NTL_all","Pop_all","Pop_U","Pop_O")
met_v_china <- c("u10","v10","t2m","blh","tcc","tcw","lsp","e")
ind_china <- paste0("is_predicted")
iv_china <- c(socio_v_china, lu_v, ind_china)
decomp_v <- c("PM25",paste0("PM25_",sector_v))

## extract the county-level source-specific PM2.5 and associated risks
PM_source_county_China <- 2001:2019%>%
  map_dfr(~ PM_source_extend_China%>%
            filter(Year == .x)%>%
            Admin_weight(County_sf_China, decomp_v)%>%
            mutate(Year = .x, CountyID = County_sf_China$CountyID, .before = 1))

## Emissions area assigned to county boundary by fraction of area
EMI_source_county_China <- 2001:2019%>%
  map_dfr(\(year){
    sector_v%>%
      purrr::map(~ EMI_data_read(sector = .x, year = year, Inv_name = "MEIC"))%>%
      reduce(c)%>%
      c(., cellSize(., unit = "m"))%>%
      Admin_weight(., County_sf_China, decompvar = names(.), weightvar = "area", is_rast = TRUE)%>%
      mutate(CountyID = County_sf_China$CountyID, Year = year, .before = 1)
  })

## Chinese panel data including GDP per capita
County_panel_China <- readRDS(paste0(dat_dir, "Economic/Chinese_panel/Annual/Chinese_county_panel_data.Rds"))%>%
  dplyr::select(PAC, County, City, Province, Year, PM25_allP, all_of(iv_china))%>%
  rename(CountyID = PAC, Pop = Pop_all)%>%
  left_join(Deflator_mod_China)%>%
  mutate(GDP_USD = GDP_pc_from_CSMAR/fx*100/deflator_usd/1e3, ## Thousand dollar
         GDP_USD_S2 = GDP_USD^2, Pop = Pop/1e3, ## Population: Thousands
         Pop_Uratio = Pop_U/(Pop_U + Pop_O), .keep = "unused")%>%
  na.omit()

County_join_China <- list(County_panel_China, 
                          PM_source_county_China,
                          EMI_source_county_China)%>%
  reduce(left_join)%>%
  na.omit()%>%
  ## rename VOC variables
  rename_at(paste0(sector_v,"_VOC"),~ gsub("_VOC","_NMVOC",.x))%>%
  filter(Year >= 2001, Year <= 2019)%>%
  arrange(CountyID, Year)

County_rejoin_China <- County_data_label(County_join_China, County_sf_China, yearseq = c(2001,2019))

fwrite(County_rejoin_China, paste0("result/Inequality/Economic_imapct_county_dataset_China.csv"))

##======Gridded GDP level and matching======
lu_v <- c("Urban_grid","Agriculture_grid","Forest_grid","Grassland_grid","Wetland_grid")
sector_v <- c("agriculture","industry","power","residential","transportation")

emi_species <- c("BC","OC","NOx","SO2","NH3","NMVOC","PM25")
emi_v_china <- c(paste0("agriculture_",emi_species),
                 paste0("residential_",emi_species),
                 paste0("industry_",emi_species),
                 paste0("power_",emi_species),
                 paste0("transportation_",emi_species))

socio_v_china <- c("GDP_USD","NTL_all","Pop_all","Pop_U","Pop_O")
met_v_china <- c("u10","v10","t2m","blh","tcc","tcw","lsp","e")
iv_china <- c(socio_v_china, lu_v, emi_v_china)

## Grid information
Grid_info_China <- fread(paste0(dat_dir, "Shp/GS(2019)1822/Grid_Additional_Information_0.1g.csv"))%>%
  dplyr::select(X_Lon,Y_Lat,GridID,ProvID,CountyID,Region)

## GDP data
Deflator_mod_China <- Deflator%>%
  filter(country == "China")%>%
  dplyr::select(-country,-deflator_2015)%>%
  rename(Year = year)

Grid_GDP_pc <- fread(paste0(dat_dir, "Economic/Chinese_Panel/Interpolation/Grid_interpolate/Grid_XGBoost_GDP_pc_0.1g_2000-2020.csv"))%>%
  dplyr::select(-Pop_all,-factor2)%>%
  left_join(Deflator_mod_China)%>%
  mutate(GDP_USD = GDP_pc/fx*100/deflator_usd/1e3, .keep = "unused") ## Thousand dollar


Grid_join_cal_China <- function(year){
  cat("======Year: ",year,"======\n")
  ##  Land use fraction
  LU_file <- paste0(dat_dir, "Others/LU-ESA CCI/Urban_Aggre_0.1g/ESA_CCL_Class_0.1g_China_",year,".tif")
  LU <- rast(LU_file)

  ##  Population data, including urban and rural
  Pop_file <- paste0(dat_dir, "Population/RU_Pop/ESA-CCL/Urban_Pop_0.00833g_China_",year,".tif")
  Pop <- rast(Pop_file)%>%
    terra::resample(LU, method = "sum") ## resample at 0.1g
  names(Pop) <- c("Pop_all","Pop_U","Pop_O")

  ##  NTL data
  NTL_file <- paste0(dat_dir, "Economic/Nighttime_Light/extended time series/ESA-CCL/Urban_NTL_0.004g_China_",year,".tif")
  NTL <- rast(NTL_file)%>%
    terra::resample(LU, method = "sum")

  ##  Emission data
  EMI <- sector_v%>%
    purrr::map(~ EMI_data_read(sector = .x, year = year, Inv_name = "MEIC"))%>%
    reduce(c)%>%
    terra::resample(LU, method = "average")

  Grid_df <- c(LU,Pop,NTL,EMI)%>%
    as.data.frame(xy = TRUE)%>%
    rename(X_Lon = x, Y_Lat = y)%>%
    rename_at(paste0(sector_v,"_VOC"),~ gsub("_VOC","_NMVOC",.x))%>%
    mutate_at(c("X_Lon","Y_Lat"), ~ round(.x,2))%>%
    left_join(Grid_info_China, .)%>%
    na.omit()%>%
    mutate(Year = year, .before = 1)%>%
    mutate_at(all_of(lu_v), ~.x/Total_grid)%>%
    dplyr::select(-Total_grid)

  Grid_GDP_pc%>%
    filter(Year == year)%>%
    left_join(Grid_df)%>%
    dplyr::select(X_Lon,Y_Lat,GridID,Year,all_of(iv_china))
}

Grid_join_China <- 2001:2019%>%
  map_dfr(~ Grid_join_cal_China(.x))

Grid_rejoin_China <- Grid_join_China%>%
  rename(Pop = Pop_all)%>%
  left_join(Grid_info_China,.)%>%
  na.omit()%>%
  Grid_data_label()%>%
  left_join(PM_source_risk_China%>%dplyr::select(-Pop))

fwrite(Grid_rejoin_China, paste0("result/Inequality/Economic_imapct_grid_dataset_China.csv"))

##======GDP per capita divisions by county and grid======
Prov_info <- fread(paste0("data/input/Province_information.csv"))%>%
  dplyr::select(ProvID,engname,Subarea)

County_info <- readRDS(paste0(dat_dir,"Shp/county_China/county_shp.Rds"))%>%
  st_drop_geometry()%>%
  rename(CountyID = PAC, ProvID = Province_code)%>%
  dplyr::select(CountyID, ProvID)%>%
  left_join(Prov_info)

County_rejoin_China <- fread(paste0("result/Inequality/Economic_imapct_county_dataset_China.csv"))
Grid_rejoin_China <- fread(paste0("result/Inequality/Economic_imapct_grid_dataset_China.csv"))
County_divide <- County_rejoin_China%>%
  dplyr::select(CountyID,GDPQ)%>%
  unique()
fwrite(County_divide, paste0("result/Inequality/County_division_2001-2019.csv"))

County_sr_divide <- County_rejoin_China%>%
  left_join(County_info)%>%
  filter(Subarea%in%c("N","E"),Year >= 2001, Year <= 2019)%>%
  group_by(CountyID)%>%
  summarise(GDP_USD = mean(GDP_USD))%>%
  ungroup()%>%
  GroupQ_mutate(group_name = "GDP_USD", group_quantile = seq(0,1,.2))%>%
  ungroup()%>%
  dplyr::select(-GDP_USD)%>%
  rename(GDPQ = GDP_USD_group)
fwrite(County_sr_divide, paste0("result/Inequality/County_division_CEC_area_2001-2019.csv"))

## Divide GDP per capita by single years
County_SY_divide <- County_rejoin_China%>%
  filter(Year >= 2001, Year <= 2019)%>%
  dplyr::select(CountyID,Year,GDP_USD)%>%
  group_by(Year)%>%
  GroupQ_mutate(group_name = "GDP_USD", group_quantile = seq(0,1,.2))%>%
  ungroup()%>%
  dplyr::select(-GDP_USD)%>%
  rename(GDPQ = GDP_USD_group)

fwrite(County_SY_divide, paste0("result/Inequality/County_single_year_division_2001-2019.csv"))

County_sr_SY_divide <- County_rejoin_China%>%
  left_join(County_info)%>%
  filter(Subarea%in%c("N","E"),Year >= 2001, Year <= 2019)%>%
  dplyr::select(CountyID,Year,GDP_USD)%>%
  group_by(Year)%>%
  GroupQ_mutate(group_name = "GDP_USD", group_quantile = seq(0,1,.2))%>%
  ungroup()%>%
  dplyr::select(-GDP_USD)%>%
  rename(GDPQ = GDP_USD_group)
fwrite(County_sr_SY_divide, paste0("result/Inequality/County_single_year_division_CEC_area_2001-2019.csv"))

## population division
Grid_pop_divide <- Grid_rejoin_China%>%
  filter(Year >= 2001, Year <= 2019)%>%
  group_by(GridID,ProvID,X_Lon,Y_Lat)%>%
  summarise(Pop = mean(Pop), GDP_USD = mean(GDP_USD))%>% ## multiyear average
  ungroup()%>%
  mutate(Pop_group = ifelse(Pop <= 500,1,2))%>%
  left_join(Prov_info)

## grid division by gdp per capita
Grid_divide <- Grid_pop_divide%>%
  group_by(Pop_group)%>%
  GroupQ_mutate(group_name = "GDP_USD", group_quantile = seq(0,1,.2))%>%
  ungroup()%>%
  dplyr::select(-Pop, -GDP_USD)%>%
  rename(GDPQ = GDP_USD_group)
fwrite(Grid_divide, paste0("result/Inequality/Grid_division_2001-2019.csv"))

Grid_sr_divide <- Grid_pop_divide%>%
  filter(Subarea%in%c("N","E"))%>%
  group_by(Pop_group)%>%
  GroupQ_mutate(group_name = "GDP_USD", group_quantile = seq(0,1,.2))%>%
  ungroup()%>%
  dplyr::select(-Pop, -GDP_USD)%>%
  rename(GDPQ = GDP_USD_group)

fwrite(Grid_sr_divide, paste0("result/Inequality/Grid_division_CEC_area_2001-2019.csv"))

## divide by gdp per capita in a single year
Grid_SY_divide <- Grid_rejoin_China%>%
  filter(Year >= 2001, Year <= 2019)%>%
  dplyr::select(GridID,X_Lon,Y_Lat,Year,Pop,GDP_USD)%>%
  mutate(Pop_group = ifelse(Pop <= 500,1,2))%>%
  group_by(Year,Pop_group)%>%
  GroupQ_mutate(group_name = "GDP_USD", group_quantile = seq(0,1,.2))%>%
  ungroup()%>%
  dplyr::select(-Pop, -GDP_USD)%>%
  rename(GDPQ = GDP_USD_group)
fwrite(Grid_SY_divide, paste0("result/Inequality/Grid_single_year_division_2001-2019.csv"))

Grid_sr_SY_divide <- Grid_rejoin_China%>%
  left_join(Prov_info)%>%
  filter(Year >= 2001, Year <= 2019,Subarea%in%c("N","E"))%>%
  dplyr::select(GridID,X_Lon,Y_Lat,Year,Pop,GDP_USD)%>%
  mutate(Pop_group = ifelse(Pop <= 500,1,2))%>%
  group_by(Year,Pop_group)%>%
  GroupQ_mutate(group_name = "GDP_USD", group_quantile = seq(0,1,.2))%>%
  ungroup()%>%
  dplyr::select(-Pop, -GDP_USD)%>%
  rename(GDPQ = GDP_USD_group)
fwrite(Grid_sr_SY_divide, paste0("result/Inequality/Grid_single_year_division_CEC_area_2001-2019.csv"))
