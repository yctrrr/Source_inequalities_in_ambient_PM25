setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
library(tidyverse)
library(data.table)
library(terra)
library(exactextractr)
library(Cairo)
library(fixest)
library(plm)
library(spdep)
library(factoextra)
library(ineq)
library(ggplot2)
library(magick)
source('D:/shaoyanchuan/codebook/function/Grid_visualization.R')
source('D:/shaoyanchuan/codebook/function/Health_calcu_core.R')
source('D:/shaoyanchuan/codebook/function/Time_series_analysis.R')
source("script/Inequality_cal_function.R")

dat_dir <- "D:/shaoyanchuan/data/"
sector <- c("agriculture","industry","power","residential","transportation")
metvar <- c("u10","v10","t2m","blh","tcc","sp","e","lsp")

base_year <- 2019
dsname <- "MEIC"
Prov_info <- fread(paste0("data/input/Province_information.csv"))%>%
  dplyr::select(ProvID,engname,Subarea)

County_sf_China <- readRDS(paste0(dat_dir,"Shp/county_China/county_shp.Rds"))%>%
  dplyr::select(-pop_2000, -pop_2010, -pop_2020)%>%
  rename(CountyID = PAC)

County_sf_rejoin_China <- fread("result/Inequality/Economic_imapct_county_dataset_China.csv", encoding = 'UTF-8')%>%
  left_join(County_sf_China%>%dplyr::select(CountyID,Province_code),.)%>%
  na.omit()%>%
  rename(ProvID = Province_code)%>%
  left_join(Prov_info)

Grid_rejoin_China <- fread(paste0("result/Inequality/Economic_imapct_grid_dataset_China.csv"))%>%
  left_join(Prov_info)

Pop_SY <- Grid_rejoin_China%>%
  filter(Year == base_year)%>%
  dplyr::select(X_Lon, Y_Lat, Subarea, Pop)

County_divide <- fread(paste0("result/Inequality/County_division_2001-2019.csv"))
County_sr_divide <- fread(paste0("result/Inequality/County_division_CEC_area_2001-2019.csv"))
County_SY_divide <- fread(paste0("result/Inequality/County_single_year_division_2001-2019.csv"))
County_sr_SY_divide <- fread(paste0("result/Inequality/County_single_year_division_CEC_area_2001-2019.csv"))

Grid_divide <- fread(paste0("result/Inequality/Grid_division_2001-2019.csv"))
Grid_sr_divide <- fread(paste0("result/Inequality/Grid_division_CEC_area_2001-2019.csv"))
Grid_SY_divide <- fread(paste0("result/Inequality/Grid_single_year_division_2001-2019.csv"))
Grid_sr_SY_divide <- fread(paste0("result/Inequality/Grid_single_year_division_CEC_area_2001-2019.csv"))

Grid_rejoin_China%>%
  na.omit%>%
  group_by(Year)%>%
  summarise_at(c(paste0("PM25_",sector),"PM25"), ~weighted.mean(.x,Pop))%>%
  ungroup()%>%
  as.data.frame()
Grid_rejoin_China%>%names

PM_source_GC_China <- fread(paste0("result/Interpretation/PM_Source_",dsname,"_",base_year,"_GEOS-Chem_China.csv"))%>%
  left_join(Pop_SY)

## load in shape files
load_shp()

##  Population-weighted source-specific disparity (%) calculation (mutate method)
#'@param group_overall the group columns used to calculate overall exposure (e.g. Year)
#'@param source_col the columns used to distinguish different groups
#'@param add_col the columns that provide additional information from the raw data
#'@param rename should certain categorical variables be renamed
#'@param overall_ineq should the overall inequality be calculated or source-specific inequality
Popw_source_ineq <- function(data, group_overall = NULL, source_col = "Sector",
                             rename = TRUE, overall_ineq = TRUE,
                             sector = c("agriculture","industry","power","residential","transportation")){
  ## Difference between source-specific and overall values
  Popw_change <- data%>%
    group_by_at(c(group_overall, source_col))%>%
    mutate(Popw_change = Popw - Popw[GDPQ == 0])%>%
    ungroup()

  if(overall_ineq == TRUE){
    ## overall inequality
    Popw_ineq <- Popw_change%>%
      group_by_at(group_overall)%>%
      mutate(Popw_overall = sum(Popw[GDPQ == 0]))%>% ## overall PM2.5
      ungroup()%>%
      mutate(Popw_ineq = Popw_change/Popw_overall*100)%>% ## inequality by percentage
      as.data.frame()%>%
      # dplyr::select(GDPQ,all_of(group_overall),
      #               all_of(source_col),Popw,Popw_ineq)%>%
      filter(GDPQ != 0)
  }else{
    ## source-specific inequality
    Popw_ineq <- Popw_change%>%
      group_by_at(c(group_overall, source_col))%>%
      mutate(Popw_ineq = Popw_change/Popw[GDPQ == 0]*100)%>%
      ungroup()%>%
      as.data.frame()%>%
      # dplyr::select(GDPQ,all_of(group_overall),
      #               all_of(source_col),Popw,Popw_ineq)%>%
      filter(GDPQ != 0)
  }
  
  if(rename == TRUE){
    Popw_ineq <- Popw_ineq%>%
      GDPQ_rename()%>%
      mutate(Sector = gsub("PM25_|_risk","",Sector)%>%
               str_to_title()%>%
               factor(.,levels = str_to_title(sector))
      )
  }else{
    Popw_ineq <- Popw_ineq
  }
  return(Popw_ineq)
}

EMI_unit_conv <- function(df, lonlat = c("X_Lon","Y_Lat"),
                          emivar = "agriculture_NH3", ts = "annual"){
  ## interval seconds
  if(ts == "annual"){
    interval_sec <- interval(paste0(year,"0101"),paste0(year,"1231"))%/%hours(1)*3600
  }

  scale_factor <- 1e3
  lonlat_rast <- df%>%dplyr::select(all_of(lonlat))%>%
    distinct()%>%mutate(ID = 1)%>%
    rast(crs = "WGS84")

  ## unit: tons
  cell_value <- cellSize(lonlat_rast, unit = "m")%>%
    as.data.frame(xy = TRUE)%>%
    setnames(c(lonlat,"area"))%>%
    mutate_at(lonlat, ~round(.x,2))%>%
    mutate(area = area*interval_sec/scale_factor)
  
  left_join(df, cell_value)%>%
    mutate_at(all_of(emivar), ~ .x*area)
}

##======Source contributions for GDP groups (County level)======
## China
look_up <- expand.grid(suffix = c(""), GDP_div = c("","GDP_single_year_"),
                       stringsAsFactors = FALSE)
look_up
for(i in 1:nrow(look_up)){
  suffix <- look_up[i,1]
  GDP_div <- look_up[i,2]

  if(GDP_div == ""){
    join_divide <- County_divide
    sr_join_divide <- County_sr_divide
  }else{
    join_divide <- County_SY_divide
    sr_join_divide <- County_sr_SY_divide
  }
  
  GDP_source_summarise_China <- 2001:2019%>%
    map_dfr(\(year) paste0("PM25_",sector, suffix)%>%
              map_dfr(~ County_sf_rejoin_China%>%
                        st_drop_geometry()%>%
                        dplyr::select(-GDPQ)%>%
                        left_join(join_divide)%>%
                        na.omit()%>%
                        filter(Year == year)%>%
                        local_popw_summarise(group_main = "GDPQ", 
                                             var = .x, prefix = FALSE)%>%
                        mutate(Sector = .x))%>%
              mutate(Year = year))

  GDP_source_summarise_sr_China <- 2001:2019%>%
    map_dfr(\(year) paste0("PM25_",sector, suffix)%>%
              map_dfr(~ County_sf_rejoin_China%>%
                        st_drop_geometry()%>%
                        dplyr::select(-GDPQ)%>%
                        left_join(sr_join_divide)%>%
                        na.omit()%>%
                        filter(Year == year)%>%
                        local_popw_summarise(group_main = "GDPQ", 
                                             var = .x, prefix = FALSE)%>%
                        mutate(Sector = .x))%>%
              mutate(Year = year))

  for(year in c(2001,2013,2017,2019)){
    GDP_source_ineq_plot_China <- GDP_source_summarise_China%>%
      filter(Year == year)%>%
      Popw_source_ineq(group_overall = NULL, source_col = "Sector")%>%
      Source_ineq_vis(x = "Sector", y = "Popw_ineq", fill = "Popw", facet = "GDPQ",
                      fill_title = expression(paste("Population-weighted PM"[2.5]," (μg/m"^"3",")")))

    CairoPNG(paste0("result/Inequality/County_analysis/",GDP_div,"County_Source_exposure",suffix,"_inequality_China_",year,".png"),height = 1500,width = 4000, res = 200)
    print(GDP_source_ineq_plot_China)
    dev.off()
  }

  GDP_source_ts_plot_China <- GDP_source_summarise_China%>%
    Popw_source_ineq(group_overall = "Year", source_col = "Sector")%>%
    Source_ineq_ts_vis(x = "Year", y = "Popw_ineq", source = "Sector", facet = "GDPQ")
  
  CairoPNG(paste0("result/Inequality/County_analysis/",GDP_div,"County_Source_exposure",suffix,"_inequality_time_series_China.png"),height = 1500,width = 5000, res = 200)
  print(GDP_source_ts_plot_China)
  dev.off()
  
  GDP_source_ts_sr_plot_China <- GDP_source_summarise_sr_China%>%
    Popw_source_ineq(group_overall = "Year", source_col = "Sector")%>%
    Source_ineq_ts_vis(x = "Year", y = "Popw_ineq", source = "Sector", facet = "GDPQ")
  
  CairoPNG(paste0("result/Inequality/County_analysis/",GDP_div,"County_Source_exposure",suffix,"_inequality_CEC_area_time_series_China.png"),height = 1500,width = 5000, res = 200)
  print(GDP_source_ts_sr_plot_China)
  dev.off()
}

##======Source contributions for GDP groups (Grid level)======
look_up <- expand.grid(suffix = c(""), GDP_div = c("","GDP_single_year_"),
                       stringsAsFactors = FALSE)
look_up
i = 1
for(i in 1:nrow(look_up)){
  suffix <- look_up[i,1]
  GDP_div <- look_up[i,2]

  if(GDP_div == ""){
    join_divide <- Grid_divide
    sr_join_divide <- Grid_sr_divide
  }else{
    join_divide <- Grid_SY_divide
    sr_join_divide <- Grid_sr_SY_divide
  }
  
  ## Inequality index for different sources derived from grids with population > 500
  GDP_source_summarise_China <- 2001:2019%>%
    map_dfr(\(year) paste0("PM25_",sector,suffix)%>%
              map_dfr(~ Grid_rejoin_China%>%
                        filter(Year == year)%>%
                        dplyr::select(-GDPQ)%>%
                        left_join(join_divide)%>%
                        na.omit()%>%
                        local_popw_summarise(group_main = "GDPQ",
                                             group_add = c("Pop_group","Year"),
                                             var = .x, prefix = FALSE)%>%
                        mutate(Sector = .x)))

  ## Inequality index for different sources derived from grids in specific regions
  GDP_source_summarise_sr_China <- 2001:2019%>%
    map_dfr(\(year) paste0("PM25_", sector, suffix)%>%
              map_dfr(~ Grid_rejoin_China%>%
                        filter(Year == year, Subarea%in%c("N","E"))%>%
                        dplyr::select(-GDPQ)%>%
                        left_join(sr_join_divide)%>%
                        na.omit()%>%
                        local_popw_summarise(group_main = "GDPQ",
                                             group_add = c("Pop_group","Year"),
                                             var = .x, prefix = FALSE)%>%
                        mutate(Sector = .x)))

  GDP_source_summarise_China%>%
    filter(Pop_group == 2, GDPQ == 0, Sector == "PM25_residential")
  8.747-8.969

  GDP_source_summarise_sr_China%>%
    filter(Pop_group == 2, GDPQ == 1, Sector == "PM25_residential")
  ## Res in 2019
  8.607-9.511
  8.61-10.35
  ## Res during 2001-2019
  20.657-9.511
  20.089-8.607
  
  # Grid_rejoin_China%>%
  #   mutate(check = PM25_agriculture_risk+PM25_industry_risk+PM25_power_risk+PM25_residential_risk+PM25_transportation_risk,
  #          check2 = PM25_agriculture+PM25_industry+PM25_power+PM25_residential+PM25_transportation,
  #          Per = check/PM25_risk,
  #          Per2 = check2/PM25)%>%
  #   dplyr::select(PM25_risk,check,Per,Per2)
  
  ## Box plot of inequality for different years
  for(year in c(2001,2013,2017,2019)){
    GDP_source_ineq_China <- GDP_source_summarise_China%>%
      filter(Pop_group == 2, Year == year)%>%
      Popw_source_ineq(group_overall = NULL, source_col = "Sector")%>%
      Source_ineq_vis(x = "Sector", y = "Popw_ineq", fill = "Popw", facet = "GDPQ",
                      fill_title = expression(paste("Population-weighted PM"[2.5]," (μg/m"^"3",")")))
    
    CairoPNG(paste0("result/Inequality/Grid_analysis/",GDP_div,"Source_exposure",suffix,"_inequality_China_",year,".png"),height = 1500,width = 4000, res = 200)
    print(GDP_source_ineq_China)
    dev.off()
  }
  
  ##  Time series change
  {
    GDP_source_ineq_China <- GDP_source_summarise_China%>%
      filter(Pop_group == 2)%>%
      Popw_source_ineq(group_overall = "Year", source_col = "Sector")

    GDP_source_ineq_ts_China <- GDP_source_ineq_China%>%
      {if(suffix == "") Source_ineq_ts_vis(., x = "Year", y = "Popw_ineq", breaks = c(2005,2015),
                                           source = "Sector", facet = "GDPQ")
        else Source_ineq_ts_vis(., x = "Year", y = "Popw_ineq", breaks = c(2005,2015), source = "Sector",
                                facet = "GDPQ", title = "Health Risk Disparity (%)")}
      # theme(axis.text.x = element_blank())
    GDP_source_ineq_China%>%
      dplyr::select(Sector, GDPQ, Year, Popw, Popw_change, Popw_ineq)%>%
      filter(Year == 2019)
    ## Res
    2.76-(-0.65)

    ## Ind
    6.98 - 2.38
    
    17.78-2.546-(10.92-0.9089)
    
    CairoPNG(paste0("result/Inequality/Grid_analysis/",GDP_div,"Source_exposure",suffix,"_inequality_time_series_China.png"),
             height = 1500, width = 4500, res = 180)
    print(GDP_source_ineq_ts_China)
    dev.off()
    
    GDP_source_ineq_sr_China <- GDP_source_summarise_sr_China%>%
      filter(Pop_group == 2)%>%
      Popw_source_ineq(group_overall = "Year", source_col = "Sector")
    
    GDP_source_ineq_sr_China%>%
      dplyr::select(Sector, GDPQ, Year, Popw, Pop, Popw_change, Popw_ineq)%>%
      filter(Year == 2019, Sector %in% c("Residential","Industry","Agriculture")
             )
    ## Res
    8.61 - 10.35
    20.09 - 8.61
    ## Ind
    11.96 - 11.03
    ## Agr
    10.90 - 12.72
    
    GDP_source_ineq_sr_China%>%
      dplyr::select(Sector, GDPQ, Year, Popw, Pop, Popw_change, Popw_ineq)%>%
      filter(Year %in% c(2001,2019), Sector %in% c("Power","Transportation"))%>%
      arrange(GDPQ,Sector)
    
    (-0.600-0.248)+(-0.125-1.534)
    # (2.527+0.707)-(0.272-2.74)
    
    GDP_source_ineq_sr_China%>%
      group_by(GDPQ,Year)%>%
      # group_by(GDPQ, Year)%>%
      summarise(Popw = sum(Popw), Popw_ineq = sum(Popw_ineq), Pop = sum(Pop)/5)%>%
      ungroup()%>%
      # filter(GDPQ == "Very Rich")%>%as.data.frame()
      # group_by(Year)%>%
      # summarise(Popw = weighted.mean(Popw, Pop))%>%
      # ungroup()%>%
      filter(Year == 2019)%>%as.data.frame()
    41.12-37.99
    3.13/41.12
    1.82/41.12
    
    GDP_source_ineq_China%>%
      filter(Year == 2019, Sector == "Industry")%>%
      summarise(Pop = sum(Pop))
    (241+118+89+70+63)/1254
    (315+120+86+65+58)/1371
    
    GDP_source_ineq_ts_sr_China <- GDP_source_ineq_sr_China%>%
      {if(suffix == "") Source_ineq_ts_vis(., x = "Year", y = "Popw_ineq", breaks = c(2005,2015), 
                                           source = "Sector", facet = "GDPQ") 
        else Source_ineq_ts_vis(., x = "Year", y = "Popw_ineq", breaks = c(2005,2015), source = "Sector",
                                facet = "GDPQ", title = "Health Risk Disparity (%)")}
      # theme(strip.background = element_blank())
    
    CairoPNG(paste0("result/Inequality/Grid_analysis/",GDP_div,"Source_exposure",suffix,"_inequality_CEC_area_time_series_China.png"),
             height = 1500, width = 4500, res = 180)
    print(GDP_source_ineq_ts_sr_China)
    dev.off()
  }
}

##======Emission contributions for GDP group at county level======
sector_sel <- "residential"
emi_species <- c("BC","OC","NOx","SO2","NH3","NMVOC","PM25")

# emi_species <- c("BC","OC","NOx","SO2","CO","NH3","VOC","PM25")
for(sector_sel in sector){
  emi_name <- paste0(sector_sel,"_",emi_species)
  colors <- RColorBrewer::brewer.pal(length(emi_name), "RdYlGn")
  colors[4] <- "#3182bd"
  
  names(County_sf_rejoin_China)
  
  EMI_by_GDP_county_China <- 2001:2019%>%
    map_dfr(\(year) emi_name%>%
              map_dfr(~ County_sf_rejoin_China%>%
                        st_drop_geometry()%>%
                        filter(Year == year)%>%
                        local_popw_summarise(group_main = "GDPQ", var = .x, prefix = FALSE)%>%
                        mutate(Sector = .x, Year = year))
    )

  EMI_by_GDP_county_ts_China <- EMI_by_GDP_county_China%>%
    Popw_source_ineq(group_overall = "Year", rename = FALSE, overall_ineq = FALSE)%>%
    GDPQ_rename()%>%
    Source_ineq_ts_vis(., colors = colors, title = "Emission Disparity (%)")
  
  CairoPNG(paste0("result/Inequality/County_analysis/Emission_exposure_inequality_",sector_sel,"_county_time_series_China.png"),height = 1500,width = 5000, res = 200)
  print(EMI_by_GDP_county_ts_China)
  dev.off()
}

##======Emission contributions for GDP groups at grid level======
Grid_divide <- fread(paste0("result/Inequality/Grid_division_2001-2019.csv"))
sector_sel <- "residential"
# emi_species <- c("BC","OC","NOx","SO2","CO","NH3","VOC","PM25")
emi_species <- c("BC","OC","NOx","SO2","NH3","NMVOC","PM25")
# emi_name <- c("agr_NH3","ind_SO2","pow_SO2","res_PM25","tra_NOx")

for(sector_sel in sector){
  if(sector_sel == "agriculture"){
    emi_name <- "agriculture_NH3"
  }else{
    emi_name <- paste0(sector_sel,"_",emi_species)
  }
  
  colors <- RColorBrewer::brewer.pal(length(emi_name), "RdYlGn")
  colors[4] <- "#3182bd"

  EMI_by_GDP_grid_China <- 2001:2019%>%
    map_dfr(\(year) emi_name%>%
              map_dfr(~ Grid_rejoin_China%>%
                        filter(Year == year)%>%
                        dplyr::select(-GDPQ)%>%
                        left_join(Grid_divide)%>%
                        na.omit()%>%
                        local_popw_summarise(group_main = "GDPQ",
                                             group_add = c("Pop_group","Year"),
                                             var = .x, prefix = FALSE)%>%
                        mutate(Sector = .x)))

  EMI_by_GDP_grid_sr_China <- 2001:2019%>%
    map_dfr(\(year) emi_name%>%
              map_dfr(~ Grid_rejoin_China%>%
                        filter(Year == year, Subarea%in%c("N","E"))%>%
                        dplyr::select(-GDPQ)%>%
                        left_join(Grid_sr_divide)%>%
                        na.omit()%>%
                        local_popw_summarise(group_main = "GDPQ",
                                             group_add = c("Pop_group","Year"),
                                             var = .x, prefix = FALSE)%>%
                        mutate(Sector = .x)))

  EMI_by_GDP_grid_ts_China<- EMI_by_GDP_grid_China%>%
    filter(Pop_group == 2)%>%
    Popw_source_ineq(group_overall = "Year",rename = FALSE, overall_ineq = FALSE)%>%
    GDPQ_rename()%>%
    Source_ineq_ts_vis(., colors = colors, title = "Emission Disparity (%)")
  
  CairoPNG(paste0("result/Inequality/Grid_analysis/Emission_exposure_inequality_",sector_sel,"_grid_time_series_China.png"),height = 1500,width = 5000, res = 200)
  print(EMI_by_GDP_grid_ts_China)
  dev.off()
  
  EMI_by_GDP_grid_sr_ts_China <- EMI_by_GDP_grid_sr_China%>%
    filter(Pop_group == 2)%>%
    Popw_source_ineq(group_overall = "Year",rename = FALSE)%>%
    GDPQ_rename()%>%
    Source_ineq_ts_vis(., colors = colors, title = "Emission Disparity (%)")
  
  CairoPNG(paste0("result/Inequality/Grid_analysis/Emission_exposure_inequality_",sector_sel,"_CEC_area_grid_time_series_China.png"),height = 1500,width = 5000, res = 200)
  print(EMI_by_GDP_grid_sr_ts_China)
  dev.off()
}

year = 2019
emi_species <- c("BC","OC","NOx","SO2","NH3","NMVOC","PM25")
# emi_species <- c("BC","OC","NOx","SO2","CO","NH3","NMVOC")
emivar <- c("agriculture_NH3",apply(expand.grid(sector%>%setdiff("agriculture"),emi_species), 1, paste, collapse="_"))
names(Grid_rejoin_China)
emi_df_long <- Grid_rejoin_China%>%
  filter(Subarea%in%c("N","E"))%>%
  dplyr::select(-GDPQ)%>%
  left_join(Grid_sr_divide)%>%
  na.omit()%>%
  EMI_unit_conv(lonlat = c("X_Lon","Y_Lat"),emivar = emivar, ts = "annual")%>%
  group_by(Year, GDPQ)%>%
  summarise_at(all_of(emivar),~ sum(.x)/1e3)%>%
  ungroup()%>%
  # separate_rows( !c("Year","GDPQ"), sep = "_")
  pivot_longer(cols = !c("Year","GDPQ"), names_to = c("Sector","species"),
               names_pattern = "([A-Za-z]+)_([A-Za-z0-9]+)")%>%
  pivot_wider(names_from = "species", values_from = "value")
fwrite(emi_df_long, "result/Inequality/Grid_analysis/Sectoral_emission_by_species.csv")
emi_df_long%>%filter(Sector == "residential")
split <- emi_df_long

##======Source contributions from GEOS-Chem for GDP groups (Grid level)======
year <- 2019

Grid_divide <- fread(paste0("result/Inequality/Grid_division_2001-2019.csv"))
GDP_source_GC_summarise_China <- paste0("PM25_",sector)%>%
    map_dfr(~ PM_source_GC_China%>%
              filter(Year == year)%>%
              left_join(Grid_divide)%>%
              na.omit()%>%
              local_popw_summarise(group_main = "GDPQ",
                                   group_add = c("Pop_group","Year"),
                                   var = .x, prefix = FALSE)%>%
              mutate(Sector = .x))
GDP_source_GC_summarise_China

GDP_source_GC_summarise_sr_China <- paste0("PM25_",sector)%>%
    map_dfr(~ PM_source_GC_China%>%
              filter(Year == year, Subarea%in%c("N","E"))%>%
              left_join(Grid_divide)%>%
              na.omit()%>%
              local_popw_summarise(group_main = "GDPQ",
                                   group_add = c("Pop_group","Year"),
                                   var = .x, prefix = FALSE)%>%
              mutate(Sector = .x))
GDP_source_GC_summarise_China
{
  GDP_source_GC_ineq_China <- GDP_source_GC_summarise_China%>%
    filter(Pop_group == 2, Year == year)%>%
    Popw_source_ineq(group_overall = NULL, source_col = "Sector")%>%
    Source_ineq_vis(x = "Sector", y = "Popw_ineq", fill = "Popw", facet = "GDPQ",
                    fill_title = expression(paste("Population-weighted PM"[2.5]," (μg/m"^"3",")")))+
    scale_y_continuous(breaks = c(-5,0,5))
  
  CairoPNG(paste0("result/Inequality/Grid_analysis/Source_GEOS-Chem_exposure_inequality_China_",year,".png"),height = 1500,width = 4000, res = 200)
  print(GDP_source_GC_ineq_China)
  dev.off()
  
  GDP_source_GC_ineq_sr_China <- GDP_source_GC_summarise_sr_China%>%
    filter(Pop_group == 2, Year == year)%>%
    Popw_source_ineq(group_overall = NULL, source_col = "Sector")%>%
    Source_ineq_vis(x = "Sector", y = "Popw_ineq", fill = "Popw", facet = "GDPQ",
                    fill_title = expression(paste("Population-weighted PM"[2.5]," (μg/m"^"3",")")))
  
  CairoPNG(paste0("result/Inequality/Grid_analysis/Source_GEOS-Chem_exposure_inequality_CEC_China_",year,".png"),height = 1500,width = 4000, res = 200)
  print(GDP_source_GC_ineq_sr_China)
  dev.off()
}

##======Source contributions for GDP groups by sensitivity tests (Grid level)======
ratio_limit_vect <- c(1,5,5)
search_dist_vect <- c(50e3,100e3,200e3)
for(v in 1:length(ratio_limit_vect)){
  ratio_limit <- ratio_limit_vect[v]
  search_dist <- search_dist_vect[v]
  oupt_dir <- paste0("result/Interpretation/SHAP_China_TAP/0.5g_deseason_",search_dist/1e3,"km_",ratio_limit,"_ratio_limit/")

  join_divide <- Grid_divide
  sr_join_divide <- Grid_sr_divide
  
  PM_source_extend_China_test <- fread(paste0(oupt_dir, "PM_Source_2001-2019_TAP_China.csv"))%>%
    dplyr::select(GridID,X_Lon,Y_Lat,Year,all_of(paste0("PM25_",sector)))

  ## Inequality index for different sources derived from grids with population > 500
  GDP_source_summarise_China_test <- 2001:2019%>%
    map_dfr(\(year) paste0("PM25_",sector)%>%
              map_dfr(~ Grid_rejoin_China%>%
                        filter(Year == year)%>%
                        dplyr::select(-GDPQ, -all_of(paste0("PM25_",sector)))%>%
                        left_join(join_divide)%>%
                        left_join(PM_source_extend_China_test)%>%
                        na.omit()%>%
                        local_popw_summarise(group_main = "GDPQ",
                                             group_add = c("Pop_group","Year"),
                                             var = .x, prefix = FALSE)%>%
                        mutate(Sector = .x)))
  
  ## Inequality index for different sources derived from grids in specific regions
  GDP_source_summarise_sr_China_test <- 2001:2019%>%
    map_dfr(\(year) paste0("PM25_", sector)%>%
              map_dfr(~ Grid_rejoin_China%>%
                        filter(Year == year, Subarea%in%c("N","E"))%>%
                        dplyr::select(-GDPQ, -all_of(paste0("PM25_",sector)))%>%
                        left_join(sr_join_divide)%>%
                        left_join(PM_source_extend_China_test)%>%
                        na.omit()%>%
                        local_popw_summarise(group_main = "GDPQ",
                                             group_add = c("Pop_group","Year"),
                                             var = .x, prefix = FALSE)%>%
                        mutate(Sector = .x)))

  ##  Time series change
  {
    GDP_source_ineq_China_test <- GDP_source_summarise_China_test%>%
      filter(Pop_group == 2)%>%
      Popw_source_ineq(group_overall = "Year", source_col = "Sector")
    
    GDP_source_ineq_ts_China_test <- GDP_source_ineq_China_test%>%
      {if(suffix == "") Source_ineq_ts_vis(., x = "Year", y = "Popw_ineq", breaks = c(2005,2015),
                                           source = "Sector", facet = "GDPQ")
        else Source_ineq_ts_vis(., x = "Year", y = "Popw_ineq", breaks = c(2005,2015), source = "Sector",
                                facet = "GDPQ", title = "Health Risk Disparity (%)")}
    # theme(axis.text.x = element_blank())

    
    CairoPNG(paste0("result/Inequality/Grid_analysis/Source_exposure_inequality_",
                    search_dist/1e3,"km_",ratio_limit,"_ratio_limit_time_series_China.png"),
             height = 1500, width = 4500, res = 180)
    print(GDP_source_ineq_ts_China_test)
    dev.off()
    
    GDP_source_ineq_sr_China_test <- GDP_source_summarise_sr_China_test%>%
      filter(Pop_group == 2)%>%
      Popw_source_ineq(group_overall = "Year", source_col = "Sector")
    
    Source_ineq_sr_China_sel <- GDP_source_ineq_sr_China_test%>%
      filter(Sector == "Residential", GDPQ == "Very Rich")%>%
      dplyr::select(-GDPQ, -Pop_group, -Sector, -Ave, -SD, -Pop)%>%
      mutate(ratio_limit = ratio_limit, search_dist = search_dist/1e3)
    
    if(v == 1){
      Source_ineq_sr_China_sel_join <- Source_ineq_sr_China_sel
    }else{
      Source_ineq_sr_China_sel_join <- full_join(Source_ineq_sr_China_sel_join, Source_ineq_sr_China_sel)
    }
    
    GDP_source_ineq_ts_sr_China_test <- GDP_source_ineq_sr_China_test%>%
      {if(suffix == "") Source_ineq_ts_vis(., x = "Year", y = "Popw_ineq", breaks = c(2005,2015), 
                                           source = "Sector", facet = "GDPQ") 
        else Source_ineq_ts_vis(., x = "Year", y = "Popw_ineq", breaks = c(2005,2015), source = "Sector",
                                facet = "GDPQ", title = "Health Risk Disparity (%)")}
    # theme(strip.background = element_blank())
    
    CairoPNG(paste0("result/Inequality/Grid_analysis/Source_exposure_inequality_",
                    search_dist/1e3,"km_",ratio_limit,"_ratio_limit_CEC_area_time_series_China.png"),
             height = 1500, width = 4500, res = 180)
    print(GDP_source_ineq_ts_sr_China_test)
    dev.off()
  }
}


ratio_limit_vect <- c(5,1,10,5,5)
search_dist_vect <- c(50e3,50e3,50e3,100e3,200e3)
v = 2
for(v in 1:length(ratio_limit_vect)){
  ratio_limit <- ratio_limit_vect[v]
  search_dist <- search_dist_vect[v]
  oupt_dir <- paste0("result/Interpretation/SHAP_China_TAP/0.5g_deseason_",search_dist/1e3,"km_",ratio_limit,"_ratio_limit/")
  
  join_divide <- Grid_divide
  sr_join_divide <- Grid_sr_divide
  
  PM_source_extend_China_test <- fread(paste0(oupt_dir, "PM_Source_2001-2019_TAP_China.csv"))%>%
    dplyr::select(GridID,X_Lon,Y_Lat,Year,all_of(paste0("PM25_",sector)))
  
  ## Inequality index for different sources derived from grids with population > 500
  GDP_source_summarise_China_test <- c(2001)%>%
    # c(2001,2013)%>%
    map_dfr(\(year) paste0("PM25_",sector)%>%
              map_dfr(~ Grid_rejoin_China%>%
                        filter(Year == year)%>%
                        dplyr::select(-GDPQ, -all_of(paste0("PM25_",sector)))%>%
                        left_join(join_divide)%>%
                        left_join(PM_source_extend_China_test)%>%
                        na.omit()%>%
                        local_popw_summarise(group_main = "GDPQ",
                                             group_add = c("Pop_group","Year"),
                                             var = .x, prefix = FALSE)%>%
                        mutate(Sector = .x)))

  ## Inequality index for different sources derived from grids in specific regions
  GDP_source_summarise_sr_China_test <- c(2001)%>%
    # c(2001,2013)%>%
    map_dfr(\(year) paste0("PM25_", sector)%>%
              map_dfr(~ Grid_rejoin_China%>%
                        filter(Year == year, Subarea%in%c("N","E"))%>%
                        dplyr::select(-GDPQ, -all_of(paste0("PM25_",sector)))%>%
                        left_join(sr_join_divide)%>%
                        left_join(PM_source_extend_China_test)%>%
                        na.omit()%>%
                        local_popw_summarise(group_main = "GDPQ",
                                             group_add = c("Pop_group","Year"),
                                             var = .x, prefix = FALSE)%>%
                        mutate(Sector = .x)))
  
  expo_change_overall <- GDP_source_summarise_China_test%>%
    filter(Pop_group == 2, GDPQ%in%c(0,5))%>%
    dplyr::select(GDPQ,Year,Sector,Popw)%>%
    pivot_wider(names_from = "Sector", values_from = "Popw")%>%
    mutate(region = "China", metric = "Popw", ratio_limit = ratio_limit, search_dist = search_dist, .after = "Year")%>%
    as.data.frame()
  
  expo_ineq_overall <- GDP_source_summarise_China_test%>%
    filter(Pop_group == 2)%>%
    Popw_source_ineq(group_overall = "Year", source_col = "Sector", rename = FALSE)%>%
    dplyr::select(GDPQ,Year,Sector,Popw_ineq)%>%
    filter(GDPQ == 5)%>%
    pivot_wider(names_from = "Sector", values_from = "Popw_ineq")%>%
    mutate(region = "China", metric = "Ineq", ratio_limit = ratio_limit, search_dist = search_dist, .after = "Year")%>%
    as.data.frame()
  
  expo_change_CEC <- GDP_source_summarise_sr_China_test%>%
    filter(Pop_group == 2, GDPQ%in%c(0,5))%>%
    dplyr::select(GDPQ,Year,Sector, Popw)%>%
    pivot_wider(names_from = "Sector", values_from = "Popw")%>%
    mutate(region = "CEC", metric = "Popw", ratio_limit = ratio_limit, search_dist = search_dist, .after = "Year")%>%
    as.data.frame()
  
  expo_ineq_CEC <- GDP_source_summarise_sr_China_test%>%
    filter(Pop_group == 2)%>%
    Popw_source_ineq(group_overall = "Year", source_col = "Sector", rename = FALSE)%>%
    dplyr::select(GDPQ,Year,Sector,Popw_ineq)%>%
    filter(GDPQ == 5)%>%
    pivot_wider(names_from = "Sector", values_from = "Popw_ineq")%>%
    mutate(region = "CEC", metric = "Ineq", ratio_limit = ratio_limit, search_dist = search_dist, .after = "Year")%>%
    as.data.frame()

  expo_change_tmp <- full_join(expo_change_overall, expo_change_CEC)
  expo_ineq_tmp <- full_join(expo_ineq_overall, expo_ineq_CEC)
  expo_df_tmp <-  full_join(expo_change_tmp, expo_ineq_tmp)
  
  if(v == 1){
    # expo_change <- expo_change_tmp
    # expo_ineq <- expo_ineq_tmp
    expo_df <- expo_df_tmp
  }else{
    # expo_change <- full_join(expo_change, expo_change_tmp)
    # expo_ineq <- full_join(expo_ineq, expo_ineq_tmp)
    expo_df <- full_join(expo_df, expo_df_tmp)
  }
}
expo_df
expo_df%>%
  arrange(region,search_dist,ratio_limit,Year,GDPQ)%>%
  dplyr::select(1:6,PM25_residential)

expo_df <- expo_df%>%
  mutate_at(paste0("PM25_",sector), ~round(.x,2)%>%as.character())
  
expo_ineq <- expo_df%>%
  filter(metric == "Ineq")%>%
  rename_at(paste0("PM25_",sector), ~paste0(.x,"_ineq"))%>%
  mutate(metric = "Popw")
expo_ineq

expo_change <- expo_df%>%
  filter(metric == "Popw")%>%
  left_join(expo_ineq)

for(i in sector){
  var <- paste0("PM25_",i)
  expo_change <- expo_change%>%
    mutate(!!sym(var) := ifelse(GDPQ == 0, !!sym(var), paste0(!!sym(var)," (",!!sym(paste0(var,"_ineq")),")")))
  
}

expo_change%>%
  dplyr::select(region,ratio_limit,search_dist,Year,GDPQ,paste0("PM25_",sector))%>%
  arrange(region,search_dist,ratio_limit,Year,GDPQ)%>%
  fwrite("result/Inequality/Grid_analysis/Source_exposure_inequality_sensitivity_test_for_richest_group.csv")

# expo_df%>%
#   arrange(region,ratio_limit,search_dist,Year,metric,GDPQ)%>%
#   fwrite("result/Inequality/Grid_analysis/Source_exposure_inequality_sensitivity_test_for_richest_group.csv")
