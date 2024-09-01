setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
library(deweather)
library(openair)
library(tidyverse)
library(lubridate)
library(terra)
library(purrr)
library(mgcv)
library(ranger)
library(furrr)
library(tictoc)
library(data.table)
library(treeshap)
library(shapviz)
library(fixest)
library(plm)
library(nlme)
source("script/SHAP_cal.R")
source("D:/shaoyanchuan/codebook/function/Grid_visualization.R")

dat_dir <- "D:/shaoyanchuan/data/"
metvar <- c("u10","v10","t2m","blh","tcc","sp","e","lsp")
sector <- c("agriculture","industry","power","residential","transportation")
speciesvar_ceds <- c("NOx","SO2","NMVOC","NH3","CO","BC","OC")
speciesvar_meic <- c("NOx","SO2","NMVOC","NH3","BC","OC","PM25")

emivar_ceds <- apply(expand.grid(speciesvar_ceds, sector), 1, paste, collapse="_")
emivar_meic <- apply(expand.grid(speciesvar_meic, sector), 1, paste, collapse="_")

##  Grid information
Grid_info_China <- fread(paste0(dat_dir,"Shp/GS(2019)1822/Grid_Information_0.1g.csv"))
Pop_China <- fread(paste0(dat_dir,"Population/RU_Pop/Urban_Pop_0.1g_2000-2020_China.csv"))%>%
  rename(X_Lon = x, Y_Lat = y, Pop = pop_all)%>%
  dplyr::select(X_Lon, Y_Lat, Year, Pop)

##======Grid-specific deweathering using regional ML======
## Extract coarse meteorological data to simulate regional transport (2.5 degree resolutions)
res <- 0.5
regional_met_China <- 2001:2019%>%
  map_dfr(~ fread(paste0(dat_dir,"ERA5/monthly_csv_",res,"X",res,"/era5-reanalysis-single-an-g",res,"-",.x,".csv")))%>%
  group_by(Year,Month,X_Lon,Y_Lat,GridID)%>%
  summarise_at(metvar,mean)%>%ungroup()

## Read in emission data of China (unit: kg/m2/s)
emi_meic_China <- 2001:2019%>%
  map_dfr(\(year) fread(paste0(dat_dir,"Air_pollution/Inventory/MEIC/Monthly/",year,"/",year,"_emission.csv")))%>%
  rename(X_Lon = Lon, Y_Lat = Lat, NMVOC = VOC)%>%
  dplyr::select(-PM10)

## regional emission data (2.5 degree)
regional_emi_meic_China <- 2001:2019%>%
  map_dfr(\(year) fread(paste0(dat_dir,"Air_pollution/Inventory/MEIC/Monthly/",year,"/",year,"_emission_g",res,".csv"))%>%
            filter(Sector %in%sector))%>%
  rename(X_Lon = Lon, Y_Lat = Lat, NMVOC = VOC)%>%
  dplyr::select(-PM10)

grid_meic_China <- distinct(regional_emi_meic_China, X_Lon, Y_Lat)%>%
  mutate(Grid = "MEIC")

regional_emi_ceds_China <- 2001:2019%>%
  map_dfr(\(year) fread(paste0(dat_dir,"Air_pollution/Inventory/CEDS/Monthly/",year,"/",year,"_emission_g",res,"_China.csv"))%>%
            filter(Sector %in%sector))%>%
  rename(X_Lon = Lon, Y_Lat = Lat, PM25 = PPM25)

regional_emi_meic_ceds_China <- regional_emi_ceds_China%>%
  distinct(X_Lon,Y_Lat)%>%
  left_join(grid_meic_China)%>%
  mutate(Grid = ifelse(is.na(Grid),"CEDS","MEIC"))%>%
  filter(Grid == "CEDS")%>%dplyr::select(-Grid)%>%
  left_join(regional_emi_ceds_China)%>%
  full_join(regional_emi_meic_China)%>%
  dplyr::select(X_Lon,Y_Lat,Year,Month,Sector,GridID,all_of(speciesvar_meic))%>%
  mutate_at(all_of(speciesvar_meic),~ifelse(is.na(.x),0,.x))%>%
  na.omit()

if(res == 0.5){
  regional_emi_meic_ceds_China <- regional_emi_meic_ceds_China%>%
    filter(!GridID %in% c(11476))
}

##======TAP in China with CEDS data======
##  Parallel setting
options(future.globals.maxSize = 10000 * 1024^3)
plan(multisession, workers = 18)
dsname <- "MEIC"
ratio_limit_vect <- c(5,1,10,5,5)
search_dist_vect <- c(50e3,50e3,50e3,100e3,200e3)

for(v in 1:length(ratio_limit_vect)){
  ratio_limit <- ratio_limit_vect[v]
  search_dist <- search_dist_vect[v]
  oupt_dir <- paste0("result/Interpretation/SHAP_China_TAP/0.5g_deseason_",search_dist/1e3,"km_",ratio_limit,"_ratio_limit")
  
  if(!dir.exists(oupt_dir)){
    dir.create(oupt_dir)
  }
  
  for(id in 1:100){
    cat(paste0("\n======Sample progress: ",id,"/100 ======\n"))
    PM_join <- readRDS(paste0(dat_dir,"Air_pollution/TAP/PM/Daily_join/",
                              "2000-2020_Daily_PM25_MET_",id,".Rds"))
    GridID_name <- unique(PM_join$GridID)

    emi_coord_sel <- Grid_info_China%>%
      filter(GridID%in%GridID_name)%>%
      dplyr::select(X_Lon,Y_Lat)
    
    emi_meic_China_sel <- emi_meic_China%>%
      filter(X_Lon < max(emi_coord_sel$X_Lon)+0.1,
             X_Lon > min(emi_coord_sel$X_Lon)-0.1,
             Y_Lat < max(emi_coord_sel$Y_Lat)+0.1,
             Y_Lat > min(emi_coord_sel$Y_Lat)-0.1)
    
    regional_emi_meic_ceds_China_sel <- regional_emi_meic_ceds_China%>%
      filter(X_Lon < max(emi_coord_sel$X_Lon)+5,
             X_Lon > min(emi_coord_sel$X_Lon)-5,
             Y_Lat < max(emi_coord_sel$Y_Lat)+5,
             Y_Lat > min(emi_coord_sel$Y_Lat)-5)
    
    regional_met_China_sel <- regional_met_China%>%
      filter(X_Lon < max(emi_coord_sel$X_Lon)+5,
             X_Lon > min(emi_coord_sel$X_Lon)-5,
             Y_Lat < max(emi_coord_sel$Y_Lat)+5,
             Y_Lat > min(emi_coord_sel$Y_Lat)-5)

    SHAP_pred_id <- GridID_name%>%
      furrr::future_map_dfr(\(gname){
        cat("======gname",gname,"======\n")
        PM_gts <- PM_join%>%
          mutate(Year = lubridate::year(Date), Month = month(Date),Day = lubridate::day(Date))%>%
          filter(Year <= 2019, Year >= 2001, GridID == gname)%>%
          rename(PM25 = PM2.5)%>%
          group_by(Year, Month)%>%
          summarise_at(c("PM25",metvar), mean)%>%
          ungroup()
        
        coords <- Grid_info_China%>%
          filter(GridID == gname)%>%
          as.numeric()%>%.[1:2]
        
        Grid_decomp(data = PM_gts, coords = coords, emi = emi_meic_China_sel,
                    regional_emi = regional_emi_meic_ceds_China_sel,
                    regional_met = regional_met_China_sel,
                    speciesvar = speciesvar_meic,
                    deseason = TRUE, ratio_limit = ratio_limit,
                    buff_distance = search_dist)%>%
          mutate(GridID = gname, .before = 1)%>%
          suppressMessages()
      }, .progress = TRUE)
    saveRDS(SHAP_pred_id, paste0(oupt_dir,"/PM_RRF_",dsname,"_SHAP_2001-2019_China_",id,".Rds"))
  }
  
  PM_SHAP_annual_China <- 1:100%>%
    map_dfr(~ read_rds(paste0(oupt_dir,"/PM_RRF_",dsname,"_SHAP_2001-2019_China_",.x,".Rds")))
  
  emivar_fm <- sector%>%
    map_dfr(~ data.frame(Sector = .x,
                         fm = paste(paste0(speciesvar_meic,"_",.x,"_SHAP"), collapse = "+")
                         ))

  speciesvar_fm <- speciesvar_meic%>%
    map_dfr(~ data.frame(Species = .x,
                         fm = paste(paste0(.x,"_",sector,"_SHAP"), collapse = "+")
                         ))%>%
    mutate(Species = ifelse(Species == "PM25","PPM25",Species))

  for(x in 1:nrow(emivar_fm)){
    PM_SHAP_annual_China <- PM_SHAP_annual_China%>%
      mutate(!!emivar_fm[x,1] := !!rlang::parse_expr(emivar_fm[x,2]))
  }
  
  for(x in 1:nrow(speciesvar_fm)){
    PM_SHAP_annual_China <- PM_SHAP_annual_China%>%
      mutate(!!speciesvar_fm[x,1] := !!rlang::parse_expr(speciesvar_fm[x,2]))
  }
  
  PM_SHAP_annual_China <- PM_SHAP_annual_China%>%
    mutate(EMI := !!rlang::parse_expr(paste(sector, collapse = "+")),
           MET := !!rlang::parse_expr(paste(paste0(metvar,"_SHAP"), collapse = "+"))
           # PM25 = PM25_rm + PM25_month_ave
    )

  fwrite(PM_SHAP_annual_China, paste0(oupt_dir,"/PM_RRF_",dsname,"_SHAP_2001-2019_China.csv"))
}


dsname <- "MEIC"
PM_SHAP_annual_China <- 1:100%>%
  map_dfr(~ read_rds(paste0("result/Interpretation/SHAP_China_TAP/PM_RRF_",dsname,"_SHAP_2001-2019_China_",.x,".Rds"))
          )

emivar_fm <- sector%>%
  map_dfr(~ data.frame(Sector = .x,
                       fm = paste(paste0(speciesvar_meic,"_",.x,"_SHAP"), collapse = "+")
                       ))

speciesvar_fm <- speciesvar_meic%>%
  map_dfr(~ data.frame(Species = .x,
                       fm = paste(paste0(.x,"_",sector,"_SHAP"), collapse = "+")
  ))%>%
  mutate(Species = ifelse(Species == "PM25","PPM25",Species))

for(x in 1:nrow(emivar_fm)){
  PM_SHAP_annual_China <- PM_SHAP_annual_China%>%
    mutate(!!emivar_fm[x,1] := !!rlang::parse_expr(emivar_fm[x,2]))
}

for(x in 1:nrow(speciesvar_fm)){
  PM_SHAP_annual_China <- PM_SHAP_annual_China%>%
    mutate(!!speciesvar_fm[x,1] := !!rlang::parse_expr(speciesvar_fm[x,2]))
}

PM_SHAP_annual_China <- PM_SHAP_annual_China%>%
  mutate(EMI := !!rlang::parse_expr(paste(sector, collapse = "+")),
         MET := !!rlang::parse_expr(paste(paste0(metvar,"_SHAP"), collapse = "+"))
         # PM25 = PM25_rm + PM25_month_ave
         )

fwrite(PM_SHAP_annual_China,paste0("result/Interpretation/SHAP_China_TAP/PM_RRF_",dsname,"_SHAP_2001-2019_China.csv"))

##======Absolute source contributions======
## Calculate absolute source contributions
year <- 2019
dsname <- "MEIC"

Pop_China <- fread(paste0(dat_dir,"Population/RU_Pop/Urban_Pop_0.1g_2000-2020_China.csv"))%>%
  rename(X_Lon = x, Y_Lat = y, Pop = pop_all)%>%
  dplyr::select(X_Lon, Y_Lat, Year, Pop)

PM_source_China <- fread(paste0(dat_dir, "CTM/geoschem/GC_output_China/source_impacts/GC_PM25_source_contribution_",
                                dsname,"_0.1g_",year,".csv"))%>%
  dplyr::select(X_Lon,Y_Lat,Sector,PM25_source)%>%
  rename(PM25 = PM25_source)%>%
  filter(Sector != "others")%>%
  mutate(Year = year)

PM_SHAP_China <- fread(paste0("result/Interpretation/SHAP_China_TAP/PM_RRF_",dsname,"_SHAP_2001-2019_China.csv"))%>%
  left_join(dplyr::select(Grid_info_China,X_Lon,Y_Lat,GridID),.)%>%
  filter(!is.na(GridID))

PM_source_extend_China <- PM_SHAP_China%>%
  dplyr::select(GridID,X_Lon,Y_Lat,Year,all_of(sector))%>%
  pivot_longer(sector, names_to = "Sector", values_to = "SHAP")%>%
  left_join(PM_source_China)%>%
  group_by(GridID, Sector)%>%
  mutate(PM25 = ifelse(Year==year, PM25, PM25[Year==year]+SHAP-SHAP[Year==year]))%>%
  ungroup()%>%
  dplyr::select(-SHAP)%>%
  pivot_wider(names_from = "Sector",values_from = "PM25",names_prefix = "PM25_")%>%
  left_join(dplyr::select(PM_SHAP_China, GridID, Year, PM25))%>%
  as.data.frame()

fwrite(PM_source_extend_China, "result/Interpretation/SHAP_China_TAP/PM_Source_2001-2019_TAP_China.csv")

PM_source_GC_China <- fread(paste0(dat_dir, "CTM/geoschem/GC_output_China/source_impacts/GC_PM25_source_contribution_MEIC_0.1g_",year,".csv"))%>%
  dplyr::select(X_Lon,Y_Lat,Sector,PM25_GC_source)%>%
  rename(PM25_GC = PM25_GC_source)%>%
  filter(Sector != "others")%>%
  mutate(Year = year)%>%
  pivot_wider(names_from = "Sector",values_from = "PM25_GC",names_prefix = "PM25_")
fwrite(PM_source_GC_China, paste0("result/Interpretation/PM_Source_",dsname,"_",year,"_GEOS-Chem_China.csv"))


##======Other tests======
## Sensitivity tests for SHAP calculations
PM_source_China <- fread(paste0(dat_dir, "CTM/geoschem/GC_output_China/source_impacts/GC_PM25_source_contribution_",
                                dsname,"_0.1g_",year,".csv"))%>%
  dplyr::select(X_Lon,Y_Lat,Sector,PM25_source)%>%
  rename(PM25 = PM25_source)%>%
  filter(Sector != "others")%>%
  mutate(Year = year)

dsname <- "MEIC"
ratio_limit_vect <- c(1,10,5,5)
search_dist_vect <- c(50e3,50e3,100e3,200e3)

for(v in 1:length(ratio_limit_vect)){
  ratio_limit <- ratio_limit_vect[v]
  search_dist <- search_dist_vect[v]
  oupt_dir <- paste0("result/Interpretation/SHAP_China_TAP/0.5g_deseason_",search_dist/1e3,"km_",ratio_limit,"_ratio_limit/")
  
  PM_SHAP_China_test <- fread(paste0(oupt_dir,"PM_RRF_",dsname,"_SHAP_2001-2019_China.csv"))%>%
    left_join(dplyr::select(Grid_info_China,X_Lon,Y_Lat,GridID),.)%>%
    filter(!is.na(GridID))

  PM_source_extend_China_test <- PM_SHAP_China_test%>%
    dplyr::select(GridID,X_Lon,Y_Lat,Year,all_of(sector))%>%
    pivot_longer(sector, names_to = "Sector", values_to = "SHAP")%>%
    left_join(PM_source_China)%>%
    group_by(GridID, Sector)%>%
    mutate(PM25 = ifelse(Year==year, PM25, PM25[Year==year]+SHAP-SHAP[Year==year]))%>%
    ungroup()%>%
    dplyr::select(-SHAP)%>%
    pivot_wider(names_from = "Sector",values_from = "PM25",names_prefix = "PM25_")%>%
    left_join(dplyr::select(PM_SHAP_China_test, GridID, Year, PM25))%>%
    as.data.frame()

  fwrite(PM_source_extend_China_test, paste0(oupt_dir,"PM_Source_2001-2019_TAP_China.csv"))
}

## population-weighted average
PM_source_extend_China <- fread(paste0("result/Interpretation/SHAP_China_TAP/PM_Source_2001-2019_TAP_China.csv"))

PM_source_extend_China%>%
  left_join(Pop_China)%>%
  mutate_at(c(paste0("PM25_",sector),"PM25"), ~ifelse(.x < 0,0,.x))%>%
  group_by(Year)%>%
  summarise_at(c(paste0("PM25_",sector),"PM25"), ~weighted.mean(.x,Pop))%>%
  ungroup()%>%
  as.data.frame()%>%na.omit()%>%
  fwrite("result/Interpretation/Population-weighted_souce_specific_PM2.5_China.csv")

PM_source_extend_China%>%
  left_join(Pop_China)%>%
  left_join(Grid_info_China)%>%
  group_by(Year,ProvID)%>%
  summarise_at(c(paste0("PM25_",sector),"PM25"), ~weighted.mean(.x,Pop))%>%
  ungroup()%>%
  as.data.frame()%>%na.omit()%>%
  fwrite("result/Interpretation/Population-weighted_souce_specific_PM2.5_by_province_China.csv")
