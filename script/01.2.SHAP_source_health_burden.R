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
source('D:/shaoyanchuan/codebook/function/Health_calcu_core.R')

dat_dir <- "D:/shaoyanchuan/data/"
metvar <- c("u10","v10","t2m","blh","tcc","sp","e","lsp")
sector <- c("agriculture","industry","power","residential","transportation")
Grid_info_China <- fread(paste0(dat_dir,"Shp/GS(2019)1822/Grid_Information_0.1g.csv"))

use_CR('MRBRT')
read_GBD(area = "China", combine = TRUE)

PM_source_extend_China <- fread("result/Interpretation/SHAP_China_TAP/PM_Source_2001-2019_TAP_China.csv")%>%
  mutate_at(paste0("PM25_",sector), ~ifelse(.x < 0,0,.x))
# PM_source_extend_China%>%
#   filter(Year == 2019)%>%
#   dplyr::select(X_Lon,Y_Lat,PM25)%>%rast%>%plot

Pop_China <- fread(paste0(dat_dir,"Population/RU_Pop/Urban_Pop_0.1g_2000-2020_China.csv"))%>%
  rename(X_Lon = x, Y_Lat = y, Pop = pop_all)%>%
  dplyr::select(X_Lon, Y_Lat, Year, Pop)

RR <- (RR_table$MEAN)%>%RR_std(format = "wide")

year = 2001
source_risk_calcu <- function(data, year = 2001){
  Mort_Age_tbl <- read_mort_age(year = year)%>%
    dplyr::select(X_Lon,Y_Lat,Year,endpoint,agegroup,Mort_Age)%>%
    pivot_wider(names_from = c("endpoint","agegroup"),
                values_from = "Mort_Age", names_prefix = "MA_")
  
  PM25_RR <- data%>%
    AP_risk(RR, Mort_Age_tbl)%>%
    rename(PM25 = concentration, PM25_risk = risk)%>%
    mutate(PM25 = as.numeric(PM25))%>%
    mutate_at(paste0("PM25_",sector),~ ifelse(.x < 0,0,.x))
  
  for(i in paste0("PM25_",sector)){
    PM25_RR[[paste0(i,"_risk")]] <- PM25_RR[[i]]/PM25_RR[["PM25"]]*PM25_RR[["PM25_risk"]]
  }
  PM25_RR
}

## PM2.5 risk calculation
PM_risk_China <- 2001:2019%>%
  map_df(function(year){
    data <- PM_source_extend_China%>%
      rename(concentration = PM25)%>%
      filter(Year == year)
    Pop_sel <- Pop_China%>%
      filter(Year == year)
    data%>%
      source_risk_calcu(year = year)%>%
      left_join(Pop_sel)
  })

fwrite(PM_risk_China, "result/Interpretation/SHAP_China_TAP/PM_Source_Health_Risk_2001-2019_TAP_China.csv")

