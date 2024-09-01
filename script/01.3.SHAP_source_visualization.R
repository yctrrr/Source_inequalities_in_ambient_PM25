setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
library(tidyverse)
library(data.table)
library(rlang)
library(stringr)
library(ggplot2)

dat_dir <- "D:/shaoyanchuan/data/"
source("script/Inequality_cal_function.R")
source('D:/shaoyanchuan/codebook/function/Grid_visualization.R')
source("D:/shaoyanchuan/codebook/function/Time_series_analysis.R")
source("D:/shaoyanchuan/codebook/function/Grid_visualization.R")

Grid_info_China <- fread(paste0(dat_dir, "Shp/GS(2019)1822/Grid_Information_0.1g.csv"))%>%
  dplyr::select(X_Lon,Y_Lat,GridID,ProvID,CountyID)
Prov_info <- fread(paste0("data/input/Province_information.csv"))%>%
  dplyr::select(ProvID,engname,Subarea)

load_shp()
speciesvar <- c("NOx","SO2","VOC","NH3","CO","PM25","PM10")
sector <- c("agriculture","industry","power","residential","transportation")
metvar <- c("u10","v10","t2m","blh","tcc","sp","e","lsp")
dsname <- "MEIC"

##  source simulations from extended SHAP values
##  remove the simulations in Taiwan
PM_source_extend_China <- fread("result/Interpretation/SHAP_China_TAP/PM_Source_2001-2019_TAP_China.csv")%>%
  left_join(Grid_info_China,.)%>%
  filter(!ProvID%in%710000)

##  source simulations from global burden of disease in 2017
PM_source_ratio_GBD_2017 <- fread(paste0(dat_dir, "CTM/GBD-MAPS/PM25_sectorial_contribution_ratio_China_2017.csv"))%>%
  mutate(Year = 2017)%>%
  left_join(Grid_info_China,.)%>%
  filter(!ProvID%in%710000)
colnames(PM_source_ratio_GBD_2017)[6:10] <- paste0('PM25_', colnames(PM_source_ratio_GBD_2017)[6:10])

PM_source_GBD_2017 <- PM_source_extend_China%>%
  filter(Year == 2017)%>%
  dplyr::select(X_Lon,Y_Lat,Year,PM25)%>%
  left_join(PM_source_ratio_GBD_2017)%>%
  mutate_at(paste0("PM25_",sector), ~.x*PM25)

##  source simulations by multiplying the total PM2.5 from GEOS-Chem
##  remove the simulations in Taiwan
base_year <- 2019
PM_source_GC_China <- fread(paste0("result/Interpretation/PM_Source_",dsname,"_",base_year,"_GEOS-Chem_China.csv"))%>%
  left_join(Grid_info_China,.)%>%
  filter(!ProvID%in%710000)

##  compare the difference between SHAP, GC and GBD simulations
PM_source_long_SHAP_China <- PM_source_extend_China%>%
  dplyr::select(-PM25)%>%
  filter(Year%in%c(2017,2019))%>%
  pivot_longer(paste0("PM25_",sector), names_to = "Sector", values_to = "PM25_SHAP")%>%
  mutate(PM25_SHAP = ifelse(PM25_SHAP < 0,0,PM25_SHAP))

PM_source_long_GC_China <- PM_source_GC_China%>%
  pivot_longer(paste0("PM25_",sector), names_to = "Sector", values_to = "PM25_GC")

PM_source_long_GBD_2017 <- PM_source_GBD_2017%>%
  dplyr::select(-PM25)%>%
  pivot_longer(paste0("PM25_",sector), names_to = "Sector", values_to = "PM25_GBD")

PM_SY <- PM_source_extend_China%>%
  filter(Year %in% c(2017,2019))%>%
  dplyr::select(X_Lon,Y_Lat,Year,PM25)

Pop_SY <- fread(paste0(dat_dir,"Population/RU_Pop/Urban_Pop_0.1g_2000-2020_China.csv"))%>%
  rename(X_Lon = x, Y_Lat = y, Pop = pop_all)%>%
  dplyr::select(X_Lon, Y_Lat, Year, Pop)%>%
  filter(Year%in%c(2017,2019))

PM_source_deviation_China <- list(PM_source_long_SHAP_China, PM_source_long_GC_China,
                                  PM_source_long_GBD_2017, PM_SY, Pop_SY)%>%
  reduce(left_join)%>%
  mutate(PM25_compare = ifelse(Year == 2017, PM25_GBD, PM25_GC),
         PM25_diff = PM25_SHAP - PM25_compare,
         PM25_SHAP_ratio = PM25_SHAP/PM25, 
         PM25_compare_ratio = PM25_compare/PM25,
         PM25_ratio_diff = PM25_SHAP_ratio - PM25_compare_ratio)

## calculate the R square across all grids
PM_source_deviation_table <- PM_source_deviation_China%>%
  group_by(Year,Sector)%>%
  summarise(R2 = cor(PM25_SHAP, PM25_compare)^2,
            error = mean(abs(PM25_diff)),
            ratio_R2 = cor(PM25_SHAP_ratio, PM25_compare_ratio),
            ratio_error = mean(abs(PM25_ratio_diff)),
            PM25_wt_SHAP = weighted.mean(PM25_SHAP, Pop),
            PM25_wt_compare = weighted.mean(PM25_compare, Pop))%>%
  ungroup()

##======Absolute source contributions=======
Abs_source_plot <- function(data, year = NULL, breaks = c(-Inf,1,3,5,10,15,20,+Inf),
                            colors, area = "China", text_table = NULL, include_total = FALSE,
                            title = expression(paste("PM"[2.5]," concentration (μg/m"^"3",")"))
                            ){
  labels <- label_gen(breaks)
  if(!is.null(year)){
    source_apportion <- data%>%
      filter(Year == year)
  }else{
    source_apportion <- data
  }
  
  if(include_total){
    source_apportion <- source_apportion%>%
      dplyr::select(X_Lon,Y_Lat,paste0("PM25_",sector),Total)
  }else{
    source_apportion <- source_apportion%>%
      dplyr::select(X_Lon,Y_Lat,paste0("PM25_",sector))
  }

  if(include_total){
    source_apportion_proj <- source_apportion%>%
      dfproj(toCRS = CRS_China, method = "average")%>%
      rename(X_Lon = x, Y_Lat = y)%>%
      pivot_longer(c(paste0("PM25_",sector),"Total"), names_to = "Sector")%>%
      mutate(Sector = gsub("*PM25_","",Sector)%>%str_to_title(),
             Sector = factor(Sector, levels = c("Total","Agriculture","Industry","Power","Residential","Transportation")))
    
  }else{
    source_apportion_proj <- source_apportion%>%
      dfproj(toCRS = CRS_China, method = "average")%>%
      rename(X_Lon = x, Y_Lat = y)%>%
      pivot_longer(paste0("PM25_",sector), names_to = "Sector")%>%
      mutate(Sector = gsub("*PM25_","",Sector)%>%str_to_title(),
             Sector = factor(Sector, levels = c("Agriculture","Industry","Power","Residential","Transportation")))
    
  }
  
  samod <- source_apportion_proj%>%
    mutate_at(c("value"), ~ cut(.,breaks,labels))

  samod_plot <- samod%>%
    grid_facet_plot(x = "X_Lon", y = "Y_Lat", value = "value",
                    facet = "Sector", facet_row = 2, facet_text = 35)%>%
    grid_plot_set(area = area, colors = colors,
                  title = title,  size_list = c(0.1, 1.5, 25, 28, 30),
                  panel_setting = "facet")+
    theme(panel.grid.major = element_blank(),
          # panel.border = element_rect(size = 1.5, colour = "black", fill = NA),
          )
  
  if(!is.null(text_table)){
    samod_plot <- samod_plot%>%
      grid_plot_add_text(text_table = text_table, cex = 8)
  }

  samod_nanhai_plot <- samod%>%
    split(f = .$Sector)%>%
    purrr::map(~ (.x%>%
                 grid_plot(x = "X_Lon", y = "Y_Lat", value = "value")%>%
                 grid_plot_set(area = area, colors = colors,
                               title = NULL,
                               panel_setting = "facet")+
                 theme(panel.grid.major = element_blank())+
                 labs(title = NULL))%>%
                 nanhai_zoom_plot(facet_df = data.frame(Sector = unique(.x$Sector)),
                                  border.line.size = 1.2)
               )

  samod_overall_plot <- (samod_plot+samod_nanhai_plot)
  return(samod_overall_plot)
}

Abs_source_facet_grid_plot <- function(data, year_seq, breaks = c(-Inf,1,3,5,10,15,20,+Inf),
                                       colors, area = "China", text_table = NULL,
                                       title = expression(paste("PM"[2.5]," concentration (μg/m"^"3",")"))){
  labels <- label_gen(breaks)

  source_apportion <- data%>%
    filter(Year %in% year_seq)

  source_apportion_proj <- year_seq%>%
    map_dfr(~ source_apportion%>%
              filter(Year == .x)%>%
              dplyr::select(X_Lon,Y_Lat,paste0("PM25_",sector))%>%
              dfproj(toCRS = CRS_China, method = "average")%>%
              rename(X_Lon = x, Y_Lat = y)%>%
              pivot_longer(paste0("PM25_",sector), names_to = "Sector")%>%
              mutate(Sector = gsub("*PM25_","",Sector)%>%str_to_title(),
                     Year = .x))

  samod <- source_apportion_proj%>%
    mutate_at(c("value"), ~ cut(.,breaks,labels))

  samod_plot <- samod%>%
    grid_facet_grid_plot(x = "X_Lon", y = "Y_Lat", value = "value",
                         row = "Year", col = "Sector", facet_text = 30)%>%
    grid_plot_set(area = area, colors = colors,
                  title = title,
                  panel_setting = "facet")+
    scale_x_continuous(breaks = c(80,100,120))+
    coord_sf(xlim = c(-2825586,2406965))+
    theme(panel.grid.major = element_blank(),
          axis.text = element_text(color = "black", size = 20))

  if(!is.null(text_table)){
    samod_plot <- samod_plot%>%
      grid_plot_add_text(text_table = text_table)
  }

  samod_nanhai_plot <- samod%>%
    split(f = ~ Sector+Year)%>%
    purrr::map(~ (.x%>%
                    grid_plot(x = "X_Lon", y = "Y_Lat", value = "value")%>%
                    grid_plot_set(area = area, colors = colors,
                                  title = NULL,
                                  panel_setting = "facet")+
                    theme(panel.grid.major = element_blank())+
                    labs(title = NULL))%>%
                 nanhai_zoom_plot(facet_df = data.frame(Sector = unique(.x$Sector),
                                                        Year = unique(.x$Year)),
                                  border.line.size = .6)
    )
  
  samod_overall_plot <- (samod_plot+samod_nanhai_plot)
  return(samod_overall_plot)
}

## China
breaks <- c(-Inf,1,3,5,10,15,25,35,+Inf)
year_seq <- c(2001,2013,2019)
labels <- label_gen(breaks)
colors <- rev(RColorBrewer::brewer.pal(length(labels), "Spectral"))

PM_samod_multiyear_China_plot <- PM_source_extend_China%>%
  mutate_at(c(paste0("PM25_",sector)),~ifelse(.x < 0,0,.x))%>%
  mutate(Total = PM25_agriculture+PM25_industry+PM25_power+PM25_residential+PM25_transportation)%>%
  Abs_source_facet_grid_plot(year_seq, breaks, colors, "China")

ggplot_oupt(PM_samod_multiyear_China_plot, filename = paste0("result/Source_change/China/PM25_source_contributions_multiyear.png"),
            export_type = "Cairo", plot_type = "local", whr = c(2500,2600,155))

year <- 2019
breaks <- c(-Inf,1,3,5,10,15,25,+Inf)
labels <- label_gen(breaks)
colors <- rev(RColorBrewer::brewer.pal(length(labels), "Spectral"))
for(year in c(2001:2019)){
  PM_samod_China_plot <- PM_source_extend_China%>%
    mutate_at(c(paste0("PM25_",sector)),~ifelse(.x < 0,0,.x))%>%
    mutate(Total = PM25_agriculture+PM25_industry+PM25_power+PM25_residential+PM25_transportation)%>%
    Abs_source_plot(year, breaks, colors, "China", include_total = TRUE)

  ggplot_oupt(PM_samod_China_plot, filename = paste0("result/Source_change/China/PM25_source_contributions_",year,".png"),
              export_type = "Cairo", plot_type = "local", whr = c(3500,1800,155))
}

PM_samod_GBD_2017_plot <- PM_source_GBD_2017%>%
  mutate_at(paste0("PM25_",sector), ~ifelse(.x < 0,0,.x))%>%
  mutate(Total = PM25_agriculture+PM25_industry+PM25_power+PM25_residential+PM25_transportation)%>%
  Abs_source_plot(2017, breaks, colors, "China")
ggplot_oupt(PM_samod_GBD_2017_plot, filename = paste0("result/Source_change/China/PM25_source_contributions_GBD_2017.png"),
            export_type = "Cairo", plot_type = "local", whr = c(3500,1800,155))

PM_samod_GC_plot <- PM_source_GC_China%>%
  mutate_at(paste0("PM25_",sector), ~ifelse(.x < 0,0,.x))%>%
  mutate(Total = PM25_agriculture+PM25_industry+PM25_power+PM25_residential+PM25_transportation)%>%
  Abs_source_plot(base_year, breaks, colors, "China")
ggplot_oupt(PM_samod_GC_plot, filename = paste0("result/Source_change/China/PM25_source_contributions_GEOS-Chem_",base_year,".png"),
            export_type = "Cairo", plot_type = "local", whr = c(3500,1800,155))

## Multiyear average PM25 during 2001-2019
PM_source_MYA_China <- PM_source_extend_China%>%
  group_by(X_Lon, Y_Lat, GridID)%>%
  summarise_at(paste0("PM25_",sector), mean)%>%
  ungroup()

PM_samod_MYA_China_plot <- PM_source_MYA_China%>%
  mutate_at(c(paste0("PM25_",sector)), ~ ifelse(.x < 0,0,.x))%>%
  mutate(Total = PM25_agriculture+PM25_industry+PM25_power+PM25_residential+PM25_transportation)%>%
  Abs_source_plot(year = NULL, breaks, colors, "China")

ggplot_oupt(PM_samod_MYA_China_plot, filename = paste0("result/Source_change/China/PM25_source_contributions_multiyear_average.png"),
            export_type = "Cairo", plot_type = "local", whr = c(3500,1800,155))


## Change in source contributions during 2001-2019
PM_source_change_China <- PM_source_extend_China%>%
  filter(Year%in%c(2001,2019))%>%
  mutate_at(paste0("PM25_",sector), ~ifelse(.x < 0,0,.x))%>%
  mutate(Total = PM25_agriculture+PM25_industry+PM25_power+PM25_residential+PM25_transportation)%>%
  group_by(X_Lon, Y_Lat, GridID)%>%
  summarise_at(paste0("PM25_",sector), ~.x[Year==2019]-.x[Year==2001])%>%
  ungroup()

breaks <- c(-Inf,-16,-12,-8,-6,-4,-2,-0.5,0.5,2,+Inf)
labels <- label_gen(breaks)
colors <- colorspace::sequential_hcl(length(labels)-3, palette = "agGrnYl")%>%
  c("#f0f0f0","#fed98e","#fd8d3c")
  # c("white","#feb24c","#d95f0e")

PM_source_change_plot <- PM_source_change_China%>%
  mutate(Total = PM25_agriculture+PM25_industry+PM25_power+PM25_residential+PM25_transportation)%>%
  Abs_source_plot(year = NULL, breaks, colors, "China",include_total = TRUE,
                  title = expression(paste("Change in PM"[2.5]," concentration (μg/m"^"3",")")))

ggplot_oupt(PM_source_change_plot, filename = paste0("result/Source_change/China/PM25_source_contributions_change_2001-2019.png"),
            export_type = "Cairo", plot_type = "local", whr = c(3500,1800,155))

##======Population-weighted change======
PW_source_PM <- fread(paste0("result/Interpretation/Population-weighted_souce_specific_PM2.5_China.csv"))%>%
  pivot_longer(cols = !Year, names_to = "Sector", values_to = "PM25")%>%
  filter(Sector!="PM25")%>%
  mutate(Sector = gsub("PM25_","",Sector)%>%str_to_title())

PW_source_PM_by_province <- fread(paste0("result/Interpretation/Population-weighted_souce_specific_PM2.5_by_province_China.csv"))%>%
  pivot_longer(cols = !c(Year,ProvID), names_to = "Sector", values_to = "PM25")%>%
  filter(Sector!="PM25")%>%
  mutate(Sector = gsub("PM25_","",Sector)%>%str_to_title())

colors <- RColorBrewer::brewer.pal(5, "RdYlGn")
colors[3] <- "#3182bd"
PW_source_PM_plot <- box_series_plot(PW_source_PM, x = "Year", y = "PM25", 
                                     fill = "Sector", errorbar = FALSE)+
  scale_fill_manual(name = NULL, values = alpha(colors, .8),
                    guide = guide_legend(frame.colour = NA, frame.linewidth = 0,
                                         title.position = "top", title.hjust = 0.5,
                                         byrow = TRUE, ticks = F))+
  labs(x = NULL, y = expression(paste("PM"[2.5]," exposure (μg/m"^"3",")")))+
  scale_y_continuous(position = 'right')+
  theme(legend.position = "right",
        legend.key.width = unit(1.5, "cm"),
        legend.key.height = unit(1.5, "cm"),
        legend.spacing.y = unit(1.2, 'cm'),
        legend.key.spacing.y = unit(1.2, 'cm'),
        legend.text = element_text(colour = "black", size = rel(2.5)),
        axis.text = element_text(colour = "black",size = rel(2.5)),
        axis.title = element_text(colour = "black",size = rel(2.8)), ##axis title size
        plot.margin = unit(c(5,5,5,10),"mm"))

ggplot_oupt(PW_source_PM_plot, filename = paste0("result/Source_change/China/PM25_source_contributions_time_series.png"),
            export_type = "Cairo", plot_type = "local", whr = c(4500,1800,140))

for(id in c("Beijing","Shanghai")){
  if(id == "Beijing"){
    provid <- 110000
  }else if (id == "Shanghai"){
    provid <- 310000
  }
  
  PW_source_PM_by_provice_plot <- PW_source_PM_by_province%>%
    filter(ProvID == provid)%>%
    box_series_plot( x = "Year", y = "PM25", fill = "Sector", errorbar = FALSE)+
    scale_fill_manual(name = NULL, values = alpha(colors, .8),
                      guide = guide_legend(frame.colour = NA, frame.linewidth = 0,
                                           title.position = "top", title.hjust = 0.5,
                                           byrow = TRUE, ticks = F))+
    labs(x = NULL, y = expression(paste("PM"[2.5]," exposure (μg/m"^"3",")")))+
    scale_y_continuous(position = 'right')+
    theme(legend.position = "right",
          legend.key.width = unit(1.5, "cm"),
          legend.key.height = unit(1.5, "cm"),
          legend.spacing.y = unit(1.2, 'cm'),
          legend.key.spacing.y = unit(1.2, 'cm'),
          legend.text = element_text(colour = "black", size = rel(2.5)),
          axis.text = element_text(colour = "black",size = rel(2.5)),
          axis.title = element_text(colour = "black",size = rel(2.8)), ##axis title size
          plot.margin = unit(c(5,5,5,10),"mm"))
  
  ggplot_oupt(PW_source_PM_by_provice_plot, filename = paste0("result/Source_change/China/PM25_source_contributions_time_series_",id,".png"),
              export_type = "Cairo", plot_type = "local", whr = c(4500,1800,140))
  
}

##======Compare with other data sets/studies======
##  Data derived from 02.1
Grid_rejoin_China <- fread(paste0("result/Inequality/Economic_imapct_grid_dataset_China.csv"))
##  Compare with GBD 2017
Grid_rejoin_China
##  Compare with population-weighted PM2.5 from different sources
Grid_rejoin_China%>%
  filter(Year%in%c(2001,2005,2010,2015,2017,2019))%>%
  group_by(Year)%>%
  summarise_at(c(paste0("PM25_",sector),"PM25"), ~ weighted.mean(.x,Pop))%>%
  mutate(anthro_sum = PM25_agriculture+PM25_residential+PM25_power+PM25_industry+PM25_transportation)%>%
  ungroup()%>%
  as.data.frame()

Pop_China <- fread(paste0(dat_dir,"Population/RU_Pop/Urban_Pop_0.1g_2000-2020_China.csv"))%>%
  rename(X_Lon = x, Y_Lat = y, Pop = pop_all)%>%
  dplyr::select(X_Lon, Y_Lat, Year, Pop)

PM_SHAP_wt_table <- PM_source_extend_China%>%
  left_join(Pop_China)%>%
  filter(Year%in%c(2001,2005,2010,2015,2017,2019))%>%
  group_by(Year)%>%
  summarise_at(c(paste0("PM25_",sector),"PM25"), ~ round(weighted.mean(.x,Pop),2))%>%
  mutate(anthro_sum = PM25_agriculture+PM25_residential+PM25_power+PM25_industry+PM25_transportation,
         Source = "SHAP")%>%
  ungroup()%>%
  as.data.frame()

PM_GBD_wt_table <- PM_source_GBD_2017%>%
  left_join(Pop_China)%>%
  filter(Year%in%c(2017))%>%
  group_by(Year)%>%
  summarise_at(c(paste0("PM25_",sector),"PM25"), ~ round(weighted.mean(.x,Pop),2))%>%
  mutate(anthro_sum = PM25_agriculture+PM25_residential+PM25_power+PM25_industry+PM25_transportation,
         Source = "GBD")%>%
  ungroup()%>%
  as.data.frame()

PM_GBD_wt_table%>%
  full_join(PM_SHAP_wt_table)%>%
  fwrite(paste0("result/Source_change/China/PM25_source_population-weighted_summary.csv"))

##  Spatial deviations between different simulations
##======Difference with GBD simulation======
PM_source_GBD_diff <- PM_source_deviation_China%>%
  filter(Year == 2017)%>%
  mutate(PM25_diff = PM25_SHAP - PM25_GBD,
         PM25_diff_percent = PM25_diff/PM25_SHAP*100,
         PM25_ratio_diff = PM25_ratio_diff*100)

PM_source_GBD_abs_diff <- PM_source_GBD_diff%>%
  dplyr::select(X_Lon, Y_Lat, Year, Sector, PM25_diff)%>%
  mutate(PM25_diff = abs(PM25_diff))%>%
  pivot_wider(names_from = "Sector", values_from = "PM25_diff")

PM_source_ratio_GBD_abs_diff <- PM_source_GBD_diff%>%
  dplyr::select(X_Lon, Y_Lat, Year, Sector, PM25_ratio_diff)%>%
  mutate(PM25_ratio_diff = abs(PM25_ratio_diff))%>%
  pivot_wider(names_from = "Sector", values_from = "PM25_ratio_diff")

text_table <- sector%>%
  map_dfr(function(i){
    df_sel <- PM_source_deviation_table%>%
      filter(Year == 2017, Sector == paste0("PM25_",i))
    data.table(label = c(paste0("~R^{2}==",formatC(df_sel$R2,digits = 2,format = "f",flag = "#")%>%str_c()%>%deparse()),
                         paste0("~MAE==",formatC(df_sel$error,digits = 2,format = "f",flag = "#")%>%str_c()%>%deparse())),
               vjust = c(.5,2), hjust = c(0,0))%>%
      mutate(Sector = str_to_title(i))
  })

# breaks <- c(-Inf,-5,-3,-1,-0.5,0.5,1,3,5,+Inf)
breaks <- c(-Inf,0.5,1,2,3,4,5,+Inf)
labels <- label_gen(breaks)
colors <- rev(colorspace::sequential_hcl(length(labels), palette = "Purple-Blue"))
# colors[5] <- "white"
PM_samod_GBD_abs_diff_plot <- PM_source_GBD_abs_diff%>%
  Abs_source_plot(2017, breaks, colors, "China", text_table = text_table,
                  title = expression(paste("Absolute differences of PM"[2.5]," (μg/m"^"3",")")))
ggplot_oupt(PM_samod_GBD_abs_diff_plot, filename = paste0("result/Source_change/China/PM25_source_contributions_absolute_diff_GBD_2017.png"),
            export_type = "Cairo", plot_type = "local", whr = c(3500,1800,155))


breaks <- c(-Inf,1,2,3,5,10,15,+Inf)
labels <- label_gen(breaks)
colors <- rev(colorspace::sequential_hcl(length(labels), palette = "Purple-Blue"))
PM_samod_ratio_GBD_abs_diff_plot <- PM_source_ratio_GBD_abs_diff%>%
  Abs_source_plot(2017, breaks, colors, "China",
                  title = expression(paste("Absolute differences of contribution ratios (%)")))

ggplot_oupt(PM_samod_ratio_GBD_abs_diff_plot, filename = paste0("result/Source_change/China/PM25_source_ratio_contributions_absolute_diff_GBD_2017.png"),
            export_type = "Cairo", plot_type = "local", whr = c(3500,1800,155))

##======Difference with GEOS-Chem simulation======
## Annual contributions
base_year <- 2019
PM_GC_annual <- 1:12%>%
  purrr::map_dfr(function(month){
    fread(paste0(dat_dir,"CTM/geoschem/GC_output_China/source_impacts/GC_baseline_PM25_",dsname,"_0.1g_",base_year,
                 formatC(month,flag = "0",width = 2),".csv"))%>%
      dplyr::select(x,y,Month,PM25)%>%rename(PM25_bs = PM25)
  })%>%
  group_by(x,y)%>%
  summarise(PM25_GC = mean(PM25_bs))%>%
  ungroup()%>%
  rename(X_Lon = x, Y_Lat = y)


PM_GC_diff <- PM_source_extend_China%>%
  filter(Year == base_year)%>%
  left_join(PM_GC_annual)%>%
  mutate(PM25_diff = PM25 - PM25_GC)

breaks <- c(-Inf,5,10,15,20,25,30,+Inf)
labels <- label_gen(breaks)
colors <- rev(colorspace::sequential_hcl(length(labels), palette = "Purple-Blue"))
{
  PM_GC_diff_proj <- PM_GC_diff%>%
    dplyr::select(X_Lon, Y_Lat, PM25_diff)%>%
    mutate(PM25_diff = abs(PM25_diff))%>%
    dfproj(toCRS = CRS_China, method = "average")%>%
    rename(X_Lon = x, Y_Lat = y)%>%
    mutate_at(c("PM25_diff"), ~ cut(.,breaks,labels))
  
  PM_GC_diff_plot <- PM_GC_diff_proj%>%
    grid_plot(x = "X_Lon", y = "Y_Lat", value = "PM25_diff")%>%
    grid_plot_set(area = "China", colors = colors,
                  title = expression(paste("Absolute differences of PM"[2.5]," (μg/m"^"3",")")),
                  panel_setting = "single")+
    theme(panel.grid.major = element_blank())

  PM_GC_diff_nanhai_plot <- PM_GC_diff_plot%>%
    nanhai_zoom_plot(border.line.size = .8)
  
  PM_GC_diff_overall_plot <- (PM_GC_diff_plot+PM_GC_diff_nanhai_plot)
  ggplot_oupt(PM_GC_diff_overall_plot, filename = paste0("result/Source_change/China/PM25_absolute_diff_GC_",base_year,".png"),
              export_type = "Cairo", plot_type = "single")
}

## Source-specific PM2.5 difference
PM_source_GC_diff <- PM_source_deviation_China%>%
  filter(Year == base_year)%>%
  mutate(PM25_diff = PM25_SHAP - PM25_GC,
         PM25_diff_percent = PM25_diff/PM25_SHAP*100)

PM_source_GC_abs_diff <- PM_source_GC_diff%>%
  dplyr::select(X_Lon, Y_Lat, Year, Sector, PM25_diff)%>%
  mutate(PM25_diff = abs(PM25_diff))%>%
  pivot_wider(names_from = "Sector", values_from = "PM25_diff")

PM_source_GC_rel_diff <- PM_source_GC_diff%>%
  dplyr::select(X_Lon, Y_Lat, Year, Sector, PM25_diff_percent)%>%
  pivot_wider(names_from = "Sector", values_from = "PM25_diff_percent")

breaks <- c(-Inf,0.5,1,2,3,4,5,+Inf)
labels <- label_gen(breaks)
colors <- rev(colorspace::sequential_hcl(length(labels), palette = "Purple-Blue"))
PM_samod_GC_abs_diff_plot <- PM_source_GC_abs_diff%>%
  Abs_source_plot(base_year, breaks, colors, "China",
                  title = expression(paste("Absolute differences of PM"[2.5]," (μg/m"^"3",")")))
ggplot_oupt(PM_samod_GC_abs_diff_plot, filename = paste0("result/Source_change/China/PM25_source_contributions_absolute_diff_GC_",base_year,".png"),
            export_type = "Cairo", plot_type = "local", whr = c(3500,1800,155))

##======Overall emission change======
meic_species <- c("BC","OC","NOx","SO2","CO","NH3","NMVOC","PM25")
species_sel <- c("BC","OC","PM25")

Abs_emi_plot <- function(data, year = NULL, emi_name, breaks = c(-Inf,1,3,5,10,15,20,+Inf),
                         colors, area = "China", text_table = NULL, facet_row = 2,
                         title = expression(paste("PM"[2.5]," concentration (μg/m"^"3",")"))){
  labels <- label_gen(breaks)
  if(!is.null(year)){
    source_apportion <- data%>%
      filter(Year == year)%>%
      dplyr::select(X_Lon,Y_Lat,all_of(emi_name))
  }else{
    source_apportion <- data%>%
      dplyr::select(X_Lon,Y_Lat,all_of(emi_name))
  }
  
  source_apportion_proj <- source_apportion%>%
    dfproj(toCRS = CRS_China, method = "average")%>%
    rename(X_Lon = x, Y_Lat = y)%>%
    pivot_longer(all_of(emi_name), names_to = "Sector")
  
  samod <- source_apportion_proj%>%
    mutate_at(c("value"), ~ cut(.,breaks,labels))
  
  samod_plot <- samod%>%
    grid_facet_plot(x = "X_Lon", y = "Y_Lat", value = "value",
                    facet = "Sector", facet_row = facet_row, facet_text = 30)%>%
    grid_plot_set(area = area, colors = colors,
                  title = title, 
                  size_list = c(0.1, 1, 25, 28, 30),
                  panel_setting = "facet")+
    scale_x_continuous(breaks = c(80,100,120))+
    coord_sf(xlim = c(-2825586,2406965))+
    theme(panel.grid.major = element_blank())
  
  if(!is.null(text_table)){
    samod_plot <- samod_plot%>%
      grid_plot_add_text(text_table = text_table, cex = 8)
  }
  
  samod_nanhai_plot <- samod%>%
    split(f = .$Sector)%>%
    purrr::map(~ (.x%>%
                    grid_plot(x = "X_Lon", y = "Y_Lat", value = "value")%>%
                    grid_plot_set(area = area, colors = colors,
                                  title = NULL,
                                  panel_setting = "facet")+
                    theme(panel.grid.major = element_blank())+
                    labs(title = NULL))%>%
                 nanhai_zoom_plot(facet_df = data.frame(Sector = unique(.x$Sector)),
                                  border.line.size = .8)
    )
  
  samod_overall_plot <- (samod_plot+samod_nanhai_plot)
  return(samod_overall_plot)
}

year_seq <- c(2001,2019)
Grid_rejoin_China <- fread(paste0("result/Inequality/Economic_imapct_grid_dataset_China.csv"))

{
  emi_name <- c("agriculture_NH3","industry_SO2","industry_BC","industry_OC","industry_PM25",
                "power_SO2","power_PM25","residential_SO2","residential_BC","residential_OC","residential_PM25",
                "transportation_NOx", "transportation_BC","transportation_OC","transportation_PM25")
  lab_name <- c("Agriculture NH3","Industry SO2","Industry BC","Industry OC","Industry PM25",
                "Power SO2","Power PM25","Residential SO2","Residential BC","Residential OC","Residential PM25",
                "Transportation NOx","Transportation BC","Transportation OC","Transportation PM25")
  
  emi_change <- Grid_rejoin_China%>%
    filter(Year%in%year_seq)%>%
    arrange(GridID, Year)%>%
    mutate_at(all_of(emi_name), ~.x - dplyr::lag(.x,1))%>%
    filter(Year == year_seq[2])%>%
    dplyr::select(X_Lon,Y_Lat,all_of(emi_name))%>%
    mutate_at(all_of(emi_name),~.x/1e-13)%>%
    rename_at(all_of(emi_name),~paste0(toupper(substr(.x,1,1)), substr(.x,2,nchar(.x)))%>%
                gsub("_"," ",.))

  breaks <- c(-Inf,-3e-11,-3e-12,-1e-12,-1e-13,1e-13,1e-12,3e-12,3e-11,+Inf)
  breaks_mod <- breaks/1e-13
  labels <- label_gen(breaks_mod)

  colors <- colorspace::diverge_hcl(length(labels), palette = "Blue-Red")
  emi_change_plot <- Abs_emi_plot(emi_change, emi_name = lab_name, breaks = breaks_mod,
                                  colors = colors, facet_row = 5,
                                  title = expression(paste("Change in emission fluxes (10"^"-13"," kg/m"^"2","/s)")))
  
  CairoPNG(paste0("result/Source_change/China/Emission_change_China_",year_seq[1],"-",year_seq[2],".png"),
           width = 3000, height = 3000, res = 140)
  print(emi_change_plot)
  dev.off()
}
