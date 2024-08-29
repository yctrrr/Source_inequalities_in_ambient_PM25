setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
library(tidyverse)
library(data.table)
library(terra)
library(exactextractr)
library(Cairo)
library(DescTools)
library(ggplot2)
library(patchwork)
library(qdapRegex)

source('D:/shaoyanchuan/codebook/function/Grid_visualization.R')
source('D:/shaoyanchuan/codebook/function/Health_calcu_core.R')
source('D:/shaoyanchuan/codebook/function/Time_series_analysis.R')
source("script/Inequality_cal_function.R")

load_shp()
dat_dir <- "D:/shaoyanchuan/data/"
# speciesvar_ceds <- c("NOx","SO2","NMVOC","NH3","CO","BC","OC")
speciesvar_meic <- c("NOx","SO2","NMVOC","NH3","BC","OC","PM25")
sector <- c("agriculture","industry","power","residential","transportation")
metvar <- c("u10","v10","t2m","blh","tcc","sp","e","lsp")

Prov_info <- fread(paste0("data/input/Province_information.csv"))%>%
  dplyr::select(ProvID,engname,Subarea)

Prov_sf_China <- polygon_sf_China%>%
  dplyr::select(1)%>%
  setnames(c("ProvID","geometry"))%>%
  left_join(Prov_info)

Prov_sf_CEC <- Prov_sf_China%>%
  filter(Subarea %in% c("N","E"))

Prov_centroids <- st_centroid(Prov_sf_China)%>%
  distinct(ProvID, .keep_all = TRUE)

Prov_centroids <- Prov_centroids%>%
  st_coordinates()%>%
  as.data.frame()%>%
  cbind(st_drop_geometry(Prov_centroids))

Prov_centroids[which(Prov_centroids$engname == "Hebei"),2] <- 4554682
Prov_centroids[which(Prov_centroids$engname == "Tianjin"),2] <- 4250000
Prov_centroids[which(Prov_centroids$engname == "Shanghai"),1] <- 1750000
Prov_centroids[which(Prov_centroids$engname == "Sichuan"),2] <- 3380000
Prov_centroids[which(Prov_centroids$engname == "Guangdong"),1] <- 895000
Prov_centroids[which(Prov_centroids$engname == "Guangdong"),2] <- 2340000
Prov_centroids[which(Prov_centroids$engname == "Shaanxi"),2] <- 3600000

City_sf_China <- readRDS(paste0(dat_dir,"Shp/city_China/city_shp.Rds"))%>%
  st_transform(crs(Prov_sf_China))%>%
  dplyr::select(1)%>%
  setnames(c("ProvID","geometry"))%>%
  left_join(Prov_info)

City_sf_CEC <- City_sf_China%>%
  filter(Subarea %in% c("N","E"))

County_sf_China <- readRDS(paste0(dat_dir,"Shp/county_China/county_shp.Rds"))%>%
  dplyr::select(-pop_2000, -pop_2010, -pop_2020)%>%
  rename(CountyID = PAC)

County_sf_rejoin_China <- fread("result/Inequality/Economic_imapct_county_dataset_China.csv", encoding = 'UTF-8')%>%
  left_join(County_sf_China%>%dplyr::select(CountyID),.)%>%
  na.omit()
Grid_rejoin_China <- fread(paste0("result/Inequality/Economic_imapct_grid_dataset_China.csv"))%>%
  left_join(Prov_info)

Pop_China <- Grid_rejoin_China%>%
  dplyr::select(X_Lon, Y_Lat, GridID, Year, Pop)

Grid_divide <- fread(paste0("result/Inequality/Grid_division_2001-2019.csv"))
Grid_sr_divide <- fread(paste0("result/Inequality/Grid_division_CEC_area_2001-2019.csv"))


speciesvar_meic <- c("NOx","SO2","NMVOC","NH3","BC","OC","PM25")
# emivar_SHAP <- paste0(speciesvar_meic, "_", sector_sel,"_SHAP")
# EMI_look_up <- expand.grid(Sector = setdiff(sector,"agriculture"), Species = speciesvar_meic)%>%
#   full_join(data.frame(Sector = "agriculture", Species = "NH3"))%>%
#   mutate(EMI_SHAP_var = paste0(Species,"_",Sector,"_SHAP"),
#          EMI_var = paste0(Sector,"_",Species))
EMI_look_up <- expand.grid(Sector = sector, Species = speciesvar_meic)%>%
  mutate(EMI_SHAP_var = paste0(Species,"_",Sector,"_SHAP"),
         EMI_var = paste0(Sector,"_",Species))
names(EMI_SHAP_China)

EMI_SHAP_China <- fread("result/Interpretation/SHAP_China_TAP/PM_RRF_MEIC_SHAP_2001-2019_China.csv")%>%
  dplyr::select(GridID,Year,all_of(EMI_look_up[,"EMI_SHAP_var"]),all_of(sector))

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

regional_plot_set <- function(base_map){
  base_map+
    xlim(420000,1680000)+ylim(2550000,4740000)+
    theme(panel.grid.major = element_blank(),
          panel.border = element_rect(colour = "black",linewidth = 2.5),
          plot.title = element_text(size = rel(2.6), hjust = .5, vjust = .4),
          plot.margin = margin(0.5,0,0,1,"cm"),
          legend.title = element_text(size = rel(2.6), hjust = .5, vjust = 1),
          legend.text = element_text(size = rel(2.1)),
          axis.text = element_blank(),
          axis.ticks = element_blank())
}

breaks_set <- function(var, year){
  if(var == "GDP_USD"){
    if(year <= 2012){
      breaks <<- c(-Inf,0.2,0.4,0.6,0.8,1,1.2,1.4,+Inf)
    }else if (year == 2019|year == 2017|year == 2013){
      breaks <<- c(-Inf,2,3,4,5,10,15,20,+Inf)
    }
    label <<- label_gen(breaks)
    legend_colors <<- rev(RColorBrewer::brewer.pal(length(label),"RdYlBu"))
    lab_title <<- "GDP per capita"
  }else if(var == "PM25_residential"){
    if(year <= 2012){
      breaks <<- c(-Inf,3,5,10,15,20,25,30,+Inf)
    }else if (year >= 2013){
      breaks <<- c(-Inf,3,5,10,15,20,25,30,+Inf)
    }
    label <<- label_gen(breaks)
    # legend_colors <<- rev(colorspace::sequential_hcl(length(label), palette = "YlOrRd"))
    legend_colors <<- rev(RColorBrewer::brewer.pal(length(labels), "Spectral"))
    lab_title <<- expression(paste("Residential PM"[2.5]," concentration"))
  }else if(var == "PM25_transportation"){
    if(year <= 2012){
      breaks <<- c(-Inf,1,2,3,4,5,6,+Inf)
    }else if (year >= 2013){
      breaks <<- c(-Inf,1,2,3,4,5,6,+Inf)
    }
  }else if(var == "residential_BC"){
    breaks <<- c(-Inf,-1e-11,-5e-12,-3e-12,-2e-12,-1e-12,-5e-13,0,+Inf)
    label <<- label_gen(breaks)
    # legend_colors <<- rev(colorspace::sequential_hcl(length(label)-3, palette = "Viridis"))%>%
    #   c("white","#E54924","#88002D")
    legend_colors <<- colorspace::diverge_hcl(length(label), palette = "Blue-Red")
    lab_title <<- expression(paste("BC emission fluxes"))
  }else if(var == "residential_OC"){
    # breaks <<- c(-Inf,4e-13,8e-13,2e-12,5e-12,1e-11,2e-11,3e-11,+Inf)
    breaks <<- c(-Inf,8e-13,2e-12,5e-12,1e-11,2e-11,5e-11,1e-10,+Inf)
    label <<- label_gen(breaks)
    colors <<- rev(colorspace::sequential_hcl(length(label), palette = "OrRd"))
    # title <<- expression(paste("OC emission fluxes"))
    # breaks <<- c(-Inf,-1e-11,-5e-12,-3e-12,-2e-12,-1e-12,-5e-13,0,+Inf)
    # label <<- label_gen(breaks)
    # colorspace::hcl_palettes("sequential")
    # legend_colors <<- rev(colorspace::sequential_hcl(length(label)-3, palette = "Viridis"))%>%
    #   c("white","#E54924","#88002D")
    lab_title <<- expression(paste("OC emission fluxes"))
  }else if(var == "residential_SO2"){
    # breaks <<- c(-Inf,4e-13,8e-13,2e-12,5e-12,1e-11,2e-11,3e-11,+Inf)
    breaks <<- c(-Inf,8e-13,2e-12,5e-12,1e-11,2e-11,5e-11,1e-10,+Inf)
    label <<- label_gen(breaks)
    legend_colors <<- colorspace::diverge_hcl(length(label), palette = "Blue-Red")
    lab_title <<- expression(paste("SO2 emission fluxes"))
  }else if(var == "residential_NOx"){
    breaks <<- c(-Inf,4e-13,8e-13,2e-12,5e-12,1e-11,2e-11,3e-11,+Inf)
    label <<- label_gen(breaks)
    legend_colors <<- colorspace::diverge_hcl(length(label), palette = "Blue-Red")
    lab_title <<- expression(paste("NO"[x]," emission fluxes"))
  }else if(var == "BC_residential_SHAP"){

    if(year >= 2013){
      breaks <<- c(-Inf,-3,-2,-1,-0.5,-0.1,0,+Inf)
      label <<- label_gen(breaks)
      legend_colors <<- colorspace::sequential_hcl(length(label)-1, palette = "agGrnYl")%>%
        c("#feb24c")
    }else{
      breaks <<- c(-Inf,-3,-2,-1,-0.5,-0.1,0.1,0.5,2,+Inf)
      label <<- label_gen(breaks)
      legend_colors <<- colorspace::sequential_hcl(length(label)-4, palette = "agGrnYl")%>%
        c("white","#ffffb2","#feb24c","#d95f0e")
    }
    
    lab_title <<- expression(paste("BC-related PM"[2.5]))
  }else if(var == "OC_residential_SHAP"){

    if(year >= 2013){
      breaks <<- c(-Inf,-3,-2,-1,-0.5,-0.1,0,+Inf)
      label <<- label_gen(breaks)
      legend_colors <<- colorspace::sequential_hcl(length(label)-1, palette = "agGrnYl")%>%
        c("#feb24c")
    }else{
      breaks <<- c(-Inf,-3,-2,-1,-0.5,-0.1,0.1,0.5,2,+Inf)
      label <<- label_gen(breaks)
      legend_colors <<- colorspace::sequential_hcl(length(label)-4, palette = "agGrnYl")%>%
        c("white","#ffffb2","#feb24c","#d95f0e")
    }
    
    lab_title <<- expression(paste("OC-related PM"[2.5]))
  }
}

##======Change in population-weighted SHAP======
## We will separately calculate the sector-species SHAP contributions and related emission change
sector_sel <- "residential"
EMI_look_up
EMI_SHAP_by_year <- sector%>%
  map_dfr(function(sector_sel){
    join_divide <- Grid_divide
    sr_join_divide <- Grid_sr_divide
    
    EMI_SHAP_var_sel <- EMI_look_up%>%
      filter(Sector == sector_sel)%>%
      .[["EMI_SHAP_var"]]
    
    EMI_SHAP_sel <- EMI_SHAP_China%>%
      dplyr::select(GridID,Year,all_of(EMI_SHAP_var_sel))%>%
      rename_at(all_of(EMI_SHAP_var_sel), ~gsub(paste0(sector_sel,"|_|SHAP"),"",.x))%>%
      mutate(Sector = sector_sel, .after = "Year")%>%
      left_join(Pop_China)%>%
      na.omit()
    EMI_SHAP_sel%>%summary
    ## population-weighted SHAP may differ from the values calculated by population-weighted PM2.5
    # EMI_SHAP_all <- EMI_SHAP_sel%>%
    #   group_by(Year, Sector)%>%
    #   summarise_at(all_of(speciesvar_meic),~weighted.mean(.x,Pop))%>%
    #   ungroup()%>%
    #   mutate(Area = "Overall", .after = "Sector",
    #          SHAP := !!rlang::parse_expr(paste(speciesvar_meic, collapse = "+")))%>%
    #   suppressMessages()
    
    
    EMI_SHAP_all <- speciesvar_meic%>%
      map_dfr(~ EMI_SHAP_sel%>%
                left_join(join_divide)%>%
                na.omit()%>%
                local_popw_summarise(group_main = "GDPQ", group_add = "Year",
                                     var = .x, prefix = FALSE)%>%
                mutate(Species = .x, .after = "Year")%>%
                suppressMessages())%>%
      mutate(Sector = sector_sel, Area = "Overall", .after = "Species")

    EMI_SHAP_sr <- speciesvar_meic%>%
      map_dfr(~ EMI_SHAP_sel%>%
                left_join(sr_join_divide)%>%
                na.omit()%>%
                local_popw_summarise(group_main = "GDPQ", group_add = "Year",
                                     var = .x, prefix = FALSE)%>%
                mutate(Species = .x, .after = "Year")%>%
                suppressMessages())%>%
      mutate(Sector = sector_sel, Area = "NE", .after = "Species")
    
    ## population-weighted SHAP for overall populations and populations living in central eastern China
    full_join(EMI_SHAP_all, EMI_SHAP_sr)
  })

fwrite(EMI_SHAP_by_year, paste0("result/Inequality/Emission_analysis/Emission_SHAP_change_by_year.csv"))

breaks <- c(-Inf,-16,-12,-8,-4,-2,-1,1,2,+Inf)
labels <- label_gen(breaks)
colors <- colorspace::sequential_hcl(length(labels)-3, palette = "agGrnYl")%>%
  c("#f0f0f0","#fed98e","#fd8d3c")
sector_sel <- "residential"
sector_sel <- "industry"
EMI_SHAP_var_sel <- EMI_look_up%>%
  filter(Sector == sector_sel)%>%
  .[["EMI_SHAP_var"]]
EMI_SHAP_var_sel
EMI_SHAP_China%>%
  left_join(Grid_sr_divide,)%>%
  na.omit%>%
  filter(Year == 2013)%>%
  dplyr::select(X_Lon,Y_Lat,Year,all_of(EMI_SHAP_var_sel))%>%
  rename_at(all_of(EMI_SHAP_var_sel), ~gsub(paste0(sector_sel,"|_|SHAP"),"",.x))%>%
  mutate(TPM25 = NOx+SO2+NMVOC+NH3+BC+OC+PM25)%>%
  dplyr::select(X_Lon,Y_Lat,Year,SO2,TPM25)%>%left_join(SHAP_var_df)%>%
  dplyr::select(-GridID,-Year)%>%
  rast%>%plot(col = colors, range = c(-16,10))

## Average emissions across the economic regions
## weighted by cell areas
cell_area <- Grid_divide%>%
  dplyr::select(X_Lon,Y_Lat)%>%
  distinct()%>%mutate(ID = 1)%>%
  rast(crs = "WGS84")%>%
  cellSize(unit = "m")%>%
  as.data.frame(xy = TRUE)%>%
  setnames(c("X_Lon","Y_Lat","cell_area"))%>%
  mutate_at(c("X_Lon","Y_Lat"), ~round(.x,2))

## unit:kg/m2/s
EMI_by_year <- sector%>%
  map_dfr(function(sector_sel){
    join_divide <- Grid_divide%>%
      left_join(cell_area)%>%
      na.omit()
    sr_join_divide <- Grid_sr_divide%>%
      left_join(cell_area)%>%
      na.omit()
    
    EMI_var_sel <- EMI_look_up%>%
      filter(Sector == sector_sel)%>%
      .[["EMI_var"]]
    speciesvar <- speciesvar_meic
    
    EMI_sel <- Grid_rejoin_China%>%
      dplyr::select(X_Lon,Y_Lat,GridID,Year,all_of(EMI_var_sel))%>%
      rename_at(all_of(EMI_var_sel), ~ gsub(paste0(sector_sel,"|_"),"",.x))%>%
      mutate(Sector = sector_sel, .after = "Year")%>%
      left_join(Pop_China)%>%
      na.omit()
    
    EMI_all <- speciesvar%>%
      map_dfr(~ EMI_sel%>%
                left_join(join_divide)%>%
                na.omit()%>%
                local_popw_summarise(group_main = "GDPQ", group_add = "Year",
                                     var = .x, pop_name = "cell_area", prefix = FALSE)%>%
                mutate(Species = .x, .after = "Year")%>%
                suppressMessages())%>%
      mutate(Sector = sector_sel, Area = "Overall", .after = "Species")
    
    EMI_sr <- speciesvar%>%
      map_dfr(~ EMI_sel%>%
                left_join(sr_join_divide)%>%
                na.omit()%>%
                local_popw_summarise(group_main = "GDPQ", group_add = "Year",
                                     var = .x, pop_name = "cell_area", prefix = FALSE)%>%
                mutate(Species = .x, .after = "Year")%>%
                suppressMessages())%>%
      mutate(Sector = sector_sel, Area = "NE", .after = "Species")
    
    full_join(EMI_all, EMI_sr)
  })

fwrite(EMI_by_year, paste0("result/Inequality/Emission_analysis/Emission_change_by_year.csv"))

EMI_SHAP_sel <- fread(paste0("result/Inequality/Emission_analysis/Emission_SHAP_change_by_year.csv"))%>%
  dplyr::select(GDPQ,Year,Species,Sector,Area,Popw,Pop)%>%
  filter(Year%in%c(2001,2013,2019))%>%
  pivot_wider(names_from = "Species", values_from = "Popw", names_prefix = "SHAP_")%>%
  mutate_at(paste0("SHAP_",speciesvar_meic), ~.x - dplyr::lag(.x,1))%>%
  group_by(GDPQ,Sector)%>%
  mutate_at(paste0("SHAP_",speciesvar_meic), ~ifelse(Year == 2001, .x[Year == 2013]+.x[Year == 2019],.x))%>%
  ungroup()%>%
  mutate(SHAP := !!rlang::parse_expr(paste(paste0("SHAP_",speciesvar_meic), collapse = "+")))
  # mutate_at(paste0("SHAP_",speciesvar_meic), ~ifelse(Year==2001,NA,.x))%>%
  # filter(Year!=2001)

EMI_sel_tmp <- fread(paste0("result/Inequality/Emission_analysis/Emission_change_by_year.csv"))%>%
  dplyr::select(GDPQ,Year,Species,Sector,Area,Popw)%>%
  filter(Year%in%c(2001,2013,2019))%>%
  pivot_wider(names_from = "Species", values_from = "Popw")%>%
  arrange(Area,Sector,GDPQ,Year)

EMI_sel <- EMI_sel_tmp%>%
  group_by(Area,Sector,GDPQ)%>%
  mutate_at(all_of(speciesvar_meic), ~ .x - dplyr::lag(.x,1))%>%
  rename_at(all_of(speciesvar_meic), ~ paste0(.x,"_change"))%>%
  left_join(EMI_sel_tmp,.)
  # mutate_at(paste0(speciesvar_meic), ~.x - dplyr::lag(.x,1))%>%
  # mutate_at(paste0(speciesvar_meic), ~ifelse(Year==2001,NA,.x))%>%
  
left_join(EMI_sel, EMI_SHAP_sel)%>%
  fwrite(paste0("result/Inequality/Emission_analysis/Emission_integrated_change_summary.csv"))

##======Source contributions and emissions in 2001 and 2019======
look_up <- expand.grid(Region = list(c("N","E")),
                       Year = list(c(2001,2013),c(2013,2019)))
# EMI_SHAP_China <- fread("result/Interpretation/SHAP_China_TAP/PM_RRF_MEIC_SHAP_2001-2019_China.csv")%>%
#   dplyr::select(GridID,Year,all_of(EMI_look_up[,"EMI_SHAP_var"]))
var = "BC_residential_SHAP"
var = "residential_BC"
var = "GDPQ"
# summary(Source_exposure_region)

Grid_rejoin_China%>%
  dplyr::select(-Pop, -GDPQ)%>%
  left_join(Grid_sr_divide,.)%>%
  filter(Subarea %in% c("N","E"))%>%
  filter(Year == 2019, GDPQ == 5)%>%
  dplyr::select(X_Lon,Y_Lat,PM25_residential, residential_OC)%>%rast%>%plot
s = "industry"
poly_sf_sel <- Prov_sf_CEC
# poly_sf_sel <- City_sf_CEC
year = 2019
for(s in c("residential","industry")){
  if(s == "residential"){
    var_vect <- c(paste0("PM25_",s),paste0(s,"_BC"),paste0(s,"_OC"),paste0(s,"_PM25"))
  }else if(s == "industry"){
    var_vect <- c(paste0("PM25_",s),paste0(s,"_SO2"),paste0(s,"_BC"),paste0(s,"_OC"),paste0(s,"_PM25"))
  }
  emi_vect <- setdiff(var_vect, paste0("PM25_",s))
  
  var_plot_set <- var_vect%>%
    # c("GDP_USD","PM25_residential","residential_BC","residential_SO2")%>%
    purrr::map(function(var){
      plot_data <- Grid_rejoin_China%>%
        dplyr::select(-Pop, -GDPQ)%>%
        left_join(Grid_sr_divide,.)%>%
        filter(Year == year, Subarea %in% c("N","E"),
               # ProvID == 410000,
               )%>%
        mutate_at(all_of(emi_vect), ~.x/1e-13)%>%
        dplyr::select(X_Lon,Y_Lat,all_of(var))%>%
        dfproj(toCRS = CRS_China, method = "average")
      
      if(var == paste0("PM25_",s)){
        if(year <= 2012){
          breaks <- c(-Inf,3,6,9,12,15,18,21,+Inf)
        }else if (year >= 2013){
          breaks <- c(-Inf,3,6,9,12,15,18,21,+Inf)
        }
        label <- label_gen(breaks)
        legend_colors <- rev(RColorBrewer::brewer.pal(length(label), "Spectral"))
        if(s == "residential"){
          lab_title <- expression(paste("Residential source"))
        }else if(s == "industry"){
          lab_title <- expression(paste("Industrial source"))
        }
        
        legend_title <- expression(paste("PM"[2.5]," concentration (μg/m"^"3",")"))
      }else{
        breaks <- c(-Inf,4e-13,8e-13,2e-12,5e-12,1e-11,3e-11,6e-11,+Inf)
        breaks <- breaks/1e-13
        label <- label_gen(breaks_mod)
        legend_colors <- rev(colorspace::sequential_hcl(length(label), palette = "YLOrRd"))
        
        if(var == paste0(s,"_PM25")){
          lab_title <- expression(paste("Primary PM"[2.5]))
        }else if(var == paste0(s,"_SO2")){
          lab_title <- expression(paste("SO"[2]))
        }else{
          lab_title <- gsub(paste0(s,"_"),"",var)
        }
        legend_title <- expression(paste("Emission fluxes (10"^"-13"," kg/m"^"2","/s)"))
      }
      
      plot_data <- plot_data%>%
        mutate_at(var, ~ cut(.x,breaks,label))
      
      base_plot <- plot_data%>%
        grid_plot(x = "x", y = "y", value = var)%>%
        grid_plot_set(area = "no", value_type = "discrete", 
                      colors = legend_colors, title = legend_title,
                      panel_setting = "single")%>%
        regional_plot_set()+
        geom_sf(data = poly_sf_sel, fill = "white", alpha = 0.1, linewidth = 0.4)+
        labs(title = lab_title)
      base_plot
    })
  
  var_plot <- var_plot_set[[1]]+
    theme(plot.margin = margin(r = 5))

  for(i in 2:(length(var_vect)-1)){
    var_plot <- var_plot+
      var_plot_set[[i]]+
      theme(legend.position = "none")
  }
  var_plot <- var_plot+
    var_plot_set[[length(var_vect)]]+
    plot_layout(nrow = 1)
  
  CairoPNG(paste0("result/Inequality/Emission_analysis/",str_to_title(s),"_exposure_spatial_distribution_CEC_",year,".png"),
           width = 5200, height = 1600, res = 140)
  print(var_plot)
  dev.off()
}

# rast(Source_exposure_region)%>%plot
# plot(Source_exposure_region$GDP_USD,Source_exposure_region$PM25_residential)


##======Regional emission change======
GDPQ_group <- 1:3
year_seq <- c(2001,2019)
s <- "residential"
poly_sf_sel <- Prov_sf_CEC
for(s in c("residential","industry")){
  if(s == "residential"){
    species <- c("BC","OC","PM25")
    # species <- c("SO2","NMVOC","CO")
    breaks <- c(-Inf,-3e-11,-1e-11,-3e-12,-1e-12,1e-12,2e-12,3e-12,+Inf)
    
  }else if(s == "industry"){
    species <- c("SO2","NMVOC","PM25")
    breaks <- c(-Inf,-1e-10,-1e-11,-5e-12,-1e-12,1e-12,1e-11,1e-10,+Inf)
  }
  
  emivar <- paste0(s,"_",species)
  Source_var_plot_set <- emivar%>%
    # c("GDP_USD","PM25_residential","residential_BC","residential_SO2")%>%
    purrr::map(function(var){
      Source_exposure_region <- Grid_rejoin_China%>%
        dplyr::select(-Pop, -GDPQ)%>%
        left_join(Grid_sr_divide,.)%>%
        filter(Year %in% year_seq,
               Subarea %in% c("N","E")
               )%>%
        dplyr::select(X_Lon,Y_Lat,Year,all_of(var))%>%
        pivot_wider(names_from = 3, values_from = 4)%>%
        mutate(change := !!rlang::sym(paste0(year_seq[2]))-!!rlang::sym(paste0(year_seq[1])), .keep = 'unused')%>%
        dfproj(toCRS = CRS_China, method = "average")
      
      # mutate(GDPQR = ifelse(GDPQ >=4 ,"Rich","Poor"))
      
      if(var == paste0(s,"_PM25")){
        # breaks <- c(-Inf,-1e-11,-5e-12,-2e-12,-1e-12,1e-12,2e-12,+Inf)
        lab_title <- expression(paste("Primary PM"[2.5]))
        legend_title <- ""
      }else{
        
        lab_title <- gsub(paste0(s,"_"),"",var)
        legend_title <- expression(paste("Emission fluxes (kg/m"^"2","/s)"))
      }

      label <- label_gen(breaks)
      legend_colors <- colorspace::diverge_hcl(length(label), palette = "Blue-Red")
      
      poly_sf_sel <- Prov_sf_China%>%
        filter(Subarea %in% c("N","E"))
      
      base_plot <- Source_exposure_region%>%
        mutate_at("change", ~ cut(.x,breaks,label))%>%
        grid_plot(x = "x", y = "y", value = "change")%>%
        grid_plot_set(area = "no",value_type = "discrete", 
                      colors = legend_colors, 
                      title = legend_title,
                      panel_setting = "single")%>%
        regional_plot_set()+
        geom_sf(data = poly_sf_sel, fill = "white", alpha = 0.1, linewidth = 0.4)+
        labs(title = lab_title)
      
      base_plot
    })
  
  Source_var_plot <- Source_var_plot_set[[1]]+
    theme(legend.position = "none",
          # plot.margin = margin(0,10,0,0)
    )+
    Source_var_plot_set[[2]]+
    theme(legend.position = "none",
          # plot.margin = margin(0,10,0,0)
    )+
    Source_var_plot_set[[3]]+
    plot_layout(nrow = 1)
  
  CairoPNG(paste0("result/Inequality/Emission_analysis/Emission_change_by_region_",s,"_sector_CEC_",year_seq[1],"-",year_seq[2],".png"),
           width = 2800, height = 1300, res = 140)
  print(Source_var_plot)
  dev.off()
}


##======Regional emission to PM2.5 change======
GDPQ_group <- 1:3
year_seq <- c(2001,2019)
s <- "residential"
poly_sf_sel <- Prov_sf_CEC
for(s in c("residential","industry")){
  if(s == "residential"){
    species <- c("BC","OC","PM25")
    # species <- c("SO2","NMVOC","CO")
    breaks <-  c(-Inf,-16,-12,-8,-6,-4,-2,-1,1,+Inf)
    label <- label_gen(breaks)
    legend_colors <- colorspace::sequential_hcl(length(label)-2, palette = "agGrnYl")%>%
      c("#f0f0f0","#fed98e")
  }else if(s == "industry"){
    species <- c("SO2","NMVOC","PM25")
    # breaks <- c(-Inf,-16,-12,-4,-2,-1,0,1,2, 4+Inf)
    breaks <-  c(-Inf,-16,-12,-8,-6,-4,-2,-1,1,+Inf)
    label <- label_gen(breaks)
    legend_colors <- colorspace::sequential_hcl(length(label)-2, palette = "agGrnYl")%>%
      c("#f0f0f0","#fed98e")
  }

  legend_title <- expression(paste("Change in PM"[2.5]," (μg/m"^"3",")"))
  SHAPvar <- paste0(species,"_",s,"_SHAP")
  
  SHAP_var_df <- EMI_SHAP_China%>%
    left_join(Grid_sr_divide,.)%>%
    filter(Year %in% year_seq,
           Subarea %in% c("N","E"))%>%
    arrange(GridID,Year)%>%
    dplyr::select(X_Lon,Y_Lat,GridID,Year,Year,all_of(SHAPvar))%>%
    group_by(GridID)%>%
    mutate_at(all_of(SHAPvar), ~.x-dplyr::lag(.x))%>%
    ungroup()%>%
    filter(Year == year_seq[2])%>%
    dplyr::select(-GridID,-Year)%>%
    dfproj(toCRS = CRS_China, method = "average")%>%
    pivot_longer(!c("x","y"))%>%
    mutate(name = gsub(paste0("_",s,"_SHAP"),"",name),
           name = factor(name, levels = species))
  
  SHAP_var_plot <- SHAP_var_df%>%
    mutate_at("value", ~ cut(.x,breaks,label))%>%
    grid_facet_plot(x = "x", y = "y", facet = "name", facet_text = 45, facet_row = 1, value = "value")%>%
    grid_plot_set(area = "no",value_type = "discrete", 
                  colors = legend_colors, 
                  title = legend_title,
                  panel_setting = "facet")%>%
    regional_plot_set()+
    geom_sf(data = poly_sf_sel, fill = "white", alpha = 0.1, linewidth = 0.4)
  
  CairoPNG(paste0("result/Inequality/Emission_analysis/Emission_change_to_PM2.5_by_region_",s,"_sector_CEC_",year_seq[1],"-",year_seq[2],".png"),
           width = 2800, height = 1300, res = 140)
  print(SHAP_var_plot)
  dev.off()
}

##======Regional sectoral PM2.5 change======
year_seq <- c(2001,2013,2019)
poly_sf_sel <- Prov_sf_CEC
s = "residential"
s = "industry"
for(s in c("residential","industry")){
  if(s == "residential"){
    breaks <-  c(-Inf,-24,-16,-8,-4,-2,-1,1,2,4,+Inf)
    label <- label_gen(breaks)
    legend_colors <- colorspace::sequential_hcl(length(label)-3, palette = "agGrnYl")%>%
      c("#f0f0f0","#fed98e","#fd8d3c")
  }else if(s == "industry"){
    # breaks <- c(-Inf,-16,-12,-4,-2,-1,0,1,2, 4+Inf)
    breaks <-  c(-Inf,-24,-16,-8,-4,-1,1,4,8,16,+Inf)
    label <- label_gen(breaks)
    legend_colors <- colorspace::sequential_hcl(length(label)-5, palette = "agGrnYl")%>%
      c("#f0f0f0","#fed98e","#fd8d3c","#fc4e2a","#b10026")
  }
  
  legend_title <- expression(paste("Change in PM"[2.5]," (μg/m"^"3",")"))
  SHAPvar <- s

  SHAP_var_df <- EMI_SHAP_China%>%
    left_join(Grid_sr_divide,.)%>%
    filter(Year %in% year_seq,
           Subarea %in% c("N","E"))%>%
    arrange(GridID,Year)%>%
    dplyr::select(X_Lon,Y_Lat,GridID,Year,all_of(SHAPvar))%>%
    group_by(GridID)%>%
    mutate_at(all_of(SHAPvar), ~.x-dplyr::lag(.x))%>%
    ungroup()%>%
    filter(Year != year_seq[1])
  
  SHAP_var_by_period <- 1:(length(year_seq)-1)%>%
    map_dfr(~ SHAP_var_df%>%
              filter(Year == year_seq[.x+1])%>%
              dplyr::select(-GridID,-Year)%>%
              dfproj(toCRS = CRS_China, method = "average")%>%
              rename(value = rlang::sym(SHAPvar))%>%
              mutate(Year = paste0(year_seq[.x],"-",year_seq[.x+1])))


  SHAP_var_plot <- SHAP_var_by_period%>%
    mutate_at("value", ~ cut(.x,breaks,label))%>%
    grid_facet_plot(x = "x", y = "y", facet = "Year", facet_text = 45, facet_row = 1, value = "value")%>%
    grid_plot_set(area = "no",value_type = "discrete", 
                  colors = legend_colors, 
                  title = legend_title,
                  panel_setting = "facet")%>%
    regional_plot_set()+
    theme(strip.text.x = element_text(
      size = 45, color = "black",margin = margin(t = 0.1, r = -.1, b = 0.3, l = -.1, "cm"), vjust = .5),
      legend.title = element_text(size = rel(2.3), hjust = .5, vjust = 1),
      legend.text = element_text(size = rel(1.8))
    )+
    geom_sf(data = poly_sf_sel, fill = "white", alpha = 0.1, linewidth = 0.4)

  CairoPNG(paste0("result/Inequality/Emission_analysis/PM2.5_change_by_region_",s,"_sector_CEC_",year_seq[1],"-",year_seq[length(year_seq)],".png"),
           width = 2000, height = 1300, res = 140)
  print(SHAP_var_plot)
  dev.off()
}
breaks <- c(-Inf,-16,-12,-8,-4,-2,-1,1,2,+Inf)
labels <- label_gen(breaks)
colors <- colorspace::sequential_hcl(length(labels)-3, palette = "agGrnYl")%>%
  c("#f0f0f0","#fed98e","#fd8d3c")

SHAP_var_df%>%filter(Year == 2013)%>%
  dplyr::select(X_Lon,Y_Lat,all_of(SHAPvar))%>%
  # rast()%>%plot
  rast%>%plot(col = colors, range = c(-16,10))

##======Time series of emission change======
colors <- colorspace::sequential_hcl(palette = "RdPu", n = 6)%>%rev()

for(s in c("residential","industry")){
  if(s == "residential"){
    species <- c("BC","OC","PM25")
    labels <- c("BC","OC",expression(paste("Primary PM"[2.5])))
  }else if(s == "industry"){
    species <- c("SO2","NMVOC","CO")
    labels <- c("SO2","NMVOC","CO")
  }
  ## emission change in central eastern China
  EMI_change_by_gdppc <- fread(paste0("result/Inequality/Emission_analysis/Emission_change_by_year.csv"))%>%
    dplyr::select(GDPQ,Year,Species,Sector,Area,Popw)%>%
    filter(Year%in%c(2001,2019), Area == "NE")%>%
    pivot_wider(names_from = "Species", values_from = "Popw")%>%
    arrange(Area,Sector,GDPQ,Year)%>%
    group_by(Area,Sector,GDPQ)%>%
    mutate_at(all_of(speciesvar_meic), ~ .x - dplyr::lag(.x,1))%>%
    ungroup()%>%
    pivot_longer(!c("GDPQ","Year","Sector","Area"), names_to = "Species")%>%
    filter(Year == 2019, Species %in% species, Sector == s)%>%
    # mutate(
    #        value = -value, 
    #        Species = paste0("Residential ",Species))%>%
    GDPQ_rename()

  # labels <- c("Residential BC","Residential OC",expression(paste("Residential PM"[2.5])))

  EMI_change_by_gdppc_plot <- box_series_plot(EMI_change_by_gdppc, x = "Species", y = "value", errorbar = FALSE,
                                            position = position_dodge(width = .9), fill = "GDPQ")+
    labs(x = NULL, y = expression(paste("Change in emissions (kg/m"^"2","/s)")))+
    scale_x_discrete(labels = labels)+
    scale_y_reverse()+
    scale_fill_manual(name = NULL, values = colors[2:6])+
    theme(legend.key.width = unit(2, "cm"),
          legend.key.height = unit(2, "cm"),  ##legend size
          legend.position = "top",
          legend.text = element_text(colour = "black", size = rel(2.8),margin = margin(l = 10,r = 10)),
          # legend.title = element_text(colour = "black", size = rel(3.2)),
          axis.text = element_text(colour = "black",size = rel(2.8)),
          axis.title = element_text(colour = "black",size = rel(3.2)),
          plot.margin = margin(0,20,0,10))
  
  
  CairoPNG(paste0("result/Inequality/Emission_analysis/Emission_change_by_gdppc_",s,"_sector_CEC_area_time_series.png"),
           width = 3000, height = 1500, res = 120)
  print(EMI_change_by_gdppc_plot)
  dev.off()
}


##======Time series of emission change to PM2.5======
colors <- colorspace::sequential_hcl(palette = "RdPu", n = 6)%>%rev()
speciesvar <- speciesvar_meic
for(s in c("residential","industry")){
  if(s == "residential"){
    species <- c("BC","OC","PM25")
    labels <- c("BC","OC",expression(paste("Primary PM"[2.5])))
  }else if(s == "industry"){
    species <- c("SO2","NMVOC","CO")
    labels <- c(expression(paste("SO"[2])),"NMVOC","CO")
  }
  ## emission change in central eastern China
  EMI_SHAP_change_by_gdppc <- fread(paste0("result/Inequality/Emission_analysis/Emission_SHAP_change_by_year.csv"))%>%
    dplyr::select(GDPQ,Year,Species,Sector,Area,Popw)%>%
    filter(Year%in%c(2001,2019), Area == "NE")%>%
    pivot_wider(names_from = "Species", values_from = "Popw")%>%
    arrange(Area,Sector,GDPQ,Year)%>%
    group_by(Area,Sector,GDPQ)%>%
    mutate_at(all_of(speciesvar), ~ .x - dplyr::lag(.x,1))%>%
    ungroup()%>%
    pivot_longer(!c("GDPQ","Year","Sector","Area"), names_to = "Species")%>%
    filter(Year == 2019, Species %in% species, Sector == s)%>%
    GDPQ_rename()

  EMI_SHAP_change_by_gdppc_plot <- box_series_plot(EMI_SHAP_change_by_gdppc, x = "Species", y = "value", errorbar = FALSE,
                                              position = position_dodge(width = .9), fill = "GDPQ")+
    labs(x = NULL, y = expression(paste("Change in PM"[2.5]," (μg/m"^"3",")")))+
    scale_x_discrete(labels = labels)+
    scale_y_reverse()+
    scale_fill_manual(name = NULL, values = colors[2:6])+
    theme(legend.key.width = unit(2, "cm"),
          legend.key.height = unit(2, "cm"),  ##legend size
          legend.position = "top",
          legend.text = element_text(colour = "black", size = rel(2.8),margin = margin(l = 10,r = 10)),
          # legend.title = element_text(colour = "black", size = rel(3.2)),
          axis.text = element_text(colour = "black",size = rel(2.8)),
          # axis.title = element_text(colour = "black",size = rel(3.2)),
          plot.margin = margin(0,20,0,10))
  
  
  CairoPNG(paste0("result/Inequality/Emission_analysis/Emission_change_to_PM2.5_by_gdppc_",s,"_sector_CEC_area_time_series.png"),
           width = 3000, height = 1500, res = 120)
  print(EMI_SHAP_change_by_gdppc_plot)
  dev.off()
}
##======Time series of sectoral PM2.5 change (by periods)======
GDP_source_summarise_sr_China <- 2001:2019%>%
  map_dfr(\(year) paste0("PM25_", sector)%>%
            map_dfr(~ Grid_rejoin_China%>%
                      filter(Year == year, Subarea%in%c("N","E"))%>%
                      dplyr::select(-GDPQ)%>%
                      left_join(Grid_sr_divide)%>%
                      na.omit()%>%
                      local_popw_summarise(group_main = "GDPQ",
                                           group_add = c("Pop_group","Year"),
                                           var = .x, prefix = FALSE)%>%
                      mutate(Sector = .x)))
GDP_source_summarise_sr_China
# colors <- colorspace::sequential_hcl(palette = "RdPu", n = 6)%>%rev()
colors <-  RColorBrewer::brewer.pal(6, "YlGn")[-6]
colors <- c("#feedde","#fed976","#feb24c","#fd8d3c","#fc4e2a")
speciesvar <- speciesvar_meic
year_seq <- c(2001,2013,2019)

for(s in c("residential","industry")){
  GDP_ineq_sr_China <- GDP_source_summarise_sr_China%>%
    filter(Pop_group == 2)%>%
    Popw_source_ineq(group_overall = "Year", source_col = "Sector")
  
  PM25_change_by_gdppc <- GDP_ineq_sr_China%>%
    filter(Sector == str_to_title(s), Year%in%year_seq)%>%
    arrange(GDPQ,Year)%>%
    group_by(GDPQ)%>%
    mutate(value = Popw - dplyr::lag(Popw))%>%
    ungroup()%>%
    filter(Year != year_seq[1])
  
  for(i in 1:(length(year_seq)-1)){
    PM25_change_by_gdppc <- PM25_change_by_gdppc%>%
      mutate(Year = ifelse(Year == year_seq[i+1],paste0(year_seq[i],"-",year_seq[i+1]),Year))
  }

  PM25_change_by_gdppc_plot <- box_series_plot(PM25_change_by_gdppc, x = "Year", y = "value", errorbar = FALSE,
                                               position = position_dodge(width = .9), fill = "GDPQ")+
    labs(x = NULL, y = expression(paste("Change in exposure (μg/m"^"3",")")))+
    # scale_x_discrete(labels = labels)+
    scale_fill_manual(name = NULL, values = colors)+
    scale_y_continuous(position = 'right')+
    theme(legend.key.width = unit(2, "cm"),
          legend.key.height = unit(2, "cm"),  ##legend size
          legend.position = "top",
          legend.text = element_text(colour = "black", size = rel(2.8),margin = margin(l = 10,r = 10)),
          # legend.title = element_text(colour = "black", size = rel(3.2)),
          axis.text = element_text(colour = "black",size = rel(2.8)),
          axis.title = element_text(colour = "black",size = rel(3.2)),
          panel.border = element_rect(colour = "black",linewidth = 3),
          plot.margin = margin(0,20,0,10))
  
  CairoPNG(paste0("result/Inequality/Emission_analysis/PM2.5_change_by_gdppc_",s,"_sector_CEC_area_time_series.png"),
           width = 2500, height = 1500, res = 120)
  print(PM25_change_by_gdppc_plot)
  dev.off()
}


##======Time series of disparity indexes in sectoral PM2.5======
## Inequality index for different sources derived from grids in specific regions
GDP_source_summarise_sr_China <- 2001:2019%>%
  map_dfr(\(year) paste0("PM25_", sector)%>%
            map_dfr(~ Grid_rejoin_China%>%
                      filter(Year == year, Subarea%in%c("N","E"))%>%
                      dplyr::select(-GDPQ)%>%
                      left_join(Grid_sr_divide)%>%
                      na.omit()%>%
                      local_popw_summarise(group_main = "GDPQ",
                                           group_add = c("Pop_group","Year"),
                                           var = .x, prefix = FALSE)%>%
                      mutate(Sector = .x)))
GDP_source_summarise_sr_China
s = "Residential"
for(s in c("Residential","Industry")){
  GDP_ineq_sr_China <- GDP_source_summarise_sr_China%>%
    filter(Pop_group == 2)%>%
    Popw_source_ineq(group_overall = "Year", source_col = "Sector")%>%
    filter(Sector == s)
  
  # GDP_res_ineq_sr_China%>%
  #   filter(Year == 2001)
  limits <- case_when(s=="Residential"~c(8,25),
                      s=="Industry"~c(10,35))
  
  GDP_ineq_sr_ts_China <- GDP_ineq_sr_China%>%
    Source_ineq_box_ts_vis(x = "Year", y = "Popw_ineq", fill = "Popw",
                           facet = "GDPQ", limits = limits, breaks = c(2005,2015),
                           fill_title = expression(paste("PM"[2.5]," exposure (μg/m"^"3",")")),
                           title = "Exposure Disparity (%)")+
    scale_x_continuous(limits = c(2000,2020), breaks = c(2005,2015))+
    theme(panel.border = element_rect(colour = "black",linewidth = 1.8),
          strip.text.x = element_text(
            size = 35, color = "black",margin = margin(t = 0.1, r = -.1, b = 0.3, l = -.1, "cm"), vjust = .5),
          legend.position = "right",
          legend.title = element_text(colour = "black", size = rel(3.2),hjust = .5, vjust = 1),
          axis.title = element_text(colour = "black",size = rel(3.2))
          )
  
  CairoPNG(paste0("result/Inequality/Emission_analysis/Exposure_inequality_",s,"_sector_CEC_area_time_series_China.png"),height = 1500,width = 5000, res = 200)
  print(GDP_ineq_sr_ts_China)
  dev.off()
}


