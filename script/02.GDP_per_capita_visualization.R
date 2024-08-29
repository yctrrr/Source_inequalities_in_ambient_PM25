setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
library(tidyverse)
library(data.table)
library(terra)
library(exactextractr)
library(Cairo)
library(ggplot2)
library(patchwork)

source('D:/shaoyanchuan/codebook/function/Sf_visualization.R')
source('D:/shaoyanchuan/codebook/function/Grid_visualization.R')
source('D:/shaoyanchuan/codebook/function/Health_calcu_core.R')
source('D:/shaoyanchuan/codebook/function/Time_series_analysis.R')
source("script/Inequality_cal_function.R")

load_shp()
dat_dir <- "D:/shaoyanchuan/data/"
sector <- c("agriculture","industry","power","residential","transportation")
metvar <- c("u10","v10","t2m","blh","tcc","sp","e","lsp")

Prov_info <- fread(paste0("data/input/Province_information.csv"))%>%
  dplyr::select(ProvID,engname,Subarea)
Prov_sf_China <- polygon_sf_China%>%
  dplyr::select(1)%>%
  setnames(c("ProvID","geometry"))%>%
  left_join(Prov_info)
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

Prov_sf_CEC <- Prov_sf_China%>%
  filter(Subarea %in% c("N","E"))
Prov_centroids_CEC <- Prov_centroids%>%
  filter(Subarea %in% c("N","E"))

County_sf_China <- readRDS(paste0(dat_dir,"Shp/county_China/county_shp.Rds"))%>%
  dplyr::select(-pop_2000, -pop_2010, -pop_2020)%>%
  rename(CountyID = PAC)
County_sf_rejoin_China <- fread("result/Inequality/Economic_imapct_county_dataset_China.csv", encoding = 'UTF-8')%>%
  left_join(County_sf_China%>%dplyr::select(CountyID),.)%>%na.omit()

Grid_rejoin_China <- fread(paste0("result/Inequality/Economic_imapct_grid_dataset_China.csv"))%>%
  left_join(Prov_info)
(Grid_rejoin_China$GDP_USD*1000)%>%summary

Grid_divide <- fread(paste0("result/Inequality/Grid_division_2001-2019.csv"))
Grid_sr_divide <- fread(paste0("result/Inequality/Grid_division_CEC_area_2001-2019.csv"))
Grid_SY_divide <- fread(paste0("result/Inequality/Grid_single_year_division_2001-2019.csv"))
Grid_sr_SY_divide <- fread(paste0("result/Inequality/Grid_single_year_division_CEC_area_2001-2019.csv"))

GDP_plot_set <- function(base_map){
  plot <- base_map+
    theme(panel.grid.major = element_blank(),
          plot.title = element_text(size = 50, hjust = .5),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.text = element_text(color = "black", size = 35),
          legend.title = element_text(color = "black", size = 40),
          legend.position = "right")
  return(plot)
}

Ineq_cal <- function(data, measurev, wt = "Pop"){
  data <- data%>%na.omit%>%as.data.frame()
  measurev%>%
    map_dfr(~ data.frame(Index = .x, GridN = nrow(data),
                         Gini = DescTools::Gini(data[,.x],data[,wt],unbiased = FALSE))
            )
}

GDPpc_facet_plot <- function(data, year_seq, area = "China", sf_data = NULL, breaks,colors,title){
  labels <- label_gen(breaks)
  
  data_mod <- data%>%
    mutate_at(c("value"), ~ cut(.,breaks,labels))
  
  if(area == "China"){
    data_mod_plot <- data_mod%>%
      grid_facet_plot(x = "X_Lon", y = "Y_Lat", value = "value",
                      facet = "Year", facet_text = 30)%>%
      grid_plot_set(area = "China", colors = colors,
                    title = NULL,
                    panel_setting = "facet")+
      scale_x_continuous(breaks = c(80,100,120))+
      coord_sf(xlim = c(-2825586,2406965), ylim = c(1800000,5940000))+
      labs(title = title)+
      theme(panel.grid.major = element_blank(),
            axis.text = element_text(color = "black", size = 25),
            plot.title = element_text(color = "black", size = 30),
            plot.margin = margin(t = 15, l = 15))
    
    nanhai_plot <- data_mod%>%
      split(f = ~ Year)%>%
      purrr::map(~ (.x%>%
                      grid_plot(x = "X_Lon", y = "Y_Lat", value = "value")%>%
                      grid_plot_set(area = "China", colors = colors,
                                    title = NULL,
                                    panel_setting = "facet")+
                      theme(panel.grid.major = element_blank())+
                      labs(title = NULL))%>%
                   nanhai_zoom_plot(facet_df = data.frame(Year = unique(.x$Year)),
                                    border.line.size = .6)
      )
    
    overall_plot <- (data_mod_plot+nanhai_plot)  
  }else{
    data_mod_plot <- data_mod%>%
      grid_facet_plot(x = "X_Lon", y = "Y_Lat", value = "value",
                      facet = "Year", facet_text = 30)%>%
      grid_plot_set(area = "non", colors = colors,
                    title = NULL,
                    panel_setting = "facet")+
      geom_sf(data = sf_data, fill = "white", alpha = 0.1, linewidth = 0.2)+
      labs(title = title)+
      theme(panel.grid.major = element_blank(),
            axis.text = element_text(color = "black", size = 25),
            plot.title = element_text(color = "black", size = 30),
            plot.margin = margin(t = 15, l = 15))
    overall_plot <- data_mod_plot
  }
  
  return(overall_plot)
}

##======Chinese population base map======
Grid_rejoin_China%>%
  dplyr::select(-GDPQ)%>%
  left_join(Grid_divide,.)%>%
  filter(Pop_group == 2)%>%
  group_by(GDPQ, Year)%>%
  summarise(Pop = sum(Pop))%>%
  ungroup()%>%
  filter(GDPQ == 5)
Pop_overall <- Grid_rejoin_China%>%
  group_by(X_Lon, Y_Lat)%>%
  summarise(Pop = mean(Pop))%>%
  ungroup()%>%
  mutate(Pop_group = ifelse(Pop <= 500,1,2))

Pop_CEC <- Pop_overall%>%
  dplyr::select(-Pop_group)%>%
  left_join(Grid_sr_divide)%>%
  na.omit()

breaks <- c(-Inf,5,10,20,40,60,80,100,+Inf)
labels <- label_gen(breaks)
colors <- rev(colorspace::sequential_hcl(length(labels), palette = "Terrain"))

Pop_overall_proj <- Pop_overall%>%
  filter(Pop_group == 2)%>%
  dplyr::select(X_Lon,Y_Lat,Pop)%>%
  dfproj(toCRS = CRS_China, method = "sum")%>%
  mutate(Pop = Pop/1e3, Pop = cut(Pop,breaks,labels))

Pop_overall_base_plot <- Pop_overall_proj%>%
  grid_plot(x = "x", y = "y", value = "Pop")%>%
  grid_plot_set(area = "no", value_type = "discrete", colors = colors,
                title = NULL,
                panel_setting = "single")+
  geom_sf(data = polygon_sf_China, fill = "white", alpha = 0.1, linewidth = 0.4)+
  geom_sf(data = Nline_sf_China, fill = NA, linewidth = 0.5)+
  xlim(-2879760, 2300893)+ylim(1800000,5940000)+
  theme(panel.grid.major = element_blank(),
        plot.title = element_text(size = 40, hjust = 28, vjust = 2.5),
        axis.text = element_text(color = "black", size = 35),
        legend.text = element_text(color = "black", size = 35),
        legend.title = element_blank(),
        legend.position = "left"
        )+
  labs(title = expression(paste("(a) Population (thousand pepole/ 100km"^"2",")")))
  # labs(title = expression(paste("Population distribution")))

Pop_overall_nanhai_plot <- (Pop_overall_base_plot+
  labs(title = NULL))%>%
  nanhai_zoom_plot(border.line.size = 1.2)

Pop_overall_plot <- (grid.draw(Pop_overall_base_plot + 
                                 geom_text(data = Prov_centroids, aes(x = X,y = Y,label = engname), size = 6)+
                                 Pop_overall_nanhai_plot))
ggplot_oupt(Pop_overall_plot, filename = paste0("result/GDP_pc/Population_distribution_average_China_2001-2019.png"),
            export_type = "Cairo", plot_type = "local", whr = c(2500,1600,160))


## population distribution for CEC
Prov_centroids_CEC
Pop_CEC_proj <- Pop_CEC%>%
  dplyr::select(X_Lon,Y_Lat,Pop)%>%
  dfproj(toCRS = CRS_China, method = "sum")%>%
  mutate(Pop = Pop/1e3, Pop = cut(Pop,breaks,labels))

Pop_CEC_plot <- Pop_CEC_proj%>%
  grid_plot(x = "x", y = "y", value = "Pop")%>%
  grid_plot_set(area = "no", value_type = "discrete", colors = colors,
                title = expression(paste("Population (thousand pepole/ 100km"^"2)")),
                panel_setting = "single")+
  geom_sf(data = Prov_sf_CEC, fill = "white", alpha = 0.1, linewidth = 0.4)+
  xlim(480000,1920000)+
  geom_text(data = Prov_centroids_CEC, aes(x = X,y = Y,label = engname), size = 10)+
  theme(panel.grid.major = element_blank(),
        plot.title = element_text(size = 42, hjust = .5, vjust = 1.2),
        axis.text = element_text(color = "black", size = 35),
        legend.position = "none")

ggplot_oupt(Pop_CEC_plot, filename = paste0("result/GDP_pc/Population_distribution_average_CEC_2001-2019.png"),
            export_type = "Cairo", plot_type = "local", whr = c(1500,1600,160))

##======GDP per capita ranking at grid level======
GDP_pc_complete <- Grid_rejoin_China%>%
  dplyr::select(X_Lon,Y_Lat,GridID,Year,GDP_USD,Pop)%>%
  left_join(Grid_divide,.)%>%
  dplyr::select(-GDPQ)
GDP_pc_complete
# GDP_pc_gini <- 2001:2019%>%
#   map_dfr(~ GDP_pc_complete%>%
#             filter(Year == .x)%>%
#             Ineq_cal("GDP_USD","Pop")%>%
#             mutate(Year = .x))
colors <- RColorBrewer::brewer.pal(6, "YlGn")[-6]
# colors <- c("#f2f0f7","#dadaeb","#fed976","#feb24c","#fd8d3c")
colors <- c("#feedde","#fed976","#feb24c","#fd8d3c","#fc4e2a")
for(year in c(20012019)){
  if(year == 20012019){
    GDP_pc_ranking <- GDP_pc_complete%>%
      filter(Pop_group == 2)%>%
      group_by(X_Lon,Y_Lat)%>%
      summarise(GDP_USD = mean(GDP_USD))%>%
      ungroup()
    GDP_pc_ranking_CEC <- GDP_pc_complete%>%
      filter(Subarea %in% c("N","E"))%>%
      group_by(X_Lon,Y_Lat)%>%
      summarise(GDP_USD = mean(GDP_USD))%>%
      ungroup()
  }else{
    GDP_pc_ranking <- GDP_pc_complete%>%
      filter(Year == year, Pop_group == 2)%>%
      dplyr::select(X_Lon,Y_Lat,GDP_USD)
    GDP_pc_ranking_CEC <- GDP_pc_complete%>%
      filter(Year == year, Subarea %in% c("N","E"))%>%
      dplyr::select(X_Lon,Y_Lat,GDP_USD)
  }

  GDP_pc_ranking_proj <- GDP_pc_ranking%>%
    dfproj(toCRS = CRS_China, method = "average")%>%
    GroupQ_mutate("GDP_USD")%>%
    rename(GDPQ = GDP_USD_group)%>%
    GDPQ_rename(var = "GDPQ")
  
  GDP_pc_ranking_CEC_proj <-  GDP_pc_ranking_CEC%>%
    dfproj(toCRS = CRS_China, method = "average")%>%
    GroupQ_mutate("GDP_USD")%>%
    rename(GDPQ = GDP_USD_group)%>%
    GDPQ_rename(var = "GDPQ")

  GDP_pc_ranking_plot <- GDP_pc_ranking_proj%>%
    grid_plot(x = "x", y = "y", value = "GDPQ")%>%
    grid_plot_set(area = "China", value_type = "discrete", colors = colors,
                  title = "GDP per capita ranking",
                  panel_setting = "single")
  
  GDP_pc_ranking_nanhai_plot <- (GDP_pc_ranking_plot+
                                   labs(title = NULL))%>%
    nanhai_zoom_plot(border.line.size = 1.2)
  
  GDP_pc_ranking_CEC_plot <- GDP_pc_ranking_CEC_proj%>%
    grid_plot(x = "x", y = "y", value = "GDPQ")%>%
    grid_plot_set(area = "no", value_type = "discrete", colors = colors,
                  title = "GDP per capita ranking",
                  panel_setting = "single")+
    geom_sf(data = Prov_sf_CEC, fill = "white", alpha = 0.1, linewidth = 0.4)

  if(year == 20012019){
    GDP_pc_ranking_plot <- GDP_pc_ranking_plot%>%
      grid_plot_set(area = "China", value_type = "discrete", colors = colors,
                    title = NULL,
                    panel_setting = "single")+
      geom_text(data = Prov_centroids, aes(x = X,y = Y,label = engname), size = 6)+
      theme(axis.text = element_text(color = "black", size = 35),
            legend.text = element_text(color = "black", size = 35),
            legend.title = element_blank(),
            legend.position = "top",
            panel.grid.major = element_blank())
    
    GDP_pc_ranking_CEC_plot <- GDP_pc_ranking_CEC_plot+
      geom_text(data = Prov_centroids_CEC, aes(x = X,y = Y,label = engname), size = 10)+
      xlim(480000,1920000)+
      theme(axis.text = element_text(color = "black", size = 35),
            legend.text = element_text(color = "black", size = 35),
            legend.title = element_blank(),
            legend.position = "none",
            panel.grid.major = element_blank())
  }else{
    GDP_pc_ranking_plot <- GDP_pc_ranking_plot%>%
      GDP_plot_set()
    GDP_pc_ranking_CEC_plot <- GDP_pc_ranking_CEC_plot%>%
      GDP_plot_set()
  }
  
  GDP_pc_ranking_plot <- (grid.draw(GDP_pc_ranking_plot + GDP_pc_ranking_nanhai_plot))
  ggplot_oupt(GDP_pc_ranking_plot,
              filename = paste0("result/GDP_pc/GDP_per_capita_spatial_ranking_China_",year,".png"),
              export_type = "Cairo", plot_type = "local", whr = c(2500,1600,160))
  
  ggplot_oupt(GDP_pc_ranking_CEC_plot,
              filename = paste0("result/GDP_pc/GDP_per_capita_spatial_ranking_CEC_",year,".png"),
              export_type = "Cairo", plot_type = "local", whr = c(1500,1600,160))
}


##======Log scale for GDP per capita at the grid level======
region = "Overall"
year = 2001
year_seq <- c(2001,2013,2019)
Grid_rejoin_China%>%names
for(region in c("Overall","CEC")){
  if(region == "Overall"){
    data <- Grid_rejoin_China%>%
      dplyr::select(-GDPQ)%>%
      left_join(Grid_divide)%>%
      filter(Pop_group == 2, Year%in%year_seq)%>%
      dplyr::select(-GDPQ,-ProvID,-engname,-Pop_group)%>%
      left_join(Grid_SY_divide)%>%
      mutate_at("GDP_USD",~ ifelse(.x <=0, 0.01, .x))%>%
      na.omit()%>%
      as.data.frame()
    
    area = "China"
    poly_sf_sel <- Prov_sf_China
    whr <- c(3000,1100,155)
    legend_title <- "(b) Log scale for GDP per capita (thousand USD/person)"

  }else if(region == "CEC"){
    data <- Grid_rejoin_China%>%
      dplyr::select(-GDPQ,-ProvID,-engname)%>%
      left_join(Grid_sr_SY_divide)%>%
      filter(Subarea%in%c("N","E"), Year%in%year_seq)%>%
      na.omit()%>%
      as.data.frame()
    area = "none"
    poly_sf_sel <- Prov_sf_CEC
    whr <- c(2000,1100,155)
    legend_title <- "(c) Log scale for GDP per capita (thousand USD/person)"

  }
  
  breaks <- c(-Inf, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, +Inf)
  label <- label_gen(breaks)
  # legend_colors <- RColorBrewer::brewer.pal(length(label),"Set1")
  legend_colors <- colorspace::sequential_hcl(length(label),"Heat")%>%rev()
  # colorspace::hcl_palettes()%>%plot
  lab_title <- "GDP per capita"
  
  gdppc_df <- year_seq%>%
    map_dfr(~ data%>%
              mutate(GDP_USD = log(GDP_USD))%>%
              filter(Year == .x)%>%
              dplyr::select(X_Lon,Y_Lat,GDPQ, GDP_USD)%>%
              dfproj(toCRS = CRS_China, method = "average")%>%
              rename(X_Lon = x, Y_Lat = y, value = GDP_USD)%>%
              # filter(GDPQ > 4)%>%
              mutate(Year = .x))
  gdppc_df%>%
    filter(Year == 2013)%>%rast()%>%plot
  # top_gdppc_df%>%filter(Year == 2001)%>%rast%>%plot
  
  gdppc_plot <- gdppc_df%>%
    GDPpc_facet_plot(year_seq = year_seq, area = area, sf_data = poly_sf_sel,
                     breaks = breaks, colors = legend_colors, title = legend_title)
  
  ggplot_oupt(gdppc_plot,
              filename = paste0("result/GDP_pc/Spatial_distribution_of_top_GDP_per_capita_",region,".png"),
              export_type = "Cairo", plot_type = "local",
              whr = whr)
  
}


##======GDP per capita change at county level======
# Deflator <- fread(paste0(dat_dir,"Economic/DOSE-global_economic_output/Deflator.csv"))%>%
#   dplyr::select(-cpi_2015,-PPP,-cpi_usd)%>%
#   filter(country == "China")%>%
#   dplyr::select(-country,-deflator_2015)%>%
#   rename(Year = year)

## We deprecated mappings with raw county data cause only few counties have data both at start and end years
for(period_index in list(c(2001,2019),c(2013,2019))){
  # period_index <- c(2001,2019)
  period_len <- period_index[2] - period_index[1]

  ##  completed county GDP per cpaita with USD dollar
  GDP_pc_county_complete <- County_sf_rejoin_China%>%
    st_drop_geometry()%>%
    dplyr::select(CountyID,Year,GDP_USD)
  count(GDP_pc_county_complete, CountyID)%>%summary

  ##  The growth rate of GDP per capita (compared with the last year)
  ##  The growth rate will be calculated as multiyear averages
  GDP_pc_rate_county <- GDP_pc_county_complete%>%
    arrange(CountyID,Year)%>%
    mutate(GDP_rate = 100*(GDP_USD - dplyr::lag(GDP_USD, 1))/dplyr::lag(GDP_USD, 1))%>%
    filter(Year%in%c(period_index[1]:period_index[2]), Year != period_index[1])%>%
    group_by(CountyID)%>%
    summarise(GDP_rate_mean = mean(GDP_rate),
              GDP_rate_sd = sd(GDP_rate))%>%
    ungroup()

  ##  Period change
  GDP_pc_change_county <- GDP_pc_county_complete%>%
    filter(Year%in%period_index)%>%
    dplyr::select(CountyID,Year,GDP_USD)%>%
    pivot_wider(names_from = Year, values_from = GDP_USD)%>%
    mutate(GDP_pc_change =  .[[3]]- .[[2]])

  ##  GDP per capita change plot
  if(period_len > 10){
    breaks <- c(-Inf,2,3,4,5,6,7,8,+Inf)
  }else{
    breaks <- c(-Inf,0.5,1,1.5,2,2.5,3,3.5,+Inf)
  }
  labels <- label_gen(breaks)
  colors <- rev(RColorBrewer::brewer.pal(length(labels),"RdYlBu"))
  
  GDP_pc_change_county_base_plot <- GDP_pc_change_county%>%
    mutate_at("GDP_pc_change", ~ cut(.x,breaks,labels))%>%
    left_join(County_sf_China,.)%>%
    st_transform(CRS_China)%>%
    na.omit()%>%
    sf_polygon_plot(value = "GDP_pc_change")%>%
    grid_plot_set(area = "China", value_type = "discrete",
                  colors = colors, title = NULL,
                  panel_setting = "single")%>%
    GDP_plot_set()+
    labs(title = "GDP per capita change")
  
  GDP_pc_change_county_nanhai_plot <- (GDP_pc_change_county_base_plot+
    labs(title = NULL))%>%
    nanhai_zoom_plot(border.line.size = 1.2)
  
  GDP_pc_change_county_plot <- (grid.draw(GDP_pc_change_county_base_plot + GDP_pc_change_county_nanhai_plot))
  
  for(v in c("GDP_rate_mean","GDP_rate_sd")){
    breaks <- case_when(v == "GDP_rate_mean"~c(-Inf,10,12,14,16,18,20,+Inf),
                        v == "GDP_rate_sd"~c(-Inf,10,15,20,25,30,35,+Inf))
    labels <- label_gen(breaks)
    
    if(v == "GDP_rate_mean"){
      colors <- rev(colorspace::sequential_hcl(length(labels), palette = "BurgYl"))
      title <- "GDP per capita growth percent"
    }else{
      colors <- rev(colorspace::sequential_hcl(length(labels), palette = "Purp"))
      title <- "Standard deviation"
    }

    GDP_pc_rate_county_base_plot <- GDP_pc_rate_county%>%
      mutate_at(all_of(v),~cut(.x,breaks,labels))%>%
      left_join(County_sf_China,.)%>%
      st_transform(CRS_China)%>%
      na.omit()%>%
      sf_polygon_plot(value = v)%>%
      grid_plot_set(area = "China",value_type = "discrete", colors = colors,
                    title = NULL,
                    panel_setting = "single")%>%
      GDP_plot_set()+
      labs(title = title)
    
    GDP_pc_rate_county_nanhai_plot <- (GDP_pc_rate_county_base_plot+
                                         labs(title = NULL))%>%
      nanhai_zoom_plot()
    
    GDP_pc_rate_county_plot <- (grid.draw(GDP_pc_rate_county_base_plot + GDP_pc_rate_county_nanhai_plot))
    
    assign(paste0(v,"_county_plot"),GDP_pc_rate_county_plot)
  }
  GDP_pc_county_plot <- GDP_pc_change_county_plot+
    plot_spacer()+
    GDP_rate_mean_county_plot+
    plot_spacer()+
    GDP_rate_sd_county_plot+
    plot_layout(widths = c(8,.5,8,.5,8), guides = "keep")

  CairoPNG(paste0("result/GDP_pc/GDP_per_capita_county_overall_plot_",
                  period_index[1],"-",period_index[2],".png"),
           width = 4800, height = 1200, res = 140)
  print(GDP_pc_county_plot)
  dev.off()
}

##======GDP per capita change at grid level======
GDP_pc_complete <- Grid_rejoin_China%>%
  dplyr::select(X_Lon,Y_Lat,GridID,Year,GDP_USD)
GDP_pc_complete

for (period_index in list(c(2001,2019),c(2013,2019))){
  ##  Period change
  # period_index <- c(2013,2019)
  period_len <- period_index[2] - period_index[1]

  GDP_pc_rate <- GDP_pc_complete%>%
    arrange(GridID,Year)%>%
    mutate(GDP_rate = 100*(GDP_USD - dplyr::lag(GDP_USD, 1))/dplyr::lag(GDP_USD, 1))%>%
    filter(Year%in%c(period_index[1]:period_index[2]), Year != period_index[1])%>%
    group_by(X_Lon,Y_Lat,GridID)%>%
    summarise(GDP_rate_mean = mean(GDP_rate),
              GDP_rate_sd = sd(GDP_rate))%>%
    ungroup()

  GDP_pc_change_proj <- GDP_pc_complete%>%
    left_join(Pop_overall)%>%
    filter(Year%in%period_index, Pop_group == 2)%>%
    dplyr::select(X_Lon,Y_Lat,Year,GDP_USD)%>%
    pivot_wider(names_from = Year, values_from = GDP_USD)%>%
    mutate(GDP_pc_change = .[[4]]- .[[3]])%>%
    dfproj(toCRS = CRS_China, method = "average")

  if(period_len > 10){
    breaks <- c(-Inf,2,3,4,5,6,7,8,+Inf)
  }else{
    breaks <- c(-Inf,0.5,1,1.5,2,2.5,3,3.5,+Inf)
  }
  labels <- label_gen(breaks)
  colors <- rev(RColorBrewer::brewer.pal(length(labels),"RdYlBu"))

  GDP_pc_change_base_plot <- GDP_pc_change_proj%>%
    mutate(GDP_pc_change = cut(GDP_pc_change,breaks,labels))%>%
    grid_plot(x = "x", y = "y", value = "GDP_pc_change")%>%
    grid_plot_set(area = "China",value_type = "discrete", colors = colors,
                  title = NULL,
                  panel_setting = "single")%>%
    GDP_plot_set()+
    labs(title = "GDP per capita change")
  
  GDP_pc_change_nanhai_plot <- (GDP_pc_change_base_plot+
    labs(title = NULL))%>%
    nanhai_zoom_plot()

  GDP_pc_change_plot <- (grid.draw(GDP_pc_change_base_plot + GDP_pc_change_nanhai_plot))

  ##  Economic growth rate
  for(v in c("GDP_rate_mean","GDP_rate_sd")){
    GDP_pc_rate_proj <- GDP_pc_rate%>%
      left_join(Pop_overall)%>%
      filter(Pop_group == 2)%>%
      dplyr::select(X_Lon,Y_Lat,GDP_rate_mean,GDP_rate_sd)%>%
      dfproj(toCRS = CRS_China, method = "average")

    breaks <- case_when(v == "GDP_rate_mean"~c(-Inf,10,12,14,16,18,20,+Inf),
                        v == "GDP_rate_sd"~c(-Inf,10,15,20,25,30,35,+Inf))
    labels <- label_gen(breaks)
    if(v == "GDP_rate_mean"){
      colors <- rev(colorspace::sequential_hcl(length(labels), palette = "BurgYl"))
      title <- "GDP per capita growth percent"
    }else{
      colors <- rev(colorspace::sequential_hcl(length(labels), palette = "Purp"))
      title <- "Standard deviation"
    }

    GDP_pc_rate_base_plot <- GDP_pc_rate_proj%>%
      mutate_at(all_of(v), ~cut(.x, breaks, labels))%>%
      grid_plot(x = "x", y = "y", value = v)%>%
      grid_plot_set(area = "China",value_type = "discrete", colors = colors,
                    title = NULL,
                    panel_setting = "single")%>%
      GDP_plot_set()+
      labs(title = title)

    GDP_pc_rate_nanhai_plot <- (GDP_pc_rate_base_plot+
      labs(title = NULL))%>%
      nanhai_zoom_plot()

    GDP_pc_rate_plot <- (grid.draw(GDP_pc_rate_base_plot + GDP_pc_rate_nanhai_plot))

    assign(paste0(v,"_plot"), GDP_pc_rate_plot)
  }
  
  GDP_pc_grid_plot <- GDP_pc_change_plot+
    plot_spacer()+
    GDP_rate_mean_plot+
    plot_spacer()+
    GDP_rate_sd_plot+
    plot_layout(widths = c(15,.5,15,.5,15), guides = "keep")

  CairoPNG(paste0("result/GDP_pc/GDP_per_capita_grid_overall_plot_",
                  period_index[1],"-",period_index[2],".png"),
           width = 4800, height = 1200, res = 140)
  print(GDP_pc_grid_plot)
  dev.off()
}

##======Other summary======
Grid_rejoin_China$Year%>%unique
Grid_ave_China <- Grid_rejoin_China%>%
  group_by(X_Lon, Y_Lat)%>%
  summarise(Pop = mean(Pop), GDP_USD = mean(GDP_USD))%>%
  ungroup()%>%
  mutate(Pop_group = ifelse(Pop <= 500,1,2))

Grid_ave_China%>%
  filter(Pop_group == 2)%>%
  cor()

##  source simulations from extended SHAP values
##  remove the simulations in Taiwan
sector <- c("agriculture","industry","power","residential","transportation")
PM_source_extend_China <- fread("result/Interpretation/SHAP_China_TAP/PM_Source_2001-2019_TAP_China.csv")%>%
  left_join(Grid_info_China,.)%>%
  filter(!ProvID%in%710000)

## Change in source contributions during 2001-2019
PM_source_change_China <- PM_source_extend_China%>%
  filter(Year%in%c(2001,2019))%>%
  mutate_at(paste0("PM25_",sector), ~ifelse(.x < 0,0,.x))%>%
  group_by(X_Lon, Y_Lat, GridID)%>%
  summarise_at(paste0("PM25_",sector), ~.x[Year==2019]-.x[Year==2001])%>%
  ungroup()

PM_source_change_China

Grid_ave_China%>%
  left_join(PM_source_change_China)%>%
  filter(Pop_group == 2)%>%cor
