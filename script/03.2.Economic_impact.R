setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
library(plm)
library(fixest)
library(lme4)
library(lmerTest)
library(optimx)
library(tidyverse)
library(data.table)
library(terra)
library(nlme)
library(gamm4)
library(furrr)
library(Epi)
source('D:/shaoyanchuan/codebook/function/Grid_visualization.R')
source('D:/shaoyanchuan/codebook/function/Health_calcu_core.R')
source('D:/shaoyanchuan/codebook/function/Time_series_analysis.R')
source('D:/shaoyanchuan/codebook/function/Sf_visualization.R')

dat_dir <- "D:/shaoyanchuan/data/"
sector <- c("agriculture","industry","power","residential","transportation")
metvar <- c("u10","v10","t2m","blh","tcc","sp","e","lsp")

County_sf_China <- readRDS(paste0(dat_dir,"Shp\\county_China\\county_shp.Rds"))%>%
  dplyr::select(-pop_2000, -pop_2010, -pop_2020)%>%
  rename(CountyID = PAC)

County_sf_rejoin_China <- fread("result/Inequality/Economic_imapct_county_dataset_China.csv", encoding = 'UTF-8')%>%
  left_join(County_sf_China%>%dplyr::select(CountyID),.)%>%
  na.omit()

Prov_info <- fread(paste0("data/input/Province_information.csv"))%>%
  dplyr::select(ProvID,engname,Subarea)

Prov_sf_China <- polygon_sf_China%>%
  dplyr::select(1)%>%
  setnames(c("ProvID","geometry"))%>%
  left_join(Prov_info)

## load dataset
load_shp()

##======Preparing data======
County_rejoin_panel_China <- County_sf_rejoin_China%>%
  st_drop_geometry()%>%
  mutate(across(c(paste0("PM25_",sector),"GDP_USD","Pop","PM25"),~ ifelse(.x <= 0, 0.01, .x)))

Grid_rejoin_China <- fread(paste0("result/Inequality/Economic_imapct_grid_dataset_China.csv"))%>%
  left_join(Prov_info)

Grid_divide <- fread(paste0("result/Inequality/Grid_division_2001-2019.csv"))
Grid_sr_divide <- fread(paste0("result/Inequality/Grid_division_CEC_area_2001-2019.csv"))
Grid_SY_divide <- fread(paste0("result/Inequality/Grid_single_year_division_2001-2019.csv"))
Grid_sr_SY_divide <- fread(paste0("result/Inequality/Grid_single_year_division_CEC_area_2001-2019.csv"))

Grid_rejoin_panel_China <- Grid_rejoin_China%>%
  dplyr::select(-GDPQ)%>%
  left_join(Grid_divide)%>%
  na.omit()%>%
  mutate(across(c(paste0("PM25_",sector),"GDP_USD","Pop","PM25"),~ ifelse(.x <=0, 0.01, .x)),
         Interv = ifelse(Year < 2013, 0, 1))%>%
  filter(Year >= 2001, Year <= 2019)

summary(Grid_rejoin_panel_China)
names(Grid_rejoin_panel_China)
# DescTools::Gini(Grid_join_China_regroup$PM25, 
#                 weights = Grid_join_China_regroup$GDP_USD*Grid_join_China_regroup$Pop)

Grid_rejoin_panel_China$ProvID%>%unique
# colnames(Grid_join_China_regroup)
Grid_rejoin_panel_China%>%is.pbalanced(index = c("GridID","Year"))


##======Random slope for year-specific regression======
tv = as.character(look_up[n,2])
iv = "GDP_USD"
logtran = TRUE
sm_reg_by_year <- function(data, tv = "PM25", iv = "GDP_USD", logtran = TRUE){
  if(logtran){
    ivtran <- paste0("log(",iv,")")
    iv_s <- paste0("Ns(",ivtran,", df = 4, ref = ",mean(log(data[[iv]])),")")
    iv_s <- paste0("s(",ivtran,")")
  }else{
    ivtran <- iv
    iv_s <- paste0("Ns(",ivtran,", df = 4, ref = ",mean(data[[iv]]),")")
    iv_s <- paste0("s(",ivtran,")")
  }

  # fm <- as.formula(paste(tv, " ~ ",iv_s," + s(Pop)"))
  fm <- as.formula(paste(tv, " ~ ",iv_s))
  weight <- data[["Pop"]]
  GAM_model <- gamm4(fm, data = data, weight = weight)
  
  pred_df_extract(GAM_model[["gam"]], data, voi_log = TRUE, iv, iv_s)
}
iv
data%>%
  dplyr::select(all_of(c(tv, iv, "Year","X_Lon","Y_Lat")))%>%
  filter(PM25 > 37.3, PM25 < 37.5)
  summary

look_up <- expand_grid(Year = c(2001,2013,2019),
                       Sector = c("PM25", paste0("PM25_",sector)),
                       Region = c("Overall","CEC"))

## nonlinear effects for central eastern China
Grid_rejoin_panel_China
options(future.globals.maxSize = 10000 * 1024^3)
plan(multisession, workers = 6)
nonlin_extract <- 1:nrow(look_up)%>%
  furrr::future_map_dfr(function(n){
    library(Epi)
    region <- as.character(look_up[n,3])
    
    if(region == "Overall"){
      data <- Grid_rejoin_panel_China%>%
        filter(Pop_group == 2, Year == as.numeric(look_up[n,1]))%>%
        as.data.frame()
    }else if(region == "CEC"){
      data <- Grid_rejoin_panel_China%>%
        filter(Pop_group == 2, Subarea%in%c("N","E"), Year == as.numeric(look_up[n,1]))%>%
        as.data.frame()
    }
    
    sm_reg_by_year(data, as.character(look_up[n,2]), iv = "GDP_USD")%>%
      mutate(Year = as.numeric(look_up[n,1]),
             Sector = gsub("PM25_","",look_up[n,2]),
             Subarea = region)
  })
t = nonlin_extract%>%filter(Year == 2001, Sector == "industry")
summary(t$beta)
fwrite(nonlin_extract, file = paste0("result/Inequality/Economic_impacts/GDP_PM25_nonlinear_regression.csv"))

## Visualization
nonlin_extract <- fread(paste0("result/Inequality/Economic_impacts/GDP_PM25_nonlinear_regression.csv"))
colors <- RColorBrewer::brewer.pal(5, "RdYlGn")
colors[3] <- "#3182bd"
nonlin_extract%>%
  filter(beta == 0)
##  Relationship plot
for(region in c("Overall","CEC")){
  nonlin_plot <- nonlin_extract%>%
    filter(Year%in%c(2001,2013,2019), Subarea==region, Sector != "PM25")%>%
    mutate(Sector = factor(Sector, levels = sector),
           Sector = str_to_title(Sector))%>%
    plot_nonlin_gam(pred_df = ., custom = TRUE)+
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Sector), alpha = 0.3,
                show.legend = FALSE)+
    geom_path(aes(y = beta, color = Sector), size = 1.5)+
    facet_wrap(~ Year, nrow = 1, scales = "free_x")+
    scale_fill_manual(name = NULL, values = colors)+
    scale_color_manual(name = NULL, values = colors,
                       guide = guide_legend(frame.colour = "white", frame.linewidth = .2,
                                            title.position = "top", title.hjust = 0.5,
                                            byrow = TRUE, ticks = F,
                                            override.aes = list(linewidth = 2)))+
    labs(x = "Log scale for GDP per capita (thousand US dollars)",
         y = expression(paste("Change in PM"[2.5]," (μg/m"^"3",")")))+
    theme(strip.text.x = element_blank(),
          panel.border = element_rect(size = 2,colour = "black", fill = NA),
          legend.key.width = )

  # ggplot_oupt(nonlin_plot,
  #             filename = paste0("result/Inequality/Economic_impacts/By_year/GDP_source_nonlinear_plot.png"),
  #             export_type = "Cairo", plot_type = "local", whr = c(4000,1000,150))
  ggplot_oupt(nonlin_plot,
              filename = paste0("result/Inequality/Economic_impacts/GDP_source_nonlinear_plot_",region,".png"),
              export_type = "ggsave", plot_type = "local", whr = c(25,7.5,200))
}

##======Density plot======
colors <- RColorBrewer::brewer.pal(5, "OrRd")
colors <- c("#feedde","#fed976","#feb24c","#fd8d3c","#fc4e2a")
for(region in c("Overall","CEC")){
  dens_df <- c(2001,2013,2019)%>%
    map_dfr(function(year){
      if(region == "Overall"){
        data_sel <- Grid_rejoin_panel_China%>%
          filter(Pop_group == 2, Year == year)%>%
          mutate(GDP_USD = log(GDP_USD))
      }else if(region == "CEC"){
        data_sel <- Grid_rejoin_panel_China%>%
          filter(Pop_group == 2, Subarea%in%c("N","E"), Year == year)%>%
          mutate(GDP_USD = log(GDP_USD))
      }

      quantiles <- quantile(data_sel$GDP_USD, probs = seq(0,1,.2))
      quantiles[1] <- -Inf
      quantiles[length(quantiles)] <- +Inf
      
      GDP_USD_dens <- density(data_sel$GDP_USD, n = 4096)
      
      ## to align with other plots
      xmin <- min(data_sel$GDP_USD)
      xmax <- max(data_sel$GDP_USD)
      ymax <- max(GDP_USD_dens$y)
      
      dens_limit <- data.frame(x = c(xmin,xmax), y = c(0,0))
      
      data.frame(x = GDP_USD_dens$x, y = GDP_USD_dens$y)%>%
        full_join(dens_limit)%>%
        filter(x >= xmin, x <= xmax)%>%
        mutate(Year = year, 
               y = y*1/ymax,
               GDPQ = findInterval(.$x, quantiles))%>%
        GDPQ_rename("GDPQ")
    })

  # plot(dens_df$x,dens_df$y)
  # data_sel%>%
  #   ggplot(aes_string(x = "GDP_USD"))+
  #   geom_density(color = "black", alpha=0.8)

  density_plot <- dens_df%>%
    ggplot(aes_string(x = "x", y = "y"))+
    # geom_density(aes(y =..scale..),color = "black", alpha=0.8)+
    geom_area(aes(group = GDPQ, fill = GDPQ), color = "black")+
    facet_wrap(~ Year, nrow = 1, scales = "free_x")+
    theme_bw(base_size = 20)+
    labs(x = NULL ,y = "Density")+
    scale_fill_manual(name = NULL, values = colors,
                      guide = guide_legend(frame.colour = "white", frame.linewidth = .1,
                                           title.position = "top", title.hjust = 0.5,
                                           byrow = TRUE, ticks = F))+
    theme(strip.text.x = element_text(size = rel(2.5), color = "black",margin = margin(.1, -.1, .1, -.1, "cm")),
          strip.background = element_blank())+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(size = 1, fill = 'white', colour = 'black'),
          panel.border = element_rect(size = 2, colour = "black", fill = NA),
          panel.spacing = unit(2, "cm"))+
    theme(
          axis.text.y = element_text(colour = "black", size = rel(2)),
          axis.text.x = element_blank(),
          # axis.text = element_text(colour = "black",size = rel(1.5)),
          axis.ticks.x = element_blank(),
          axis.title = element_text(colour = "black",size = rel(2)),
          plot.margin = margin(0, 5, 0, 15))+
    theme(legend.key.width=unit(1, "cm"),
          legend.key.height=unit(1, "cm"),  ##legend size
          legend.spacing.y = unit(.8, 'cm'),
          legend.key.spacing.y = unit(.8, 'cm'),
          legend.key = element_blank(),
          legend.position = "right",
          legend.text = element_text(colour = "black",size = rel(1.5)),
          legend.title = element_blank(),
          legend.background = element_blank())
  
  ggplot_oupt(density_plot,
              filename = paste0("result/Inequality/Economic_impacts/GDP_by_year_density_plot_",region,".png"),
              export_type = "ggsave", plot_type = "local", whr = c(25,7.5,200))
}

##======Spatial distribution of the top rich groups======
# data = top_gdppc_df
# year_seq = c(2001,2013,2019)
# area = "none"
# sf_data = poly_sf_sel
# title = legend_title
# colors = legend_colors
GDPpc_facet_plot <- function(data,year_seq,area = "China",sf_data = NULL,breaks,colors,title){
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
            axis.text = element_text(color = "black", size = 20),
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
            axis.text = element_text(color = "black", size = 20),
            plot.title = element_text(color = "black", size = 30),
            plot.margin = margin(t = 15, l = 15))
    overall_plot <- data_mod_plot
  }

  return(overall_plot)
}

Abs_source_facet_grid_plot2 <- function(data, year_seq, breaks = c(-Inf,1,3,5,10,15,20,+Inf),
                                       colors, area = "China", sf_data = NULL, text_table = NULL,
                                       title = expression(paste("PM"[2.5]," concentration (μg/m"^"3",")"))){
  labels <- label_gen(breaks)

  
  samod <- data%>%
    mutate_at(c("value"), ~ cut(.,breaks,labels))
  if(area == "China"){
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
  }else{
    samod_plot <- samod%>%
      grid_facet_grid_plot(x = "X_Lon", y = "Y_Lat", value = "value",
                           row = "Year", col = "Sector", facet_text = 30)%>%
      grid_plot_set(area = area, colors = colors,
                    title = title,
                    panel_setting = "facet")+
      geom_sf(data = sf_data, fill = "white", alpha = 0.1, linewidth = 0.2)+
      scale_x_continuous(breaks = c(110,120))+
      theme(panel.grid.major = element_blank(),
            axis.text = element_text(color = "black", size = 20))
    samod_overall_plot <- samod_plot
  }

  return(samod_overall_plot)
}

region = "Overall"
year = 2001
year_seq <- c(2001,2013,2019)

for(region in c("Overall","CEC")){
  if(region == "Overall"){
    data <- Grid_rejoin_panel_China%>%
      filter(Pop_group == 2, Year%in%year_seq)%>%
      dplyr::select(-GDPQ,-ProvID,-engname,-Subarea,-Pop_group)%>%
      left_join(Grid_SY_divide)%>%
      na.omit()%>%
      as.data.frame()
    area = "China"
    poly_sf_sel <- Prov_sf_China
  }else if(region == "CEC"){
    data <- Grid_rejoin_panel_China%>%
      filter(Pop_group == 2, Subarea%in%c("N","E"), Year%in%year_seq)%>%
      dplyr::select(-GDPQ,-ProvID,-engname,-Subarea,-Pop_group)%>%
      left_join(Grid_sr_SY_divide)%>%
      na.omit()%>%
      as.data.frame()
    area = "none"
    poly_sf_sel <- Prov_sf_China%>%
      filter(Subarea %in% c("N","E"))
  }
  
  data%>%
    filter(Year == 2001, log(GDP_USD) > 1)%>%
    dplyr::select(X_Lon,Y_Lat,GDP_USD,PM25_residential)%>%
    rast%>%plot
  
  breaks <- c(-Inf, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, +Inf)
  label <- label_gen(breaks)
  # legend_colors <- RColorBrewer::brewer.pal(length(label),"Set1")
  legend_colors <- colorspace::sequential_hcl(length(label),"Heat")%>%rev()
  # colorspace::hcl_palettes()%>%plot
  lab_title <- "GDP per capita"
  legend_title <- "Log scale for GDP per capita (thousand US dollars/person)"
  
  top_gdppc_df <- year_seq%>%
    map_dfr(~ data%>%
              mutate(GDP_USD = log(GDP_USD))%>%
              filter(Year == .x)%>%
              dplyr::select(X_Lon,Y_Lat,GDPQ, GDP_USD)%>%
              dfproj(toCRS = CRS_China, method = "average")%>%
              rename(X_Lon = x, Y_Lat = y, value = GDP_USD)%>%
              # filter(GDPQ > 4)%>%
              mutate(Year = .x))

  # top_gdppc_df%>%filter(Year == 2001)%>%rast%>%plot
  
  top_gdppc_plot <- top_gdppc_df%>%
    GDPpc_facet_plot(year_seq = year_seq, area = area, sf_data = poly_sf_sel,
                     breaks = breaks, colors = legend_colors, title = legend_title)
  
  ggplot_oupt(top_gdppc_plot,
              filename = paste0("result/Inequality/Economic_impacts/Spatial_distribution_of_top_GDP_per_capita_",region,".png"),
              export_type = "Cairo", plot_type = "local",
              whr = c(3000,1100,155))
  
  
  breaks <- c(-Inf,1,3,5,10,15,25,35,+Inf)
  labels <- label_gen(breaks)
  colors <- rev(RColorBrewer::brewer.pal(length(labels), "Spectral"))
  
  PM_samod_multiyear_df <- year_seq%>%
    map_dfr(~ data%>%
              mutate_at(c(paste0("PM25_",sector)),~ifelse(.x < 0,0,.x))%>%
              filter(Year == .x)%>%
              dplyr::select(X_Lon,Y_Lat,GDPQ,paste0("PM25_",sector))%>%
              dfproj(toCRS = CRS_China, method = "average")%>%
              filter(GDPQ > 4)%>%
              dplyr::select(-GDPQ)%>%
              rename(X_Lon = x, Y_Lat = y)%>%
              pivot_longer(paste0("PM25_",sector), names_to = "Sector")%>%
              mutate(Sector = gsub("*PM25_","",Sector)%>%str_to_title(),
                     Year = .x))
  # PM_samod_multiyear_df%>%
  #   filter(Year == 2001, Sector == "Industry")%>%dplyr::select(1:2,4)%>%rast%>%plot
  PM_samod_multiyear_plot <- PM_samod_multiyear_df%>%
    Abs_source_facet_grid_plot2(year_seq, breaks, colors, area = area, sf_data = poly_sf_sel)
  
  ggplot_oupt(PM_samod_multiyear_plot, filename = paste0("result/Inequality/Economic_impacts/Spatial_distribution_of_top_GDP_per_capita_source_contribution_",region,".png"),
              export_type = "Cairo", plot_type = "local", whr = c(2600,2800,155))
  
}
