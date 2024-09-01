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

dat_dir <- "D:/shaoyanchuan/data/"
source('D:/shaoyanchuan/codebook/function/Grid_visualization.R')
source('D:/shaoyanchuan/codebook/function/Health_calcu_core.R')
source('D:/shaoyanchuan/codebook/function/Time_series_analysis.R')
source('D:/shaoyanchuan/codebook/function/Sf_visualization.R')

## load dataset
load_shp()
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

##======Random slope for year-specific regression======
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

look_up <- expand_grid(Year = c(2001,2013,2019),
                       Sector = c("PM25", paste0("PM25_",sector)),
                       Region = c("Overall","CEC"))

## nonlinear effects for central eastern China
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
fwrite(nonlin_extract, file = paste0("result/Inequality/Economic_impacts/GDP_PM25_nonlinear_regression.csv"))

## Visualization
nonlin_extract <- fread(paste0("result/Inequality/Economic_impacts/GDP_PM25_nonlinear_regression.csv"))
colors <- RColorBrewer::brewer.pal(5, "RdYlGn")
colors[3] <- "#3182bd"

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


  density_plot <- dens_df%>%
    ggplot(aes_string(x = "x", y = "y"))+
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
    theme(axis.text.y = element_text(colour = "black", size = rel(2)),
          axis.text.x = element_blank(),
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
