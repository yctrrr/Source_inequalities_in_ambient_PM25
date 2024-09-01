##  Function for inequality calculation  ##
##  Version 1 # Yanchuan Shao 5th May 2024  ##
library(tidyverse)
library(ggplot2)

## Local group summary
#'@param group_main the inequality category that should be analyzed (e.g. income)
#'@param group_add the column that indicates additional differences and how the overall average is defined (e.g. time)
#'@param pop_name the column name of population
#'@param var variable of interest
#'@param prefix should the variable name be prefixed in the output?
local_popw_summarise <- function(data, group_main = "GDPQ", group_add = NULL,
                                 pop_name = "Pop", var = "EMI", prefix = TRUE, ...){
  if(isTRUE(prefix)){
    popw_name <- paste0(var,"_Popw")
    ave_name <- paste0(var,"_Ave")
    sd_name <- paste0(var,"_SD")
  }else{
    popw_name <- paste0("Popw")
    ave_name <- paste0("Ave")
    sd_name <- paste0("SD")
  }

  average <- data%>%
    group_by_at(all_of(group_add))%>%
    summarise(!!rlang::sym(popw_name) := weighted.mean(!!rlang::sym(var), !!rlang::sym(pop_name)),
              !!rlang::sym(ave_name) := mean(!!rlang::sym(var)),
              !!rlang::sym(sd_name) := sd(!!rlang::sym(var)),
              # !!rlang::sym(paste0(var,"_Gini")) := DescTools::Gini(!!rlang::sym(var),Pop,unbiased = FALSE),
              Pop = sum(!!rlang::sym(pop_name)),
              ...)%>%ungroup()%>%
    as.data.frame()%>%suppressMessages()
  
  ## the main group category will be renamed as 0
  for(i in c(group_main)){
    average <- average%>%
      mutate(!!rlang::sym(i) := 0)
  }

  data%>%
    group_by_at(all_of(c(group_main, group_add)))%>%
    summarise(!!rlang::sym(popw_name) := weighted.mean(!!rlang::sym(var), !!rlang::sym(pop_name), na.rm = TRUE),
              !!rlang::sym(ave_name) := mean(!!rlang::sym(var)),
              !!rlang::sym(sd_name) := sd(!!rlang::sym(var)),
              # !!rlang::sym(paste0(var,"_Gini")) := DescTools::Gini(!!rlang::sym(var),Pop,unbiased = FALSE),
              Pop = sum(!!rlang::sym(pop_name)),
              )%>%
    ungroup()%>%
    as.data.frame()%>%
    ##  Prevent NaN return in Population-weighted values
    mutate(!!rlang::sym(popw_name) := ifelse(Pop == 0, !!rlang::sym(ave_name), !!rlang::sym(popw_name)))%>%
    full_join(average)%>%suppressMessages()
}


## Generate group category for data frame by variable quantile
GroupQ_mutate <- function(data, group_name, group_quantile = seq(0,1,.2)){
  data%>%
    mutate(!!sym(paste0(group_name,"_group")) := 
             cut(!!sym(group_name), quantile(!!sym(group_name), group_quantile),
                 include.lowest = TRUE, labels = 1:(length(group_quantile)-1))%>%
             as.numeric()
    )
}

## Rename the GDP category
GDPQ_rename <- function(data, var = "GDPQ"){
  data%>%
    mutate_at(all_of(var), ~ ifelse(.x == 1, "Very Poor",
                                    ifelse(.x == 2, "Poor",
                                           ifelse(.x == 3, "Middle",
                                                  ifelse(.x == 4, "Rich", "Very Rich")))))%>%
    mutate_at(all_of(var), ~ factor(.x, levels = c("Very Poor","Poor","Middle","Rich", "Very Rich")))
}

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

## Visualize the inequality by different GDP groups
Source_ineq_vis <- function(data, x = "Sector", y = "Popw_ineq", fill = "Popw",
                            facet = "GDPQ", limits = NULL, 
                            fill_title = expression(paste("Population-weighted PM"[2.5]," (μg/m"^"3",")")),
                            title = "Exposure Disparity (%)"){
  ggplot(data, aes_string(x = x, y = y))+
    geom_bar(aes_string(fill = fill), stat = "identity", color = "black")+
    facet_wrap(facet, nrow = 1)+
    coord_flip()+
    theme_bw()+
    labs(x = NULL, y = title)+
    scale_x_discrete(limits = rev)+
    scale_fill_distiller(name = fill_title,
                         limits = limits,
                         palette = "YlOrRd", direction = 1,
                         guide = guide_colorbar(frame.colour = "black",
                                                frame.linewidth = .5,
                                                title.position = "left", ticks = F,
                                                title.vjust = 1))+
    theme(panel.background = element_rect(size = 1,fill='white', colour='black'),
          panel.border = element_rect(size = 1,colour = "black", fill=NA))+
    theme(axis.text = element_text(colour = "black",size = rel(2.5)), ##axis text size
          axis.title = element_text(colour = "black",size = rel(3)),
          plot.margin = unit(c(5,10,5,5),"mm"),
          strip.text.x = element_text(size = 35, color = "black",margin = margin(.1, -.1, .1, -.1, "cm")),
          strip.background = element_blank())+ ##axis title size
    theme(legend.key.width=unit(4.5, "cm"),
          legend.key.height=unit(1, "cm"),  ##legend size
          legend.position = "top",
          legend.text = element_text(colour = "black",size = rel(2.5)),
          legend.title = element_text(colour = "black",size = rel(2.5)),
          legend.ticks = element_blank(),
          legend.background = element_blank()
          # legend.box.background = element_rect(size= 0.5,colour = "black")
    ) ##legend text size
}

## Visualize the inequality along the time series
Source_ineq_ts_vis <- function(data, x = "Year", y = "Popw_ineq", breaks = c(2005,2015),
                               source = "Sector", facet = "GDPQ",
                               title = "Exposure Disparity (%)", colors = NULL){
  if(is.null(colors)){
    colors <- RColorBrewer::brewer.pal(5, "RdYlGn")
    colors[3] <- "#3182bd"
  }

  if(is.null(facet)){
    ggplot(data, aes_string(x = x, y = y, color = source))+
      geom_point(size = 2, show.legend = FALSE)+
      geom_line(linewidth = 1.6)+
      geom_abline(intercept = 0, slope = 0, linewidth = .6)+
      labs(x = NULL, y = title)+
      scale_x_continuous(breaks = breaks)+
      scale_color_manual(name = NULL, values = colors, 
                         guide = guide_legend(frame.colour = "white", frame.linewidth = .1,
                                              title.position = "top", title.hjust = 0.5,
                                              byrow = TRUE, ticks = F)
                         )+
      theme(panel.background = element_rect(size = 1, fill='white', colour='black'),
            panel.border = element_rect(size = 2,colour = "black", fill=NA),
            )+
      theme(axis.text = element_text(colour = "black",size = rel(3)), ##axis text size
            axis.title = element_text(colour = "black",size = rel(3.5)),
            plot.margin = unit(c(5,0,5,0),"mm"),
            strip.text.x = element_text(size = 35, color = "black",margin = margin(.5, -.1, .1, -.5, "cm")),
            strip.background = element_blank())+ ##axis title size
      theme(legend.key.width=unit(1, "cm"),
            legend.key.height=unit(1, "cm"),  ##legend size
            legend.spacing.y = unit(.8, 'cm'),
            legend.key.spacing.y = unit(.8, 'cm'),
            legend.key = element_blank(),
            legend.position = "right",
            legend.text = element_text(colour = "black",size = rel(3)),
            legend.title = element_text(colour = "black",size = rel(3)),
            legend.background = element_blank())
  }else{
    ggplot(data, aes_string(x = x, y = y, color = source))+
      geom_point(size = 2, show.legend = FALSE)+
      geom_line(linewidth = 1.6)+
      geom_abline(intercept = 0, slope = 0, linewidth = .6)+
      facet_wrap(facet, nrow = 1)+
      labs(x = NULL, y = title)+
      scale_x_continuous(breaks = breaks)+
      scale_color_manual(name = NULL, values = colors,
                         guide = guide_legend(frame.colour = "white", frame.linewidth = .1,
                                              title.position = "top", title.hjust = 0.5,
                                              byrow = TRUE, ticks = F)
                         )+
      theme(panel.background = element_rect(size = 1,fill='white', colour='black'),
            panel.border = element_rect(size = 2,colour = "black", fill=NA),
            )+
      theme(axis.text = element_text(colour = "black",size = rel(3)), ##axis text size
            axis.title = element_text(colour = "black",size = rel(3.5)),
            plot.margin = unit(c(5,0,5,0),"mm"),
            strip.text.x = element_text(size = 35, color = "black",margin = margin(.5, -.1, .1, -.5, "cm")),
            strip.background = element_blank())+ ##axis title size
      theme(legend.key.width=unit(1, "cm"),
            legend.key.height=unit(1, "cm"),  ##legend size
            legend.spacing.y = unit(.8, 'cm'),
            legend.key.spacing.y = unit(.8, 'cm'),
            legend.key = element_blank(),
            legend.position = "right",
            legend.text = element_text(colour = "black",size = rel(3)),
            legend.title = element_text(colour = "black",size = rel(3.5)),
            legend.background = element_blank())
  }
}

## Visualize the inequality for specific sources using box plot
# x = "Year"
# y = "Popw_ineq"
# fill = "Popw"
# facet = "GDPQ"
# title = "Exposure Disparity (%)"
# limits = c(8,20)
# fill_title = expression(paste("Population-weighted PM"[2.5]," (μg/m"^"3",")"))
Source_ineq_box_ts_vis <- function(data, x = "Year", y = "Popw_ineq", fill = "Popw",
                                   facet = "GDPQ", limits = c(8,20), breaks = c(2005,2015),
                                   fill_title = expression(paste("Population-weighted PM"[2.5]," (μg/m"^"3",")")),
                                   title = "Exposure Disparity (%)"){
  ggplot(data, aes_string(x = x, y = y, fill = fill))+
    geom_bar(aes_string(fill = fill), stat = "identity")+
    facet_wrap(paste0(facet), nrow = 1)+
    labs(x = NULL, y = title)+
    scale_x_continuous(breaks = breaks)+
    scale_fill_distiller(name = fill_title,
                         limits = limits,
                         # palette = "BuPu", 
                         palette = "Spectral",
                         direction = -1,
                         guide = guide_colorbar(frame.colour = "black",
                                                frame.linewidth = .5,
                                                title.position = "right", ticks = FALSE,
                                                title.vjust = 1))+
    theme(panel.background = element_rect(size = 1, fill='white', colour='black'),
          panel.border = element_rect(size = 1, colour = "black", fill=NA))+
    theme(axis.text = element_text(colour = "black",size = rel(3)), ##axis text size
          axis.title = element_text(colour = "black",size = rel(3.5)),
          plot.margin = unit(c(5,0,5,0),"mm"),
          strip.text.x = element_text(size = 35, color = "black",margin = margin(.1, -.1, .1, -.1, "cm")),
          strip.background = element_blank(),
          )+ ##axis title size
    theme(legend.key.width = unit(1, "cm"),
          legend.key.height = unit(3, "cm"),  ##legend size
          legend.spacing.y = unit(.8, 'cm'),
          legend.key.spacing.y = unit(.8, 'cm'),
          legend.key = element_blank(),
          legend.position = "right",
          legend.ticks = element_blank(),
          legend.text = element_text(colour = "black",size = rel(3)),
          legend.title = element_text(colour = "black",size = rel(3.5), angle = -90),
          legend.background = element_blank())
}

pred_df_extract <- function(model, data,  voi_log = FALSE, voi, voi_fm){
  # predictions
  pred <- predict(model, se.fit = TRUE)
  pred$fit%>%summary
  
  (pred$fit+model$coefficients[1])%>%head
  pred$fit%>%head
  if(isTRUE(voi_log)){
    indp_val <- log(data[, voi])
  }else{
    indp_val <- data[, voi]
  }
  
  # polish predictions
  pred_df <- tibble(pred$fit[,voi_fm], pred$se.fit[,voi_fm], indp_val)%>%
    as.data.frame()%>%
    janitor::clean_names()%>%
    dplyr::rename(beta = 1,
                  se = 2,
                  indp_var = 3)%>%
    arrange(indp_var)%>%
    mutate(lci = beta - 1.96 * se,
           uci = beta + 1.96 * se)
  return(pred_df)
}

## plot the variable relationship derived from GAM
plot_nonlin_gam <- function(model = NULL, data = NULL, voi_log = FALSE, 
                            voi = NULL, voi_fm = NULL, 
                            pred_df = NULL, custom = FALSE,
                            base_family = "", color = NULL){
  
  if(is.null(pred_df)){
    pred_df <- pred_df_extract(model, data, voi_log, voi, voi_fm)
  }

  # plot
  if(!custom){
    if(is.null(color)) color <- "#99d8c9"
    base_plot <- pred_df %>%
      ggplot(aes(x = indp_var))+
      geom_ribbon(aes(ymin = lci, ymax = uci), fill = color, alpha = 0.3)+
      geom_path(aes(y = beta), size = 1.2)+
      geom_hline(yintercept = 0, linetype = 2, size =0.5)+
      theme_bw(base_size = 20, base_family = base_family)+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(size = 1, fill = 'white', colour = 'black'),
            panel.border = element_rect(size = 1, colour = "black", fill = NA))+
      theme(axis.text = element_text(colour = "black",size = rel(2)),
            axis.title = element_text(colour = "black",size = rel(3)),
            plot.margin = margin(20, 50, 0, 0)) ##axis title size
  }else{
    base_plot <- pred_df %>%
      ggplot(aes(x = indp_var))+
      geom_hline(yintercept = 0, linetype = 2, size =0.5)+
      theme_bw(base_size = 20, base_family = base_family)+
      theme(strip.text.x = element_text(size = rel(2), color = "black",margin = margin(.1, -.1, .1, -.1, "cm")),
            strip.background = element_blank())+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(size = 1, fill = 'white', colour = 'black'),
            panel.border = element_rect(size = 1, colour = "black", fill = NA),
            panel.spacing = unit(2, "cm"))+
      theme(axis.text = element_text(colour = "black",size = rel(1.5)), ##axis text size
            axis.title = element_text(colour = "black",size = rel(2)),
            plot.margin = unit(c(5,0,5,0),"mm")) ##axis title size
  }
  plot <- base_plot+
    theme(legend.key.width=unit(1.2, "cm"),
          legend.key.height=unit(1, "cm"),  ##legend size
          legend.spacing.y = unit(.8, 'cm'),
          legend.key.spacing.y = unit(.8, 'cm'),
          legend.key = element_blank(),
          legend.position = "right",
          legend.text = element_text(colour = "black",size = rel(1.5)),
          legend.title = element_text(colour = "black",size = rel(1.5)),
          legend.background = element_blank())
  plot
}


plot_linear <- function(model, data, one_basis, voi,
                        base_family = "", color = "#99d8c9"){
  # predict
  pred <- crosspred(basis = one_basis, model = model, 
                    at = data[[voi]], cen = mean(data[[voi]]))
  
  fit <- as_tibble(pred$matfit, rownames = NA) %>%
    tibble::rownames_to_column(., "indp_var") %>%
    rename(beta = lag0)
  
  # ii Extract 95% CI  
  lci <- as_tibble(pred$matlow, rownames = NA) %>%
    tibble::rownames_to_column(., "indp_var") %>%
    rename(lci = lag0) 
  
  uci <- as_tibble(pred$mathigh, rownames = NA) %>%
    tibble::rownames_to_column(., "indp_var") %>%
    rename(uci = lag0)
  
  # iii Combine fit and se 
  table_full <- left_join(fit, lci) %>%
    left_join(., uci) %>%
    dplyr::mutate(indp_var = as.numeric(indp_var))
  
  # iv Plot
  plot <- table_full %>%
    ggplot(aes(x = indp_var)) +
    theme_minimal(base_size = 15) +
    geom_ribbon(aes(ymin = lci, ymax = uci), fill= color, alpha = 0.3) +
    geom_path(aes(y = beta), size = 0.7) +
    geom_hline(yintercept = 0, linetype = 2, size =0.5) +
    theme(plot.margin = unit(c(5,5,5,5),"mm"))
  
  plot
}
