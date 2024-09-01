library(tidyr)
## Deseasonalized function
#' @param group_var columns used to calculate time-average values
#' @param cal_var variables need to be deseasonalized
#' @param add_var columns to identify the final output
# cal_var = c("PM2.5",metvar, emivar)
# group_var = c("Month")
# add_var = "Year"
deseason_fun <- function(data, group_var = c("Month","GridID"), cal_var, add_var = c("Date")){
  average_data <- data%>%
    group_by_at(vars(all_of(group_var)))%>%
    summarise_at(cal_var,list(de = mean))%>%ungroup()
  
  desdata <- map_dfc(cal_var, ~ left_join(data, average_data)%>%
                       transmute(!!str_c(.x,"_rm") := 
                                   !!rlang::sym(.x) - !! rlang::sym(str_c(.x, "_de"))))%>%
    bind_cols(data%>%dplyr::select(all_of(c(group_var,add_var))), .)
}


##======Main calculation function======
##  Calculate SHAP value for time series data
#' @param data observed time series data
#' @param fm formula used in shap calculation
#' @param tvar time variable
#' @param obsvar observation variable
#' @param group whether to group emission and meteorological variable
# data <- Monthly_observed_data
# sectoral_emi <- emi_ratio_by_sector_join
# tvar = c("Year","Month")
# fm = NULL
# mtry = "d3"
# obsvar = "PM25"
EMI_SHAP_ts <- function(data, sectoral_emi,
                        tvar = c("Year","Month"), obsvar = "PM25", fm = NULL,
                        mtry = "d3", sector, speciesvar, metvar){
  if(is.null(fm)){fm <- as.formula(paste0(obsvar," ~ ."))}
  emivar <- apply(expand.grid(speciesvar, sector), 1, paste, collapse="_")
  
  X <- data%>%dplyr::select(-all_of(tvar), -obsvar)
  Y <- data%>%dplyr::select(all_of(tvar), obsvar)
  
  if(mtry == "d3"){
    mtry_num <- length(X)/3
  }else{
    mtry_num <- NULL
  }
  
  rf_model <- data%>%dplyr::select(-all_of(tvar))%>%
    ranger::ranger(fm, data = ., mtry = mtry_num, num.threads = 1)
  
  shap <-  ranger.unify(rf_model, as.matrix(X))%>%
    treeshap::treeshap(X)
  
  shap_df <- shap$shaps%>%
    rename_all(~ gsub("_rm","",.x))
  
  ## calculate the change in shap values for individual species
  species_impact <- shap_df%>%
    dplyr::select(all_of(starts_with(paste0(speciesvar,"_"))))%>%
    cbind(Y%>%dplyr::select(all_of(tvar)),.)%>%
    pivot_longer(!c("Year","Month"),values_to = "SHAP")%>%
    tidyr::extract("name","Species",regex = "(.*)_",remove = FALSE)%>%
    tidyr::extract("name","GridID",regex = "_(.*)",remove = FALSE)%>%
    dplyr::select(-name)%>%
    arrange(GridID,Species,Month,Year)%>%
    group_by(GridID,Species,Month)%>%
    mutate(SHAP = SHAP - dplyr::lag(SHAP),
           SHAP = ifelse(is.na(SHAP),0,SHAP),
           GridID = as.numeric(GridID))%>%
    ungroup()

  ## assign the contribution of species for different sectors
  emi_impact <- sectoral_emi%>%
    left_join(species_impact)%>%
    filter(Year!=2000)%>%
    mutate(SHAP = SHAP*value)%>%
    # arrange(GridID,Species,)%>%
    group_by(GridID,Sector,Species,Year)%>%
    summarise(SHAP = mean(SHAP))%>%
    ungroup()%>%
    group_by(GridID,Sector,Species)%>%
    mutate(SHAP = cumsum(SHAP))%>%
    ungroup()%>%
    pivot_wider(names_from = c("Species","Sector","GridID"), values_from = "SHAP")

  met_impact <- shap_df%>%dplyr::select(starts_with(paste0(metvar,"_")))%>%
    cbind(Y%>%dplyr::select(Year),.)%>%
    group_by(Year)%>%
    summarise_all(mean)%>%
    ungroup()%>%
    mutate_at(vars(-"Year"),~.x - dplyr::lag(.x))%>%
    mutate_at(vars(-"Year"),~ifelse(is.na(.x),0,.x)%>%cumsum())
  
  for(.x in metvar){
    met_impact[paste0(.x,"_SHAP")] <- met_impact%>%dplyr::select(starts_with(paste0(.x,"_")))%>%rowSums()
  }
  
  for(.x in emivar){
    emi_impact[paste0(.x,"_SHAP")] <- emi_impact%>%dplyr::select(starts_with(paste0(.x,"_")))%>%rowSums()
    # shap_df[paste0(.x,"_local")] <- shap_df%>%dplyr::select(all_of(paste0(.x,"_rm")))%>%rowSums()
    # if(deseason){
    #   emi_impact[paste0(.x,"_local_SHAP")] <- emi_impact%>%dplyr::select(all_of(paste0(.x,"_0_rm")))%>%rowSums()
    # }else{
    #   emi_impact[paste0(.x,"_local_SHAP")] <- emi_impact%>%dplyr::select(all_of(paste0(.x,"_0")))%>%rowSums()
    # }
    # shap_df <- shap_df%>%dplyr::select(all_of(c(paste0(metvar,"_SHAP"),paste0(emivar,"_SHAP"))))
  }
  met_impact <- met_impact%>%dplyr::select(Year,all_of(c(paste0(metvar,"_SHAP"))))
  
  emi_impact <- emi_impact%>%dplyr::select(Year,all_of(c(paste0(emivar,"_SHAP"))))
  
  shap_rejoin_df <- met_impact%>%
    left_join(emi_impact)
  
  baseline <- predict(rf_model, X)$predictions%>%mean
  resid <- Y[,obsvar] - (baseline + rowSums(shap_df))

  shap_ts <- Y%>%
    mutate(baseline = baseline, resid = deframe(resid), OOB_R2 = rf_model$r.squared)%>%
    group_by(Year)%>%
    summarise_all(mean)%>%
    ungroup()%>%
    dplyr::select(-Month)%>%
    left_join(shap_rejoin_df)

  return(shap_ts)
}


## sum up the contributions of species
#' @param data data frame
species_row_sum <- function(data,
                            speciesvar = c("NOx","SO2","NMVOC","NH3","BC","OC","PM25"),
                            sector = c("agriculture","industry","power","residential","transportation"),
                            sector_ratio = FALSE){
  for(x in speciesvar){
  data <- data%>%
    mutate(!!x := !!rlang::parse_expr(paste(paste0(x,"_",sector),collapse = "+")))
  if(sector_ratio){
    data <- data%>%
      mutate_at(paste0(x,"_",sector), ~.x/!!rlang::sym(x),1)%>%
      ## prevent NA/Nan values
      mutate_at(paste0(x,"_",sector), ~ifelse(is.na(.x)|is.nan(.x),0,.x))
    }else{
    data <- data%>%
      mutate_at(paste0(x,"_",sector), ~ifelse(is.na(.x)|is.nan(.x),0,.x))
    }
  # lag_x <- dplyr::lag(data[[x]],1)

  }
  return(data)
}

## calculate the sectoral species ratio change
#' @param data
#' @param group_bar: Already-combi
species_ratio_change <- function(data, group_var,
                                 speciesvar = c("NOx","SO2","NMVOC","NH3","BC","OC","PM25")){
  
  emi_ratio <- emi_data%>%
    arrange(Sector,Month,Year)%>%
    group_by(Sector,Month)%>%
    mutate_at(all_of(speciesvar),~.x-dplyr::lag(.x,1))%>%
    ungroup()
  
  emi_ratio_all <- emi_ratio%>%
    group_by(Year,Month)%>%
    summarise_at(all_of(speciesvar),sum)%>%
    ungroup()%>%
    rename_at(all_of(speciesvar),~paste0(.x,"_all"))
  ## constrain the maximum absolute value of sectoral emission change ratio no larger than ratio_limit
  emi_ratio_by_sector <- speciesvar%>%
    map_dfc(~ emi_ratio%>%
              left_join(emi_ratio_all)%>%
              group_by(Year,Month)%>%
              mutate(ratio = !!rlang::sym(.x)/!!rlang::sym(paste0(.x,"_all")),
                     ratio_max = max(abs(ratio)),
                     ratio_flag = ifelse(ratio_max>ratio_limit,1,0),
                     ratio_constrain = ifelse(ratio_max>ratio_limit,ratio_limit/ratio_max,1),
                     ratio = ifelse(ratio_flag == 0,ratio,ifelse(abs(ratio) == ratio_max,
                                                                 (ratio-1)*ratio_constrain+1,ratio*ratio_constrain))
              )%>%
              ungroup()%>%
              transmute(!!.x := ratio)%>%
              mutate_at(.x, ~ifelse(is.na(.x)|is.nan(.x),0,.x))
    )%>%
    cbind(emi_ratio%>%dplyr::select(Year,Month,Sector),.)%>%
    mutate(GridID = 0, .before = 1)%>%
    suppressMessages()
  
}

## calculate the contributions of emissions and meteorological on gridded PM2.5
#' @param gname GridID
#' @param data: Already-combined PM2.5 and meteorological data
#' @param emi: Emission data (Don't need to be detrended)
#' @param regional_emi: The regional emissions data. It should be destrend for current version
#' @param regional_met: The regional meteorological data. It should be destrend for current version
#' @param buff_distance
#' @param deseason: Weather the variable should be deseasonalized
# data <- PM_gts
# emi <- emi_meic_China_sel
# regional_emi <- regional_emi_meic_ceds_China
# regional_met <- regional_met_China
# by_species = TRUE
# buff_distance = 250e3
# ratio_limit = 10
# sector = c("agriculture","industry","power","residential","transportation")
# speciesvar = speciesvar_meic
# speciesvar = c("NOx","SO2","NMVOC","NH3","BC","OC","PM25")
Grid_decomp <- function(data, coords, emi, regional_emi, regional_met,
                        metvar = c("u10","v10","t2m","blh","tcc","sp","e","lsp"),
                        speciesvar = c("NOx","SO2","NMVOC","NH3","BC","OC","PM25"),
                        sector = c("agriculture","industry","power","residential","transportation"),
                        deseason = FALSE, annual_ratio = FALSE, ratio_limit = 10,
                        buff_distance = 250e3){
  ## Extract emission data and join with observation data
  emivar <- apply(expand.grid(speciesvar, sector), 1, paste, collapse="_")

  ## find the closest grid of emission data
  emi_coords_filter <- emi%>%
    dplyr::distinct(X_Lon,Y_Lat)%>%
    filter(X_Lon <= coords[1]+0.2, X_Lon >= coords[1]-0.2,
           Y_Lat <= coords[2]+0.2, Y_Lat >= coords[2]-0.2)

  if(nrow(emi_coords_filter) > 0){
    emi_coords <- geosphere::distGeo(emi_coords_filter, coords)%>%
      which.min()%>%emi_coords_filter[.,]%>%as.numeric()

    emi_data <- emi%>%
      filter(X_Lon == emi_coords[1], Y_Lat == emi_coords[2])%>%
      dplyr::select(Year,Month,Sector,all_of(speciesvar))
      # pivot_wider(id_cols = c("Year","Month"),
      #             names_from = c("Sector"),
      #             values_from = speciesvar)%>%
      # species_row_sum(sector_ratio = FALSE)
  
    emi_by_species <- emi_data%>%
      group_by(Year,Month)%>%
      summarise_at(all_of(speciesvar),sum)%>%
      ungroup()%>%
      rename_at(all_of(speciesvar),~paste0(.x,"_0"))
    
    emi_ratio <- emi_data%>%
      arrange(Sector,Month,Year)%>%
      group_by(Sector,Month)%>%
      mutate_at(all_of(speciesvar),~.x-dplyr::lag(.x,1))%>%
      ungroup()
    
    emi_ratio_all <- emi_ratio%>%
      group_by(Year,Month)%>%
      summarise_at(all_of(speciesvar),sum)%>%
      ungroup()%>%
      rename_at(all_of(speciesvar),~paste0(.x,"_all"))
    ## constrain the maximum absolute value of sectoral emission change ratio no larger than ratio_limit
    emi_ratio_by_sector <- speciesvar%>%
      map_dfc(~ emi_ratio%>%
                left_join(emi_ratio_all)%>%
                group_by(Year,Month)%>%
                mutate(ratio = !!rlang::sym(.x)/!!rlang::sym(paste0(.x,"_all")),
                       ratio_max = max(abs(ratio)),
                       ratio_flag = ifelse(ratio_max>ratio_limit,1,0),
                       ratio_constrain = ifelse(ratio_max>ratio_limit,ratio_limit/ratio_max,1),
                       ratio = ifelse(ratio_flag == 0,ratio,ifelse(abs(ratio) == ratio_max,
                                                                   (ratio-1)*ratio_constrain+1,ratio*ratio_constrain))
                )%>%
                ungroup()%>%
                transmute(!!.x := ratio)%>%
                mutate_at(.x, ~ifelse(is.na(.x)|is.nan(.x),0,.x))
      )%>%
      cbind(emi_ratio%>%dplyr::select(Year,Month,Sector),.)%>%
      mutate(GridID = 0, .before = 1)%>%
      suppressMessages()
    emi_ratio_by_sector
    ## extract regional meteorological features within 500km buffer circle
    regional_met_coords_filter <- regional_met%>%distinct(X_Lon,Y_Lat)
    
    nest_met <- ((geosphere::distGeo(regional_met_coords_filter, coords))<buff_distance)%>%
      which()%>%regional_met_coords_filter[.,]%>%
      left_join(regional_met)%>%
      dplyr::select(-X_Lon,-Y_Lat)%>%
      # pivot_wider(id_cols = c("Year","Month"),names_from = "GridID",
      #             values_from = paste0(metvar,"_rm"))%>%
      pivot_wider(id_cols = c("Year","Month"),names_from = "GridID",
                  values_from = paste0(metvar))%>%
      suppressMessages()
    
    
    ## extract regional emissions
    regional_emi_coords_filter <- regional_emi%>%distinct(X_Lon,Y_Lat)
    
    nest_emi <- ((geosphere::distGeo(regional_emi_coords_filter, coords))<buff_distance)%>%
      which()%>%regional_emi_coords_filter[.,]%>%
      left_join(regional_emi)%>%
      dplyr::select(-X_Lon,-Y_Lat)

    nest_emi_by_species <- nest_emi%>%
      group_by(Year,Month,GridID)%>%
      summarise_at(all_of(speciesvar),sum)%>%
      ungroup()%>%
      pivot_wider(id_cols = c("Year","Month"),names_from = c("GridID"),
                  values_from = paste0(speciesvar))%>%
      suppressMessages()

    nest_emi_ratio <- nest_emi%>%
      arrange(GridID,Sector,Month,Year)%>%
      group_by(GridID,Sector,Month)%>%
      mutate_at(all_of(speciesvar),~.x-dplyr::lag(.x,1))%>%
      ungroup()
  
    nest_emi_ratio_all <- nest_emi_ratio%>%
      group_by(GridID,Year,Month)%>%
      summarise_at(all_of(speciesvar),sum)%>%
      ungroup()%>%
      rename_at(all_of(speciesvar),~paste0(.x,"_all"))
  
    nest_emi_ratio_by_sector <- speciesvar%>%
      map_dfc(~ nest_emi_ratio%>%
                left_join(nest_emi_ratio_all)%>%
                group_by(Year,Month,GridID)%>%
                mutate(ratio = !!rlang::sym(.x)/!!rlang::sym(paste0(.x,"_all")),
                       ratio_max = max(abs(ratio)),
                       ratio_flag = ifelse(ratio_max>ratio_limit,1,0),
                       ratio_constrain = ifelse(ratio_max>ratio_limit,ratio_limit/ratio_max,1),
                       ratio = ifelse(ratio_flag == 0,ratio,ifelse(abs(ratio) == ratio_max,
                                                                   (ratio-1)*ratio_constrain+1,ratio*ratio_constrain))
                )%>%
                ungroup()%>%
                transmute(!!.x := ratio)%>%
                mutate_at(.x, ~ifelse(is.na(.x)|is.nan(.x),0,.x))
              )%>%
      cbind(nest_emi_ratio%>%dplyr::select(GridID,Year,Month,Sector),.)%>%
      suppressMessages()

    # nest_emi_ratio_by_sector%>%
    #   group_by(GridID,Year,Month)%>%
    #   summarise_at(all_of(speciesvar),sum)%>%summary
    
    ## Emission ratio by sector: used for linear decomposition of SHAP values
    emi_ratio_by_sector_join <- full_join(emi_ratio_by_sector,nest_emi_ratio_by_sector)%>%
      pivot_longer(!c("GridID","Year","Sector","Month"),names_to = "Species", values_to = "value")%>%
      suppressMessages()
    
    # emi_ratio_by_sector_join%>%
    #   group_by(Year,Month,GridID,Species)%>%
    #   summarise_at("value",sum)%>%
    #   ungroup%>%summary

    
    Monthly_observed_data <- list(data, emi_by_species, nest_met, nest_emi_by_species)%>%
      reduce(left_join)%>%suppressMessages()
  
    PM25_ave <- mean(data$PM25)
    
    ## Join Monthly observed data
    if(deseason){
      cal_var <- names(Monthly_observed_data)%>%
        setdiff(c("Year","Month"))
     
      Monthly_deseason_data <- deseason_fun(Monthly_observed_data,
                                            group_var = c("Month"),
                                            cal_var = cal_var, add_var = c("Year"))%>%
        suppressMessages() 
    
      ## Calculate SHAP values
      SHAP_pre <- EMI_SHAP_ts(Monthly_deseason_data, sectoral_emi = emi_ratio_by_sector_join,
                              tvar = c("Year","Month"), obsvar = "PM25_rm", 
                              fm = NULL, mtry = "d3",
                              sector = sector, speciesvar = speciesvar, metvar = metvar)%>%
        mutate(PM25 = PM25_rm + PM25_ave, .after = "Year")%>%
        suppressMessages()

      # fm = paste(names(SHAP_pre)[7:length(SHAP_pre)],collapse = "+")
      # SHAP_pre%>%
      #   mutate(t := !!rlang::parse_expr(fm)+resid,.after = "Year")%>%
      #   mutate_at(c("t","PM25_rm"),~.x-dplyr::lag(.x))
    }else{
      ## Calculate SHAP values
      SHAP_pre <- EMI_SHAP_ts(Monthly_observed_data, sectoral_emi = emi_ratio_by_sector_join,
                              tvar = c("Year","Month"), obsvar = "PM25", 
                              fm = NULL, mtry = "d3",
                              sector = sector, speciesvar = speciesvar, metvar = metvar)%>%
        suppressMessages() 
    }
  
  }else{
    # SHAP_pre <- data.frame(PM25_rm = NA)
    SHAP_pre <- data.frame(Year = NA)
  }
  return(SHAP_pre)
}
