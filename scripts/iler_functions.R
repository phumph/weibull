# PTH 2020-JAN-29

# -------------------- #
# function definitions #
# -------------------- #

do_replicate_dates <- function(df) {
  stopifnot("flowers" %in% names(df))
  idx <- row.names(df)
  rep_vec <- df$flowers
  idx_rep <- rep(idx, rep_vec)
  df_rep <- df[idx_rep, ]
  stopifnot(dim(df_rep)[1] == sum(df$flowers))
  return(df_rep)
}


preprocess_data <- function(dat, min_years = 19, replicate_dates = FALSE) {

  stopifnot(all(c("doy", "species", "year", "plot") %in% names(dat)))
  
  # ensure columns are correct data type
  dat$plot    <- as.character(dat$plot)
  dat$species <- as.character(dat$species)
  dat$doy     <- as.numeric(as.character(dat$doy))
  dat$flowers <- as.numeric(dat$flowers)
  dat$year    <- as.numeric(as.character(dat$year))
  dat$habitat <- as.character(dat$habitat)
  
  # determine species present with sufficient numbers of years' worth of data to include
  dat %>%
    dplyr::group_by(species) %>%
    dplyr::summarise(n_years = length(unique(year)), .groups = "drop") %>%
    dplyr::filter(n_years >= min_years) %>%
    dplyr::select(species) ->
    species_to_include

  # filter data to species present in >= min_years
  # AND ensure all records have > 0 flowers
  # AND ensure flowers is int
  dat %>%
    dplyr::filter(species %in% species_to_include$species,
                  flowers > 0) %>%
    dplyr::mutate(flowers = round(flowers)) ->
    dat_filt
  
  # ensure each row is unique
  dat_filt <- unique(dat_filt)
  
  # date-replication step
  # each row is replicated as many times as there were counts of flowers
  # for a given observation day for a species-plot-year combination.
  # Note that hard-coded is that col index 4 equals 'flowers'.
  if (replicate_dates == TRUE) {
    dat_filt <- do_replicate_dates(dat_filt)
  }

  # re-order data.frame
  dat_filt <-
    dat_filt %>%
    dplyr::arrange(species, plot, year, doy)

  return(dat_filt)
}


calc_weibull_pearse <- function(dat, by_plot = FALSE, k_min = 10, K = NULL) {

  if (by_plot == TRUE) {
    factor_list <- list(dat$species, dat$year, dat$plot)
    names(factor_list) <- c("species", "year", "plot")
  } else if (by_plot == FALSE) {
    factor_list <- list(dat$species, dat$year)
    names(factor_list) <- c("species", "year")
  } else {
    stop("Value for by_plot needs to be T/F.", call. = FALSE)
  }

  theta <- as.numeric(tapply(dat$doy,
                             factor_list,
                             function(x) {
                               if(length(x) > k_min) {
                                 weib.limit(x, k = K, upper = FALSE)
                               } else {
                                 NA
                               }
                             }))
  
  theta.end <- as.numeric(tapply(dat$doy,
                                 factor_list,
                                 function(x) {
                                   if(length(x) > k_min) {
                                     weib.limit(x, k = K, upper = TRUE)
                                   } else {
                                     NA
                                   }
                                 }))

  # first day of flowers observed
  min <- as.numeric(
    tapply(dat$doy,
           factor_list,
           function(x) if(length(x) > k_min) min(x) else NA))
  
  # last day of flowers observed
  max <- as.numeric(
    tapply(dat$doy,
           factor_list,
           function(x) if(length(x) > k_min) max(x) else NA))
  
  # mean day of year for flowers observed
  mean <- as.numeric(
    tapply(dat$doy,
           factor_list,
           function(x) if(length(x) > k_min) mean(x, na.rm = TRUE) else NA))

  # median day of year for flowers observed
  median <- as.numeric(
    tapply(dat$doy,
           factor_list,
           function(x) if(length(x) > k_min) median(x, na.rm = TRUE) else NA))

  
  if (!"flowers" %in% names(dat)) {
    # obsolete and incorrect abundance function from Pearse et al.
    abundance <- as.numeric(
      tapply(dat$doy,
             factor_list,
             function(x) if(length(x) > k_min) max(table(x)) else NA))
  } else {
    
    factor_list_2 <- factor_list
    factor_list_2[["doy"]] <- dat$doy
    
    dat %>%
      dplyr::select(flowers) %>%
      aggregate(by = factor_list_2, FUN = sum) ->
      dat2
    
    # maximum observed flower abundance for a single day
    abundance <- as.numeric(
      tapply(dat2$flowers,
             list(dat2$species, dat2$year),
             function(x) if(length(x) > k_min) max(x) else NA))  
  }
  
  species <- as.character(
    tapply(dat$species,
           factor_list, unique))

  year <- as.numeric(
    tapply(dat$year,
           factor_list, unique))

  results <- data.frame(
      theta_onset = theta,
      theta_end = theta.end,
      min,
      max,
      mean,
      species,
      year,
      median,
      abundance)

  return(results)
}

# write script for producing weibull estimates
calc_weibull_pearse_obs <- function(dat, by_plot = FALSE, k_min = 10) {

  if (by_plot == TRUE){
    factor_list <- list(dat$species, dat$plot, dat$year)
  } else if (by_plot == FALSE) {
    factor_list <- list(dat$species, dat$year)
  } else {
    stop("Value for by_plot needs to be T/F.", call. = FALSE)
  }

  theta <- as.numeric(tapply(dat$doy,
                             factor_list,
                             function(x) {
                               if(length(x) > k_min) {
                                 theta.hat(x, k = 30, late = FALSE)
                               } else {
                                 NA
                               }
                             }))

  theta.end <- as.numeric(tapply(dat$doy,
                                 factor_list,
                                 function(x) {
                                   if(length(x) > k_min) {
                                     theta.hat(x, k = 30, late = TRUE)
                                   } else {
                                     NA
                                   }
                                 }))

  min <- as.numeric(
    tapply(dat$doy,
           factor_list,
           function(x) if(length(x) > k_min) min(x) else NA))

  max <- as.numeric(
    tapply(dat$doy,
           factor_list,
           function(x) if(length(x) > k_min) max(x) else NA))

  mean <- as.numeric(
    tapply(dat$doy,
           factor_list,
           function(x) if(length(x) > k_min) mean(x, na.rm = TRUE) else NA))

  median <- as.numeric(
    tapply(dat$doy,
           factor_list,
           function(x) if(length(x) > k_min) median(x, na.rm = TRUE) else NA))

  abundance <- as.numeric(
    tapply(dat$doy,
           factor_list,
           function(x) if(length(x) > k_min) max(table(x)) else NA))

  species <- as.character(
    tapply(dat$species,
           factor_list, unique))

  year <- as.numeric(
    tapply(dat$year,
           factor_list, unique))

  results <- na.omit(data.frame(theta, theta.end, min, max, mean, species, year, median, abundance))

  return(results)
}

# function to transform year column, as in Pearse et al. (2017)
rescale_year <- function(x, year_col = 'year', uniques = TRUE, fill_missing = TRUE) {

  the_years <- x[[year_col]]

  # transform years vector if conditions are true
  if (uniques == TRUE) {
    the_years <- unique(x[[year_col]])
    if (fill_missing == TRUE) {
      the_years = seq(min(the_years), max(the_years), by = 1)
    }
  }

  # re-scale
  year_z <- (the_years - mean(the_years, na.rm = T)) / sd(the_years, na.rm = T)
  the_years_df <- data.frame(the_years, year_z) %>% unique()

  names(the_years_df)[1] <- year_col

  x %>%
    dplyr::left_join(the_years_df, by = 'year') ->
    x2

  return(x2)
}


coef_extract <- function(model, time_term, model_name, gsub_str = NA) {

  coefs <-
    model %>%
    broom::tidy() %>%
    dplyr::filter(grepl(paste0(':',time_term), term)) %>%
    dplyr::mutate(model = model_name,
                  term = sapply(term, function(x) gsub(paste0(':',time_term), '', x)))

  if (!is.null(gsub_str)) {
    coefs$term <- sapply(coefs$term, function(x) {
      gsub(gsub_str, '', x)
    })
  }

  names(coefs)[!grepl('term|model', names(coefs))] <- paste0(names(coefs)[!grepl('term|model', names(coefs))],
                                                             '_',
                                                             model_name)
  names(coefs)[1] <- 'species'
  return(coefs)
}


run_lms <- function(x) {

  onset.model.all <- lm(theta_onset ~ species * year_z + log(abundance), data = x)
  #peak.model.all  <- lm(em50        ~ species * year_z + log(abundance), data = x)
  peak.model.all  <- lm(mean        ~ species * year_z + log(abundance), data = x)
  end.model.all   <- lm(theta_end   ~ species * year_z + log(abundance), data = x)
  

  # extract the coefficients
  onset_coefs <-
    onset.model.all %>%
    coef_extract(time_term  = 'year_z',
                 model_name = 'onset',
                 gsub_str   = 'species')

  peak_coefs <-
    peak.model.all %>%
    coef_extract(time_term  = 'year_z',
                 model_name = 'mean',
                 gsub_str   = 'species')

  end_coefs <-
    end.model.all %>%
    coef_extract(time_term  = 'year_z',
                 model_name = 'end',
                 gsub_str   = 'species')

  all_coefs <-
    dplyr::full_join(
      dplyr::select(onset_coefs, -model),
      dplyr::select(peak_coefs, -model),
      by = 'species') %>%
    dplyr::full_join(
      dplyr::select(end_coefs, -model),
      by = 'species')

  return(all_coefs)

}


run_deming <- function(x) {

  d.fvp <- deming(estimate_mean ~ estimate_onset,
                  xstd = std.error_onset,
                  ystd = std.error_mean,
                  data = x)

  # first vs. end
  d.fve <- deming(estimate_end ~ estimate_onset,
                  xstd = std.error_onset,
                  ystd = std.error_end,
                  data = x)

  # peak vs. end
  d.pve <- deming(estimate_end ~ estimate_mean,
                  xstd = std.error_mean,
                  ystd = std.error_end,
                  data = x)

  # compile and export confidence limits on slope and intercept estimates
  list(d.fvp,
       d.fve,
       d.pve) -> dl1

  names(dl1) <- c('first_v_peak',
                  'first_v_end',
                  'peak_v_end')

  dl1.coefs <- lapply(dl1, function(x) data.frame(x$coefficients,x$ci))

  for (i in seq_along(dl1.coefs)) {
    dl1.coefs[[i]]$model <- names(dl1)[i]
  }

  dl1.coefs <-
    dl1.coefs %>%
    do.call(rbind, .)

  names(dl1.coefs)[1] <- c('estimate')
  dl1.coefs$term <- 'slope'
  dl1.coefs$term[grep('Intercept', row.names(dl1.coefs))] <- 'intercept'

  return(dl1.coefs)
}


plot_coefs <- function(x, y, plot_title) {

  x %>%
    ggplot(aes(x = estimate_onset,
               y = estimate_mean,
               xmin = estimate_onset - std.error_onset,
               xmax = estimate_onset + std.error_onset,
               ymin = estimate_mean - std.error_mean,
               ymax = estimate_mean + std.error_mean)) +
    geom_errorbar(alpha = 0.5, lwd = 0.25, col = 'dodgerblue') +
    geom_errorbarh(alpha = 0.5, lwd = 0.25, col = 'dodgerblue') +
    geom_point(col = 'dodgerblue') +
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(size = 8, face = 'bold')) +
    geom_abline(intercept = 0, slope = 1, lty = 2, alpha = 0.5, col = 'gray40') +
    geom_abline(intercept = y$estimate[y$term == 'intercept' & y$model == 'first_v_peak'],
                slope = y$estimate[y$term == 'slope' & y$model == 'first_v_peak'],
                lty = 1,
                alpha = 0.5,
                col = 'dodgerblue',
                lwd = 1) +
    coord_cartesian(xlim = c(-30,15), ylim = c(-30,15)) +
    ggtitle(paste0(plot_title)) ->
    first_v_peak_plot

  x %>%
    ggplot(aes(x = estimate_onset,
               y = estimate_end,
               xmin = estimate_onset - std.error_onset,
               xmax = estimate_onset + std.error_onset,
               ymin = estimate_end - std.error_end,
               ymax = estimate_end + std.error_end)) +
    geom_errorbar(alpha = 0.5, lwd = 0.25, col = 'red') +
    geom_errorbarh(alpha = 0.5, lwd = 0.25, col = 'red') +
    geom_point(col = 'red') +
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(size = 8, face = 'bold')) +
    geom_abline(intercept = 0, slope = 1, lty = 2, alpha = 0.5, col = 'gray40') +
    geom_abline(intercept = y$estimate[y$term == 'intercept' & y$model == 'first_v_end'],
                slope = y$estimate[y$term == 'slope' & y$model == 'first_v_end'],
                lty = 1,
                alpha = 0.5,
                col = 'red',
                lwd = 1) +
    coord_cartesian(xlim = c(-30,15), ylim = c(-30,15)) +
    ggtitle(paste0(plot_title)) ->
    first_v_end_plot

  ggpubr::ggarrange(plotlist = list(first_v_peak_plot, first_v_end_plot),
                    nrow = 2,
                    ncol = 1) -> deming_plot

  return(deming_plot)
}

plot_rescaled_year <- function(x, the_title) {

  years <- seq(min(x$year), max(x$year),by = 1)
  years_z <- (years - mean(years)) / sd(years)
  y <- data.frame(year = years, year_z_correct = years_z)

  x %>%
    left_join(y, by = 'year') %>%
    dplyr::select(year, year_z, year_z_correct) %>%
    unique() %>%
    ggplot(aes(x = year_z_correct, y = year_z)) +
    geom_point(alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0) +
    theme_bw() +
    coord_cartesian(xlim = c(-2,2), ylim = c(-2,2)) +
    ggtitle(paste0(the_title)) +
    theme(plot.title = element_text(size = 8, face = 'bold'))
}


filter_species_list <- function(x, spp_to_remove) {

  x <- x[!x$species %in% spp_to_remove, ]

  return(x)
}
