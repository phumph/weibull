
rm(list = ls(all = TRUE))

# -------------- #
# load libraries #
# -------------- #

library(reshape2)
library(deming)
library(car)
library(devtools)
library(broom)
library(plyr)
library(lme4)
library(phest)
library(lmerTest)
library(dplyr)
library(stringr)
library(ggplot2)

# -------------------- #
# function definitions #
# -------------------- #

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


################################################################################
##                                                                            ##
##  ~~  REPEAT Pearse MODEL STRUCTURE WITH CaraDonna WEIBULL ESTIMATES  ~~    ##
##                                                                            ##
################################################################################

source(file.path("./2018-07-09_headers.R")) # from Will Pearse

data <- na.omit(read.csv("../data/rmbl_pnas_data.csv"))

## rename some headers; make all lowercase
names(data)[4] <- "flowers"
names(data) <- tolower(names(data))

# subset data to include only years after 1973
data <- subset(data, year > '1973')

## Determine whether each species has a record for each year
## sum number of years for which this is true
## and return TRUE if this number if at least 19.
## Returns FALSE if fewer than 19 years on which at leats one observation was made.
counts <- apply(
  with(data,
       tapply(doy,
              list(species, year),
              function(x) length(x) > 0)),
  1,
  function(x) sum(x, na.rm = TRUE) >= 19)

## subset initial data.frame to include only species passing this threshold.
## in this case, the only species present initially are those passing this threshold.
data <- data[data$species %in% names(counts)[counts == TRUE],]

## re-order dataset:
data <- data[order(data$species, data$plot, data$year, data$doy), ]

## ensure factors are actually stored as characters or numeric vectors
data$species <- as.character(data$species)
data$doy     <- as.numeric(as.character(data$doy))
data$year    <- as.numeric(as.character(data$year))

# set min number of observations to warrant inclusion of plant species in a given year
k_min <- 10

# estimate first and last flowering dates with Weibull estimator (requires 11 dates for each species in each year; otherwise get an NA)
onset   <- as.numeric(with(data,
                           tapply(doy,
                                  list(species, year),
                                  function(x) if(length(x) > k_min) theta.hat(x, k = 30) else NA)))

end     <- as.numeric(with(data,
                           tapply(doy,
                                  list(species, year),
                                  function(x) if(length(x) > k_min) theta.hat(x, k = 30, TRUE) else NA)))

min.doy <- as.numeric(with(data, tapply(doy, list(species, year), function(x) min(x))))
max.doy <- as.numeric(with(data, tapply(doy, list(species, year), function(x) max(x))))

# create dataframe
species <- as.character(with(data, tapply(species, list(species, year), unique)))
year    <- as.numeric(with(data, tapply(year, list(species, year), unique)))
results <- data.frame(onset, end, species, year, min.doy, max.doy)
results <- results[results$species %in% names(Filter(function(x) x > k_min, table(results$species))), ] #Amy spot checked to make sure code did what it is supposed to do

# clean up
rm(onset, end, species, year, min.doy, max.doy)

# three unreasonable end estimates changed to NA's in the file below (AcMi 1976, DuHo 1976, CoLi 1996); see Appendix 1 for details
sort.results <- results[order(-results$end), ] # remove
sort.results$end[1] <- NA
sort.results$end[2] <- NA
sort.results$end[4] <- NA
results.corrected <- sort.results

# add flower abundance data to corrected.results dataframe
sum.data <- ddply(data,
                  .(species, year, doy), 
                  summarize, 
                  sum.flws = sum(flowers)) # sum flowers across plots on each date to get peak #flowers for each species in each year

peak.abundance.data <- ddply(sum.data,
                             .(year, species),
                             summarize,
                             peak.flowers = max(sum.flws))

results.corrected <- results.corrected[order(results.corrected$year, results.corrected$species), ] # now both dataframes are sorted by year then species

sum.data <-
  data %>%
  dplyr::group_by(species, year, doy) %>%
  dplyr::summarise(sum.flws = sum(flowers)) %>%
  dplyr::group_by(year, species) %>%
  dplyr::summarise(peak.flowers = max(sum.flws))

merge.results <-
  results.corrected %>%
  dplyr::left_join(sum.data, by = c('species', 'year'))

#head(merge.results)

# Remove 8 species with fewer than 19 years of data from species list in CaraDonna et al. (2014)
# ("Maianthemum stellatum" "Gentiana parryi" "Lomatium dissectum" "Dodecatheon pulchellum"
#  "Hydrophyllum fendleri" "Descurainia richardsonii" "Pedicularis bracteosa" "Sedum rosea")

spp_to_remove <- c("Maianthemum stellatum",
                   "Gentiana parryi",
                   "Lomatium dissectum",
                   "Dodecatheon pulchellum",
                   "Hydrophyllum fendleri",
                   "Descurainia richardsonii",
                   "Pedicularis bracteosa",
                   "Sedum rosea")

species.subset.n.flws <-
  merge.results %>%
  dplyr::filter(!species %in% spp_to_remove)

# now add em50 data to this dataframe (em50 = day on which 50% of flowers were counted)
em50.wide <- read.csv("../data/EM50_RMBL.csv", header=T)
em50.long <- melt(em50.wide, id.vars ='year')
names(em50.long)[2] <- "species"
names(em50.long)[3] <- "em50"

em50.long$species <- as.character(em50.long$species)
em50.long$year <- as.numeric(em50.long$year)
em50.long$species <- gsub("\\.", " ", em50.long$species)

species.subset.n.flws$species   <- gsub("Potentilla gracilis v. pulcherrima", "Potentilla gracilis v pulcherrima", species.subset.n.flws$species)
species.subset.n.flws$species   <- gsub("Lupinus sp.", "Lupinus sp", species.subset.n.flws$species)
species.subset.n.flws$year.1    <- NULL
species.subset.n.flws$species.1 <- NULL

all.data <- merge(species.subset.n.flws, em50.long, by = c("species", "year"))

# convert year, as in Pearse et al.
all.data$orig.year <- all.data$year
all_years <- unique(all.data$year)
all_years_scaled_df <- data.frame(year = all_years,
                                  year_z = as.numeric(scale(all_years)))

# join centered/scaled year back with all.data
all.data <- dplyr::left_join(all.data, 
                             all_years_scaled_df,
                             by = 'year')

# run the models
onset.model.all <- lm(onset ~ species * year_z + log(peak.flowers), data = all.data)
peak.model.all  <- lm(em50  ~ species * year_z + log(peak.flowers), data = all.data)
end.model.all   <- lm(end   ~ species * year_z + log(peak.flowers), data = all.data)

# extract the coefficients
onset_coefs <- 
  onset.model.all %>%
  coef_extract(time_term  = 'year_z',
               model_name = 'onset',
               gsub_str   = 'species')

peak_coefs <- 
  peak.model.all %>%
  coef_extract(time_term  = 'year_z',
               model_name = 'em50',
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

### repeat Deming regressions with shifts from model with all species and flower abundance covariate ###

# first vs. peak
d.fvp <- deming(estimate_em50 ~ estimate_onset,
                xstd = std.error_onset,
                ystd = std.error_em50,
                data = all_coefs)

# first vs. end
d.fve <- deming(estimate_end ~ estimate_onset,
                xstd = std.error_onset,
                ystd = std.error_end,
                data = all_coefs)

# peak vs. end
d.pve <- deming(estimate_end ~ estimate_em50,
                xstd = std.error_em50,
                ystd = std.error_end,
                data = all_coefs)

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

write.table(dl1.coefs,
            file = '../output/deming_reg_coefs_Iler.csv',
            quote = F,
            sep = ',',
            row.names = F,
            col.names = T)

# produce plots of the Deming regression results:
all_coefs %>%
  ggplot(aes(x = estimate_onset,
             y = estimate_em50,
             xmin = estimate_onset - std.error_onset,
             xmax = estimate_onset + std.error_onset,
             ymin = estimate_em50 - std.error_em50,
             ymax = estimate_em50 + std.error_em50)) +
  geom_errorbar(alpha = 0.5, lwd = 0.25, col = 'dodgerblue') +
  geom_errorbarh(alpha = 0.5, lwd = 0.25, col = 'dodgerblue') +
  geom_point(col = 'dodgerblue') +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_abline(intercept = 0, slope = 1, lty = 2, alpha = 0.5, col = 'gray40') +
  geom_abline(intercept = dl1.coefs$estimate[dl1.coefs$term == 'intercept' & dl1.coefs$model == 'first_v_peak'],
              slope = dl1.coefs$estimate[dl1.coefs$term == 'slope' & dl1.coefs$model == 'first_v_peak'],
              lty = 1,
              alpha = 0.5,
              col = 'dodgerblue',
              lwd = 1) +
  coord_cartesian(xlim = c(-25,15), ylim = c(-25,10)) ->
  first_v_peak_plot

all_coefs %>%
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
  theme(panel.grid = element_blank()) +
  geom_abline(intercept = 0, slope = 1, lty = 2, alpha = 0.5, col = 'gray40') +
  geom_abline(intercept = dl1.coefs$estimate[dl1.coefs$term == 'intercept' & dl1.coefs$model == 'first_v_end'],
              slope = dl1.coefs$estimate[dl1.coefs$term == 'slope' & dl1.coefs$model == 'first_v_end'],
              lty = 1,
              alpha = 0.5,
              col = 'red',
              lwd = 1) +
  coord_cartesian(xlim = c(-25,15), ylim = c(-25,10)) ->
  first_v_end_plot

ggpubr::ggarrange(plotlist = list(first_v_peak_plot, first_v_end_plot),
                  nrow = 2,
                  ncol = 1,
                  labels = c('a','b')) -> iler_deming_plot

# save xy plot of onset versus peak and last flowering with stderr based on Weibull estimates
ggsave(iler_deming_plot,
       filename = file.path('../output/iler_deming_plot.png'),
       device = 'png',
       width = 5,
       height = 8,
       units = 'in',
       dpi = 300)

##################################################################################
##                                                                              ##
##         ~~  CODE FOR PEARSE ET AL. (2017), FROM W. PEARSE  ~~                ##
##                                                                              ##
##################################################################################

#Modelling Rocky Mountain data

#! Comments added by Will Pearse 2018-07-09. This is the analysis
#! script, with some additional comments (to explain) and some file
#! load path/name changes.
#! Comments added afterwards have "#!" at the start of them, not "#"

#! This became the released code in the manuscript. Original path/filename: "../code/headers.R"
source("./2018-07-09_headers.R")

#Load data
#! This came from a file called "phenology data_1973_2012 from Amy
#! Iler.xlsx".
#data <- na.omit(read.csv("amy_iler_sheet2.csv"))
data <- na.omit(read.csv("../data/rmbl_pnas_data.csv"))
names(data) <- tolower(names(data))
names(data)[4] <- "n.flowers"
counts <- apply(with(data, tapply(doy, list(species,year), function(x) length(x) > 0)), 1, function(x) sum(x, na.rm=TRUE) >= 19)
data <- data[data$species %in% names(counts)[counts==TRUE],]
data <- data[order(data$species, data$plot, data$year, data$doy),]
data <- data[data$n.flowers >= 1,]
data$n.flowers <- round(data$n.flowers)
data <- do.call(rbind, apply(data, 1, function(x) sapply(x[-4], function(y) rep(y, x[4]))))

#! This is where Pearse et al. replicate each entry
#! calculate by how much:
fold_replication <- dim(data)[1] / dim(unique(data))[1] #! factor of 23-fold longer data.frame

data <- 
  data %>%
  as.data.frame() %>%
  dplyr::arrange(species, plot, year, doy)

data$doy <- as.numeric(as.character(data$doy))
data$year <- as.numeric(as.character(data$year))
data$species <- as.character(data$species)

#! dedup data and run through same routine:
data_dedup <- unique(data)

#Calculate everything
theta <- as.numeric(with(data, tapply(doy, list(species, plot, year), function(x) if(length(x) > 10) theta.hat(x,k=30) else NA)))
theta.end <- as.numeric(with(data, tapply(doy, list(species, plot, year), function(x) if(length(x) > 10) theta.hat(x,k=30,TRUE) else NA)))
min <- as.numeric(with(data, tapply(doy, list(species, plot, year), function(x) if(length(x) > 10) min(x) else NA)))
max <- as.numeric(with(data, tapply(doy, list(species, plot, year), function(x) if(length(x) > 10) max(x) else NA)))
mean <- as.numeric(with(data, tapply(doy, list(species, plot, year), function(x) if(length(x) > 10) mean(x,na.rm=TRUE) else NA)))
median <- as.numeric(with(data, tapply(doy, list(species, plot, year), function(x) if(length(x) > 10) median(x,na.rm=TRUE) else NA)))
abundance <- as.numeric(with(data, tapply(doy, list(species, plot, year), function(x) if(length(x) > 10) max(table(x)) else NA)))
species <- as.character(with(data, tapply(species, list(species, plot, year), unique)))
year <- as.numeric(with(data, tapply(year, list(species, plot, year), unique)))
results <- na.omit(data.frame(theta, theta.end, min, max, mean, species, year, median, abundance))
results$orig.year <- results$year
results$year <- as.numeric(scale(results$year))
results_dup <- results[results$species %in% names(Filter(function(x) x > 10, table(results$species))),]
rm(theta,min,max,mean,median,year,species,theta.end)

#! RE-GENERATE RESULTS WITH DATA_DEDUP
theta <- as.numeric(with(data_dedup, tapply(doy, list(species, plot, year), function(x) if(length(x) > 10) theta.hat(x,k=30) else NA)))
theta.end <- as.numeric(with(data_dedup, tapply(doy, list(species, plot, year), function(x) if(length(x) > 10) theta.hat(x,k=30,TRUE) else NA)))
min <- as.numeric(with(data_dedup, tapply(doy, list(species, plot, year), function(x) if(length(x) > 10) min(x) else NA)))
max <- as.numeric(with(data_dedup, tapply(doy, list(species, plot, year), function(x) if(length(x) > 10) max(x) else NA)))
mean <- as.numeric(with(data_dedup, tapply(doy, list(species, plot, year), function(x) if(length(x) > 10) mean(x,na.rm=TRUE) else NA)))
median <- as.numeric(with(data_dedup, tapply(doy, list(species, plot, year), function(x) if(length(x) > 10) median(x,na.rm=TRUE) else NA)))
abundance <- as.numeric(with(data_dedup, tapply(doy, list(species, plot, year), function(x) if(length(x) > 10) max(table(x)) else NA)))
species <- as.character(with(data_dedup, tapply(species, list(species, plot, year), unique)))
year <- as.numeric(with(data_dedup, tapply(year, list(species, plot, year), unique)))
results <- na.omit(data.frame(theta, theta.end, min, max, mean, species, year, median, abundance))
results$orig.year <- results$year
results$year <- as.numeric(scale(results$year))
results_dedup <- results[results$species %in% names(Filter(function(x) x > 10, table(results$species))),]
rm(theta,min,max,mean,median,year,species,theta.end)

#! compare results_dup with results_dedup
#! another effect is that when year is re-scaled it is done incorrectly because of the duplicated rows
#! Check how many estimates per species-year were generated by the routine above:
results_dup %>%
  dplyr::group_by(species, orig.year) %>%
  dplyr::summarise(n_theta = n()) %>%
  dplyr::arrange(species, orig.year) ->
  n_thetas_dup

results_dedup %>%
  dplyr::group_by(species, orig.year) %>%
  dplyr::summarise(n_theta = n()) %>%
  dplyr::arrange(species, orig.year) ->
  n_thetas_dedup

#! they're totally different numbers of estimates.
#! I need to check whether, if done correctly, it makes a huge difference.
#! combine data.frames to compare
results_dup %>%
  dplyr::group_by(species, orig.year) %>%
  dplyr::summarise(n_theta = n(),
                   mean_theta = mean(theta, na.rm = T),
                   sd_theta   = sd(theta, na.rm = T)) ->
  results_dup2

results_dedup %>%
  dplyr::group_by(species, orig.year) %>%
  dplyr::summarise(n_theta = n(),
                   mean_theta = mean(theta, na.rm = T),
                   sd_theta   = sd(theta, na.rm = T)) ->
  results_dedup2

res_comp <- 
  dplyr::full_join(results_dup2,results_dedup2, by = c('species','orig.year'))

# plot:
res_comp %>%
  ggplot(aes(x = mean_theta.x, y = mean_theta.y)) +
  geom_point(alpha = 0.33) +
  geom_abline(slope = 1, intercept = 0) +
  coord_cartesian(xlim = c(-100,300), ylim = c(-100,300)) +
  theme_bw()

# plot difference
res_comp %>%
  dplyr::mutate(theta_diff = mean_theta.y - mean_theta.x) ->
  res_comp

# re-plot
res_comp %>%
  ggplot(aes(x = theta_diff)) +
  geom_histogram(bins = 100, alpha = 0.66) +
  theme_bw() +
  scale_x_continuous(limits = c(-50, 50)) +
  geom_vline(xintercept = 0)

# systematic bias towards later first flowering date with the pseudo-replication


#Modelling
mean.model  <- lm(mean  ~ species * year + log(abundance), data = results)
theta.model <- lm(theta ~ species * year + log(abundance), data = results)
theta.end.model <- lm(theta.end ~ species * year + log(abundance), data=results)

#...logging doesn't do much to model fit, but it makes abundance
#   matter more (--> conservative?) and it looks like it should be logged

# OLD coefficient code (Amy added annotation here)
coefs_P1 <- with(results, as.data.frame(cbind(
  summary(lm(theta ~ species:orig.year + log(abundance)))$coefficients[-1:-2,1:2],
  summary(lm(mean ~ species:orig.year + log(abundance)))$coefficients[-1:-2,1:2],
  summary(lm(theta.end ~ species:orig.year + log(abundance)))$coefficients[-1:-2,1:2]
)))

# CORRECTED coefficient code (Amy added annotation here)
coefs_P2 <- with(results, as.data.frame(cbind(
  summary(lm(theta ~ species+species:orig.year + log(abundance)-1))$coefficients[-1:-66,1:2],
  summary(lm(mean ~ species+species:orig.year + log(abundance)-1))$coefficients[-1:-66,1:2],
  summary(lm(theta.end ~ species+species:orig.year + log(abundance)-1))$coefficients[-1:-66,1:2]
)))

coefs_P1$species <- unique(results$species)
coefs_P2$species <- unique(results$species)

names(coefs_P1)[1:6] <- c("start", "start.se", "bulk", "bulk.se", "end", "end.se")
names(coefs_P2)[1:6] <- c("start", "start.se", "bulk", "bulk.se", "end", "end.se")

#...undo the scaling of year...
#coefs$start <- coefs$start / sd(orig.year)

# Deming Regressions (Amy added this part)

# first vs. peak
d.fvp2 <- deming(coefs_P2$bulk ~ coefs_P2$start, xstd = coefs_P2$start.se, ystd = coefs_P2$bulk.se)

# first vs. end
d.fve2 <- deming(coefs_P2$end ~ coefs_P2$start, xstd = coefs_P2$start.se, ystd = coefs_P2$end.se)

# peak vs. first
d.pvf2 <- deming(coefs_P2$start ~ coefs_P2$bulk, xstd = coefs_P2$bulk.se, ystd = coefs_P2$start.se)

# peak vs. end
d.pve2 <- deming(coefs_P2$end ~ coefs_P2$bulk, xstd = coefs_P2$bulk.se, ystd = coefs_P2$end.se)

list(d.fvp2,
     d.fve2,
     d.pvf2,
     d.pve2) -> dl2

names(dl2) <- c('first_v_peak',
                'first_v_end',
                'peak_v_first',
                'peak_v_end')

dl2.coefs <- lapply(dl2, function(x) data.frame(x$coefficients,x$ci))

for (i in seq_along(dl2.coefs)) {
  dl2.coefs[[i]]$model <- names(dl2)[i]
}

dl2.coefs <- 
  dl2.coefs %>%
  do.call(rbind, .)

names(dl2.coefs)[1] <- c('estimate')
dl2.coefs$term <- 'slope'
dl2.coefs$term[grep('Intercept', row.names(dl2.coefs))] <- 'intercept'

write.table(dl2.coefs,
            file = '../output/deming_reg_coefs_Pearse.csv',
            quote = F,
            sep = ',',
            row.names = F,
            col.names = T)

# plotting the coefficient estimates from Deming regressions against one another:
all_coefs <- 
  dplyr::left_join(dl1.coefs,dl2.coefs, by = c('model','term'))

# plot coefficient estimates from both approaches
all_coefs %>%
  dplyr::filter(term == 'slope') %>%
  ggplot() +
  geom_abline(slope = 1, intercept = 0, col = 'gray40', lty = 1, lwd = 0.25) +
  geom_hline(yintercept = 1, lty = 1, col = 'gray40', lwd = 0.25) +
  geom_vline(xintercept = 1, lty = 1, col = 'gray40', lwd = 0.25) +
  geom_errorbarh(aes(y = estimate.y,
                    xmin = lower.0.95.x,
                    xmax = upper.0.95.x, col = model), alpha = 0.5) +
  geom_errorbar(aes(x = estimate.x,
                    ymin = lower.0.95.y,
                    ymax = upper.0.95.y, col = model), width = 0, alpha = 0.5) +
  geom_point(aes(x = estimate.x, y = estimate.y, col = model), alpha = 0.5) +
  coord_cartesian(xlim = c(0.4,2.1),
                  ylim = c(0.4,2.1)) +
  xlab('Iler estimates') +
  ylab('Pearse estimates') +
  theme_bw() -> p1

ggsave(p1,
       filename = 'coef_comparison.png',
       width = 4,
       height = 3,
       units = 'in',
       device = 'png',
       dpi = 300
)

################################################################################
##                                                                            ##
##  ~~  COMPARING THETA ESTIMATES BTWN Pearse et al. and Iler re-analysis  ~~ ##
##                                                                            ##
################################################################################

all_comb <- dplyr::left_join(
  dplyr::select(all.data,species,orig.year,theta = onset, theta.end = end, onset.obs = min.doy, end.obs = max.doy),
  dplyr::select(results,species,orig.year,theta,theta.end, onset.obs = min, end.obs = max),
  by = c('species','orig.year')) %>%
  dplyr::filter(theta.y > 0)

# plot
all_comb %>%
  ggplot() +
  geom_point(aes(x = theta.x, y = theta.y), alpha = 0.1, col = 'black') +
  theme_bw() +
  coord_cartesian(xlim = c(50,260),
                  ylim = c(50,260)) +
  xlab(TeX("$\\theta_{start}$  Iler")) +
  ylab(TeX("$\\theta_{start}$  Pearse")) ->
  theta_start_plot

all_comb %>%
  ggplot() +
  geom_point(aes(x = theta.end.x, y = theta.end.y), alpha = 0.1, col = 'black') +
  theme_bw() +
  coord_cartesian(xlim = c(100,350),
                  ylim = c(100,350)) +
  xlab(TeX("$\\theta_{end}$  Iler")) +
  ylab(TeX("$\\theta_{end}$  Pearse")) ->
  theta_end_plot

ggpubr::ggarrange(plotlist = list(theta_start_plot,theta_end_plot), labels = c('a','b'),
                  ncol = 2) %>%
  ggsave(filename = '../output/theta_comps.png',
         width = 6,
         height = 3,
         units = 'in',
         device = 'png',
         dpi = 300)

# # now plot difference between estimated and observed onset and end
# all_comb %>%
#   ggplot() +
#   geom_histogram(aes(x = onset.obs.x - theta.x), bins = 100) +
#   theme_bw() +
#   xlab(TeX("$\\theta_{start} - t_{obs}$  Iler"))
# 
# all_comb %>%
#   ggplot() +
#   geom_histogram(aes(x = onset.obs.y - theta.y), bins = 100) +
#   theme_bw() +
#   xlab(TeX("$\\theta_{start} - t_{obs}$  Pearse"))

