
#title: "Weibull"
#output: html_document
#author: Parris T Humphrey
#date: 2020-FEB-05

# Amy edits to Parris's R Markdown file (easier for me to run this in a 'regular' script)
# August 2020

#####################
##                 ##
##  PRELIMINARIES  ##
##                 ##
#####################

setwd("~/Dropbox/AMY/Research/Pearse Response")

library(reshape2)
library(deming)
library(car)
library(devtools)
library(broom)
library(plyr)
library(lme4)
#library(phest) not available for my version of R; should be fine since Parris loaded the functions from github in #a separate script
library(lmerTest)
#library(dplyr) # not sure why these three are loaded, because they are part of the tidyverse; knitr too?
#library(stringr)
#library(ggplot2)
library(knitr)
library(kableExtra)
library(ggpubr)
library(tidyverse)

source("iler_functions.R")
source("weibull.R")


####################
##                ##
##      DATA      ##
##                ##
####################

data <- na.omit(read.csv("amy_iler_sheet2.csv"))
# data <- na.omit(read.csv("amy_iler_sheet2_duplicates corrected.csv")) # some duplicated rows from the original conversion from spreadsheets to R dataframe
# results are the same if we use this file

# rename some headers; make all lowercase
names(data)[4] <- "flowers"
names(data) <- tolower(names(data))

# subset data to include only years after 1973
data <- subset(data, year > '1973')

## Process data in two ways:
# first WITH flower-fold pseudoreplication of each species-year-plot-doy in data
# second WITHOUT this step
data_dup   <- preprocess_data(data, min_years = 19, pseudoreplicate = TRUE)
data_dedup <- preprocess_data(data, min_years = 19, pseudoreplicate = FALSE)

## Generate Weibull estimates
# run estimates with and without plot-level aggregation using replicated and non-replicated input data
# use Pearse et al. (2017) default k_min

# replicated input data
thetas_plot_dup     <- calc_weibull_pearse(data_dup, by_plot = TRUE,  k_min = 10) # plot is gone as a variable, but there is still an estimate for each plot
thetas_noplot_dup   <- calc_weibull_pearse(data_dup, by_plot = FALSE, k_min = 10)

# non-replicated input data
thetas_plot_dedup   <- calc_weibull_pearse(data_dedup, by_plot = TRUE,  k_min = 10) 
thetas_noplot_dedup <- calc_weibull_pearse(data_dedup, by_plot = FALSE, k_min = 10)


#####################
##                 ##
##     ANALYSIS    ##
##                 ##
#####################

##### (1) What are the effects of replication (Parris's 3a; Figure S2)
### panel (a)
# compare actual estimates by taking means
thetas_plot_dup %>%
  dplyr::group_by(species, year) %>%
  dplyr::summarise(n_theta = n(),
                   mean_theta = mean(theta_onset, na.rm = T),
                   sd_theta   = sd(theta_onset, na.rm = T)) ->
  thetas_plot_dup2

thetas_plot_dedup %>%
  dplyr::group_by(species, year) %>%
  dplyr::summarise(n_theta = n(),
                   mean_theta = mean(theta_onset, na.rm = T),
                   sd_theta   = sd(theta_onset, na.rm = T)) ->
  thetas_plot_dedup2

res_comp_a <-
  dplyr::full_join(thetas_plot_dup2, thetas_plot_dedup2, by = c('species','year'))

# plot:
res_comp_a %>%
  ggplot(aes(x = mean_theta.x, y = mean_theta.y)) +
  geom_point(alpha = 0.33) +
  geom_abline(slope = 1, intercept = 0) +
  coord_cartesian(xlim = c(-100,300), ylim = c(-100,300)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab('estimated onset (replicated dates)') + 
  ylab('estimated onset (non-replicated dates)') ->
  xyplot_a

# plot difference
res_comp_a %>%
  dplyr::mutate(theta_diff = mean_theta.y - mean_theta.x) ->
  res_comp_a

# re-plot
res_comp_a %>%
  ggplot(aes(x = theta_diff)) +
  geom_histogram(bins = 100, alpha = 0.66) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(limits = c(-30, 30), breaks = seq(-30,30,10)) +
  scale_y_continuous(limits = c(0, 300)) +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = mean(res_comp_a$theta_diff, na.rm = T), col = 'darkorange2', lty = 2) +
  xlab('non-replicated - replicated dates (estimated)') ->
  diff_plot_a

plot_a <- ggpubr::ggarrange(plotlist = list(xyplot_a, diff_plot_a), ncol = 2, align = 'hv')

### panel (b)
res_comp_b <-
  dplyr::full_join(thetas_noplot_dup, thetas_noplot_dedup, by = c('species','year'))

res_comp_b %>%
  ggplot(aes(x = theta_onset.x, y = theta_onset.y)) +
  geom_point(alpha = 0.33) +
  geom_abline(slope = 1, intercept = 0) +
  coord_cartesian(xlim = c(-100,300), ylim = c(-100,300)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab('estimated onset (replicated dates)') + 
  ylab('estimated onset (non-replicated dates)') ->
  xyplot_b

# plot difference
res_comp_b %>%
  dplyr::mutate(theta_diff = theta_onset.y - theta_onset.x) ->
  res_comp_b

# re-plot
res_comp_b %>%
  ggplot(aes(x = theta_diff)) +
  geom_histogram(bins = 100, alpha = 0.66) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(limits = c(-30, 30), breaks = seq(-30,30,10)) +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = mean(res_comp_b$theta_diff, na.rm = T), col = 'darkorange2', lty = 2) +
  xlab('non-replicated - replicated dates (estimated)') ->
  diff_plot_b

plot_b <- ggpubr::ggarrange(plotlist = list(xyplot_b, diff_plot_b), ncol = 2, align = 'hv')

ggpubr::ggarrange(plotlist = list(plot_a, plot_b), nrow = 2, ncol = 1, align = 'hv', labels = c('a','b'))


##### (2) What is the effect of plot-level versus year-level estimates? (Parris's 3c)

### panel a:
# bring in min
thetas_plot_dup %>%
  dplyr::group_by(species, year) %>%
  dplyr::summarise(min_obsv = min(min)) ->
  the_mins

res_comp_a %>%
  dplyr::left_join(the_mins, by = c('species','year')) %>%
  dplyr::mutate(diff_from_min.x = mean_theta.x - min_obsv, # theta.x is from replicated dates above
                diff_from_min.y = mean_theta.y - min_obsv) -> # theta.y is from non-replicated dates above
  res_comp2
view(res_comp2)
#write.csv(res_comp2, file="res_comp2.csv")

# now plot:
res_comp2 %>%
  ggplot(aes(x = diff_from_min.x)) +
  geom_histogram(bins = 100, alpha = 0.66) +
  theme_bw() +
  scale_x_continuous(limits = c(-30, 30)) +
  geom_vline(xintercept = 0) +
  ggtitle("replicated dates") +
  xlab("estimated onset - observed onset") +
  theme(plot.title = element_text(size = 8, face = 'bold'),
        panel.grid = element_blank()) ->
  diff_from_min_plot.x

res_comp2 %>%
  ggplot(aes(x = diff_from_min.y)) +
  geom_histogram(bins = 100, alpha = 0.66) +
  theme_bw() +
  scale_x_continuous(limits = c(-30, 30)) +
  geom_vline(xintercept = 0) +
  ggtitle("non-replicated dates") + 
  xlab("estimated onset - observed onset") +
  theme(plot.title = element_text(size = 8, face = 'bold'),
        panel.grid = element_blank()) ->
  diff_from_min_plot.y

diff_plot_a <- ggpubr::ggarrange(plotlist = list(diff_from_min_plot.x, diff_from_min_plot.y), ncol = 2, align = 'hv', labels = c('a','b'))

### panel b
thetas_noplot_dup %>%
  dplyr::group_by(species, year) %>%
  dplyr::summarise(min_obsv = min(min)) ->
  the_mins

res_comp_b %>%
  dplyr::left_join(the_mins, by = c('species','year')) %>%
  dplyr::mutate(diff_from_min.x = theta_onset.x - min_obsv,
                diff_from_min.y = theta_onset.y - min_obsv) ->
  res_comp2 # should this be named something different to make sure it does not override the previous res_comp2?
view(res_comp2)

# now plot:
res_comp2 %>%
  ggplot(aes(x = diff_from_min.x)) +
  geom_histogram(bins = 100, alpha = 0.66) +
  theme_bw() +
  scale_x_continuous(limits = c(-30, 30)) +
  geom_vline(xintercept = 0) +
  ggtitle("replicated dates") + 
  xlab("estimated onset - observed onset") +
  theme(plot.title = element_text(size = 8, face = 'bold'),
        panel.grid = element_blank()) ->
  diff_from_min_plot.x # these should be re-named too

res_comp2 %>%
  ggplot(aes(x = diff_from_min.y)) +
  geom_histogram(bins = 100, alpha = 0.66) +
  theme_bw() +
  scale_x_continuous(limits = c(-30, 30)) +
  geom_vline(xintercept = 0) +
  ggtitle("non-replicated dates") +
  xlab("estimated onset - observed onset") +
  theme(plot.title = element_text(size = 8, face = 'bold'),
        panel.grid = element_blank()) ->
  diff_from_min_plot.y

diff_plot_b <- ggpubr::ggarrange(plotlist = list(diff_from_min_plot.x, diff_from_min_plot.y), ncol = 2, align = 'hv')

ggpubr::ggarrange(plotlist = list(diff_plot_a, diff_plot_b), nrow = 2, ncol = 1, align = 'hv', labels = c('a','b'))


##### (3) What is the effect of scaling the year variable improperly (Parris's 3d)

## re-scaling year

# first, done the wrong way:
thetas_plot_dup_A   <- rescale_year(thetas_plot_dup,   uniques = FALSE, fill_missing = FALSE)
thetas_noplot_dup_A <- rescale_year(thetas_noplot_dup, uniques = FALSE, fill_missing = FALSE)

thetas_plot_dedup_A   <- rescale_year(thetas_plot_dedup,   uniques = FALSE, fill_missing = FALSE)
thetas_noplot_dedup_A <- rescale_year(thetas_noplot_dedup, uniques = FALSE, fill_missing = FALSE)

# second, done the correct way:
thetas_plot_dup_B   <- rescale_year(thetas_plot_dup,   uniques = TRUE, fill_missing = TRUE)
thetas_noplot_dup_B <- rescale_year(thetas_noplot_dup, uniques = TRUE, fill_missing = TRUE)

thetas_plot_dedup_B   <- rescale_year(thetas_plot_dedup,   uniques = TRUE, fill_missing = TRUE)
thetas_noplot_dedup_B <- rescale_year(thetas_noplot_dedup, uniques = TRUE, fill_missing = TRUE)

pz1 <- plot_rescaled_year(thetas_plot_dup_A, 'dup, by plot')
pz2 <- plot_rescaled_year(thetas_noplot_dup_A, 'dup, by year')
pz3 <- plot_rescaled_year(thetas_plot_dedup_A, 'dedup, by plot')
pz4 <- plot_rescaled_year(thetas_noplot_dedup_A, 'dedup, by year')

ggpubr::ggarrange(plotlist = list(pz1, pz2, pz3, pz4), ncol = 2, nrow = 2,
                  labels = c('a','b','c','d'))

# The intervals between years is not preserved in any of the datasets that rescale year without first 
# taking the unique values and filling in any missing entries. This will produce errors in the relationship 
# between onset and year in the regression moels that follow and will confound our inference of the relationship 
# between onset, peak, or end-of-flowering versus time.


#####################
##                 ##
##   RE-ANALYSIS   ##
##                 ##
#####################

thetas_pearse <- thetas_plot_dup_A
thetas_iler   <- thetas_noplot_dedup_B

# AMY wrote the code below to get rid of species with fewer than 19 years of data
thetas_iler <- thetas_iler %>% 
  dplyr::group_by(species) %>%
  dplyr::mutate(n_years = n()) %>% 
  dplyr::filter(n_years >= 19)

thetas_pearse <- thetas_pearse %>% 
  dplyr::group_by(species) %>%
  dplyr::mutate(n_years = n()) %>% 
  dplyr::filter(n_years >= 19)

unique(thetas_iler$species) # 52 species
unique(thetas_pearse$species) # 60 species

# How close are estimated values to observed?
min(thetas_iler$theta_onset)
max(thetas_iler$theta_onset)
summary_thetas_iler <- thetas_iler %>% 
  dplyr::mutate(onset_diff = min-theta_onset, end_diff = theta_end-max) %>% 
  dplyr::group_by(species) %>% 
  dplyr::summarise(mean_onset_diff = mean(onset_diff), mean_end_diff = mean(end_diff))
mean(summary_thetas_iler$mean_onset_diff)
sd(summary_thetas_iler$mean_onset_diff)/sqrt(52)
mean(summary_thetas_iler$mean_end_diff)
sd(summary_thetas_iler$mean_end_diff)/sqrt(52)

min(thetas_pearse$theta_onset)

# add em50 dataframe (em50 = day on which 50% of flowers were counted)
em50.wide <- read.csv("EM50_RMBL.csv", header = T)
em50.long <- melt(em50.wide, id.vars = 'year')
names(em50.long)[2] <- "species"
names(em50.long)[3] <- "em50"

em50.long$species <- as.character(em50.long$species)
em50.long$year    <- as.numeric(em50.long$year)
em50.long$species <- gsub("\\.", " ", em50.long$species)

thetas_iler %>%
  dplyr::left_join(em50.long, by = c('species','year')) ->
  thetas_iler

thetas_pearse %>%
  dplyr::left_join(em50.long, by = c('species','year')) ->
  thetas_pearse


# run linear regressions using onset, peak, and end estimates
lms_iler   <- run_lms(thetas_iler)
lms_pearse <- run_lms_p(thetas_pearse)

# run Deming regression and capture coefficient output
deming_iler   <- run_deming(lms_iler)
deming_pearse <- run_deming_p(lms_pearse) # Are these the numbers in the correction? NO.


# write output
write.table(deming_iler,
            file = '../deming_reg_coefs_Iler.csv',
            quote = F,
            sep = ',',
            row.names = F,
            col.names = T)

write.table(deming_pearse,
            file = '../deming_reg_coefs_Pearse.csv',
            quote = F,
            sep = ',',
            row.names = F,
            col.names = T)

##### Results tables

deming_iler %>%
  dplyr::arrange(term, model) %>%
  dplyr::bind_rows(
    dplyr::arrange(deming_pearse, term, model)
  ) %>%
  dplyr::filter(term == 'slope') %>%
  dplyr::mutate(estimate = round(estimate,3),
                lower.0.95 = round(lower.0.95,3),
                upper.0.95 = round(upper.0.95,3)) %>%
  kable() %>%
  kable_styling(bootstrap_options = c('condensed','striped'), full_width = FALSE) %>%
  pack_rows("Iler models", 1, 3) %>%
  pack_rows("Pearse models", 4, 6)

##### Results plots

# produce plots representing deming regression against one another
iler_deming_plot <- plot_coefs(lms_iler, deming_iler, 'Iler')
pearse_deming_plot <- plot_coefs(lms_pearse, deming_pearse, 'Pearse')

# save xy plot of onset versus peak and last flowering with stderr based on Weibull estimates
ggsave(iler_deming_plot,
       filename = file.path('../output/iler_deming_plot.png'),
       device = 'png',
       width = 5,
       height = 8,
       units = 'in',
       dpi = 300)

ggsave(pearse_deming_plot,
       filename = file.path('../output/pearse_deming_plot.png'),
       device = 'png',
       width = 5,
       height = 8,
       units = 'in',
       dpi = 300)

ggpubr::ggarrange(plotlist = list(iler_deming_plot, pearse_deming_plot),
                  labels = c('a','b'),
                  ncol = 2,
                  align = 'hv')