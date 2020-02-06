# 2020-JAN-30
# PTH

library(dplyr)
library(ggplot2)

# load data; put in pwd
data <- na.omit(read.csv("../data/rmbl_pnas_data.csv"))

names(data)[4] <- "flowers"
names(data) <- tolower(names(data))

## ============= ##
## two concerns: ##
## ============= ##

## 1. In original data, there appears to be duplicated rows:
(diff_in_rows <- nrow(data) - nrow(unique(data)))
# 30 entries are duplicated exactly

## 2. Many species have multiple flower counts per year-plot-doy combination.
# This is unexpected if each species gets a total count per plot per day per year.

data %>%
  dplyr::group_by(species, year, plot, doy, habitat) %>%
  dplyr::summarise(n_rows = n()) %>%
  dplyr::filter(n_rows > 1) %>%
  dplyr::arrange(desc(n_rows)) -> tmp1

tmp1 %>%
  dplyr::group_by(n_rows) %>%
  dplyr::summarise(count_of_dups = n()) ->
  count_of_dups

count_of_dups
# 145 entries with 2 rows
# 15 entries with 3 rows
# 1 entry with 12 rows

# let's take a closer look at the top offender:
data %>%
  dplyr::filter(species == tmp1$species[1],
                year    == tmp1$year[1],
                plot    == tmp1$plot[1],
                doy     == tmp1$doy[1])

# there are 12 rows of data with the same plot-year-doy for this species and it is very confusing.
# what's going on?