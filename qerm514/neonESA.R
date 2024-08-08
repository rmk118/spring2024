# Ruby Krasnow
# NEON/EFI workshop - ESA Annual Meeting 
# 2024-08-08

# https://www.neonscience.org/resources/learning-hub/tutorials/neon-beetle-forecasting

# List of packages required:
packages <- c("tidyverse", "lubridate", "tsibble", "fable", "fabletools", "remotes", "neon4cast", "score4cast")

# Load packages into session
lapply(packages, require, character.only = TRUE)
rm(packages)

# To download the NEON site information table:
neon_site_info <- read_csv("https://www.neonscience.org/sites/default/files/NEON_Field_Site_Metadata_20231026.csv")

# choose site
my_site = "OSBS"

# date where we will start making predictions
forecast_startdate <- "2022-01-01" #fit up through 2021, forecast 2022 data

# date where we will stop making predictions
forecast_enddate <- "2025-01-01"

# beetle targets are here
url <- "https://sdsc.osn.xsede.org/bio230014-bucket01/challenges/targets/project_id=neon4cast/duration=P1W/beetles-targets.csv.gz"

# read in the table
targets <- read_csv(url) %>%
  mutate(datetime = as_date(datetime)) %>%  # set proper formatting
  dplyr::filter(site_id == my_site,  # filter to desired site
                datetime < "2022-12-31") # excluding provisional data

# initial viz
targets %>% 
  as_tsibble(index = datetime, key = variable) %>%
  autoplot() +
  facet_grid(variable ~ ., scales = "free_y") + 
  theme_bw() +
  theme(legend.position = "none")

targets_train <- targets %>%
  filter(datetime < forecast_startdate) %>%
  pivot_wider(names_from = variable, values_from = observation) %>%
  as_tsibble(index = datetime)

# specify and fit models

# Using a log(x + 1) transform on the abundance data
mod_fits <- targets_train %>% 
  tsibble::fill_gaps() %>% # gap fill the data for the random walk model
  fabletools::model(
    mod_mean = fable::MEAN(log1p(abundance)), # generate a forecast from the historical mean plus the standard deviation of the historical data
    mod_naive = fable::NAIVE(log1p(abundance))) # random walk model, requires gap filling. Generates a forecast from the current observation plus random process noise

# make a forecast
fc_null <- mod_fits %>%
  fabletools::forecast(h = "3 years") 

# visualize the forecast
fc_null %>% 
  autoplot(targets_train) +
  facet_grid(.model ~ ., scales = "free_y") +
  theme_bw()

# Get climate data
path_to_clim_data <- "https://data.cyverse.org/dav-anon/iplant/projects/NEON/ESA2024/forecasting_beetles_workshop/modeled_climate_2012-2050_OSBS_CMCC_CM2_VHR4.csv"

clim_long <- read_csv(path_to_clim_data)  %>%
  filter(datetime <= forecast_enddate)

# make a tsibble object
clim_long_ts <- clim_long %>%
  as_tsibble(index = datetime, 
             key = c(variable, model_id))

# make wide
clim_wide <- clim_long %>%
  select(-unit) %>%
  pivot_wider(names_from = variable, values_from = prediction)

# visualize climate data
clim_long_ts %>% 
  ggplot(aes(datetime, prediction)) + 
  geom_line() +
  facet_grid(variable ~ ., scales = "free_y") +
  geom_vline(xintercept = lubridate::as_date(forecast_startdate),
             lty = 2) + 
  theme_bw() +
  theme(legend.position = "none")

# subset into past and future datasets, based on forecast_startdate
clim_past <- clim_wide %>%
  filter(datetime < forecast_startdate,
         datetime > "2012-01-01")

clim_future <- clim_wide %>%
  filter(datetime >= forecast_startdate,
         datetime <= forecast_enddate)

# combine target and climate data into a training dataset
targets_clim_train <- targets_train %>%
  left_join(clim_past)


# specify and fit model
mod_fit_candidates <- targets_clim_train %>%
  fabletools::model(
    mod_temp = fable::TSLM(log1p(abundance) ~ temperature_2m_mean),
    mod_precip = fable::TSLM(log1p(abundance) ~ precipitation_sum),
    mod_temp_precip = fable::TSLM(log1p(abundance) ~ temperature_2m_mean + precipitation_sum))

# look at fit stats and identify the best model using AIC, where the lowest AIC is the best model
glance(mod_fit_candidates) %>%
  select(`.model`, AIC)

# visualize model fit
# augment reformats model output into a tsibble for easier plotting
fabletools::augment(mod_fit_candidates) %>%
  ggplot(aes(x = datetime)) +
  geom_line(aes(y = abundance, lty = "Obs"), color = "dark gray") +
  geom_line(aes(y = .fitted, color = .model, lty = "Model")) +
  facet_grid(.model ~ .) +
  theme_bw()

# focus on temperature model for now, it has the lowest AIC
mod_best_lm <- mod_fit_candidates %>% select(mod_temp)
report(mod_best_lm)

# make a forecast
# filter "future" climate data to target climate model
fc_best_lm <- mod_best_lm %>%
  fabletools::forecast(
    new_data = 
      clim_future %>%
      as_tsibble(index = datetime)) 

# visualize the forecast
fc_best_lm %>% 
  autoplot(targets_train) +
  facet_grid(.model ~ .) + 
  theme_bw()

# update dataframe of model output for submission
fc_best_lm_efi <- fc_best_lm %>% 
  mutate(site_id = my_site) %>% #efi needs a NEON site ID
  neon4cast::efi_format() %>%
  mutate(
    project_id = "neon4cast",
    model_id = "bet_abund_example_tslm_temp",
    reference_datetime = forecast_startdate,
    duration = "P1W")

head(fc_best_lm_efi)

# scoring the forecast
# filter to 2022 because that is the latest release year

# 2023 is provisional and most sites do not yet have data reported

targets_2022 <- targets %>% 
  dplyr::filter(
    datetime >= "2022-01-01", 
    datetime < "2023-01-01",
    variable == "abundance",
    observation > 0)

# list of target site dates for filtering mod predictions
target_site_dates_2022 <- targets_2022 %>%
  select(site_id, datetime) %>% distinct()

# filter model forecast data to dates where we have observations
mod_results_to_score_lm <- fc_best_lm_efi %>%
  left_join(target_site_dates_2022,.) %>%
  dplyr::filter(!is.na(parameter))

# score the forecasts
mod_scores <- score(
  forecast = mod_results_to_score_lm,
  target = targets_2022) 

head(mod_scores)

# get scores for the mean and naive models

# the fc_null object has scores from both models
# note: we need to add site_id back in for efi_format() to work
fc_null_efi <- fc_null %>% 
  mutate(site_id = my_site) %>% #efi needs a NEON site ID
  neon4cast::efi_format() 

# filter to dates where we have target data from 2022
mod_results_to_score_null <- fc_null_efi %>%
  left_join(target_site_dates_2022,.) %>%
  dplyr::filter(!is.na(parameter))

# score the forecasts for those dates
mod_null_scores <- score(
  forecast = mod_results_to_score_null,
  target = targets_2022) 

# stack the scores for our best_lm and the null models
# forcing reference_datetime to be the same type in both tables
# so they will stack

all_mod_scores <- bind_rows(
  mod_null_scores %>% mutate(
    reference_datetime = as.character(reference_datetime)), 
  mod_scores %>% mutate(
    reference_datetime = as.character(reference_datetime)))

mod_results_to_score_lm <- mod_results_to_score_lm %>%
  select(site_id,datetime,parameter,model_id,family,variable,prediction)

rbind(mod_results_to_score_null, mod_results_to_score_lm) %>% 
  ggplot(., aes(datetime, prediction, color = model_id, group=interaction(parameter, model_id))) +
  geom_line(lwd = 1)+
  geom_point(aes(datetime, observation), color = "black", size = 6, inherit.aes = F, data = targets_2022)+
  ylab("Abundance")+
  theme_classic()

all_mod_scores %>%
  ggplot(aes(datetime, crps, color = model_id)) +
  geom_line() +
  theme_bw() +
  ggtitle("Forecast scores over time at OSBS for 2022")
