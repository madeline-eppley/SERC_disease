#### MarineGEO environmental data  for Severn, South, and West Rivers ####
## Madeline G. Eppley (MGE)
## data source: Maryland DNR/Eyes on the Bay - Fixed Station Monthly Monitoring
## salinity and temperature variables for ordinal logistic models

library(dplyr)
library(lubridate)
library(readr)
library(ggplot2)
library(tidyr)
library(patchwork)

setwd("~/Desktop/marineGEO_oysterdisease")

# this script will largely resemble the data processing for the seascape samples - I'm going to use my previous script structure 
download_date <- "01-21-2026"
source_description <- "Maryland DNR/Eyes on the Bay - Fixed Station Monthly Monitoring"
collection_type <- "monthly_fixed"

# station info df
stations <- data.frame(
  Station = c("WT7.1", "WT8.1", "WT8.3"),
  River = c("Severn", "South", "West"),
  lat = c(39.0068, 38.9493, 38.8501),
  lon = c(-76.5046, -76.5464, -76.5330),
  stringsAsFactors = FALSE)

## data download
# MGE downloaded raw data from EPA Chesapeake Bay Program DataHub on 1/21/26
# https://datahub.chesapeakebay.net/WaterQuality
# selected stations: WT7.1, WT8.1, WT8.3
# selected parameters: Salinity, Water Temperature
# selected date range: 01/01/2022 to 10/31/2025

raw_data <- read_csv("/Users/madelineeppley/Desktop/WaterQualityWaterQualityStation.csv", show_col_types = FALSE)
dim(raw_data) # 1666 x 30
colnames(raw_data)
head(raw_data[1:5])
head(raw_data[5:10])
head(raw_data[10:15])
head(raw_data[15:20]) # parameter col has salinity, MeasureValue has salinity measurement
head(raw_data[20:25]) # unit has PPT
print(table(raw_data$Parameter)) #833 observations of salinity and temp
print((table(raw_data$Station))) # a little bit uneven, WT8.3 has only 280 data points 

# clean data
cleaned_data <- raw_data %>%
  filter(Parameter %in% c("WTEMP", "SALINITY"),
         Station %in% c("WT7.1", "WT8.1", "WT8.3"))
dim(cleaned_data) #1666 x 30 - good, we retained all data

# date and time
cleaned_data <- cleaned_data %>%
  mutate(
    datetime = mdy_hms(paste(SampleDate, SampleTime)),
    year = year(datetime),
    month = month(datetime))

# add river names instead of station number to bind later with disease data
cleaned_data <- cleaned_data %>%
  left_join(stations %>% select(Station, River), by = "Station")
head(cleaned_data$River)

# split the temperature and salinity info into different data sets 
sal_data <- cleaned_data %>%
  filter(Parameter == "SALINITY") %>%
  rename(salinity = MeasureValue) %>%
  select(River, Station, datetime, year, month, salinity, Latitude, Longitude)

temp_data <- cleaned_data %>%
  filter(Parameter == "WTEMP") %>%
  rename(temp = MeasureValue) %>%
  select(River, Station, datetime, year, month, temp, Latitude, Longitude)

as.character(min(sal_data$datetime))
as.character(max(sal_data$datetime))
as.character(min(temp_data$datetime))
as.character(max(temp_data$datetime))

min(sal_data$salinity) # 2.1
max(sal_data$salinity) # 30.6
min(temp_data$temp) # 4.31
max(temp_data$temp) # 16.09

# get some summary stats for temp and salinity
sal_summary <- sal_data %>%
  group_by(River) %>%
  summarise(min = min(salinity), max = max(salinity), mean = mean(salinity), n = n())
print(sal_summary)

temp_summary <- temp_data %>%
  group_by(River) %>%
  summarise(min = min(temp), max = max(temp), mean = mean(temp), n = n())
print(temp_summary)

# add in oyster collection dates info
oyster_collection_dates <- data.frame(
  River = c("Severn", "South", "West", "Severn", "South", "West", 
            "Severn", "South", "West", "Severn", "South", "West"),
  Year = c(2022, 2022, 2022, 2023, 2023, 2023, 2024, 2024, 2024, 2025, 2025, 2025),
  collection_date = as.Date(c(
    "2022-10-18", "2022-10-18", "2022-10-18",
    "2023-08-22", "2023-08-22", "2023-08-22",
    "2024-09-27", "2024-09-27", "2024-09-27",
    "2025-08-25", "2025-08-25", "2025-08-25")))

# get quantiles for salinity in the 3 months before collection
sal_quantiles_precollection <- data.frame()
for (i in 1:nrow(oyster_collection_dates)) {
  river <- oyster_collection_dates$River[i]
  year <- oyster_collection_dates$Year[i]
  coll_date <- oyster_collection_dates$collection_date[i]
  start_date <- coll_date - months(3)
  
  river_data <- sal_data %>%
    filter(River == river, year %in% c(year - 1, year),
           datetime >= start_date, datetime <= coll_date)
  
  if (nrow(river_data) > 0) {
    sal_quantiles_precollection <- bind_rows(sal_quantiles_precollection,
                                             data.frame(
                                               River = river,
                                               Year = year,
                                               collection_date = coll_date,
                                               sal_q10 = quantile(river_data$salinity, 0.10),
                                               sal_q90 = quantile(river_data$salinity, 0.90),
                                               sal_n = nrow(river_data)))}}

temp_quantiles_precollection <- data.frame()
for (i in 1:nrow(oyster_collection_dates)) {
  river <- oyster_collection_dates$River[i]
  year <- oyster_collection_dates$Year[i]
  coll_date <- oyster_collection_dates$collection_date[i]
  start_date <- coll_date - months(3)
  
  river_data <- temp_data %>%
    filter(River == river, year %in% c(year - 1, year),
           datetime >= start_date, datetime <= coll_date)
  
  if (nrow(river_data) > 0) {
    temp_quantiles_precollection <- bind_rows(temp_quantiles_precollection,
                                              data.frame(
                                                River = river,
                                                Year = year,
                                                collection_date = coll_date,
                                                temp_q10 = quantile(river_data$temp, 0.10),
                                                temp_q90 = quantile(river_data$temp, 0.90),
                                                temp_n = nrow(river_data)))}}


env_quantiles_precollection <- sal_quantiles_precollection %>%
  left_join(temp_quantiles_precollection, by = c("River", "Year", "collection_date"))

sal_quantiles_summer <- data.frame()
for (yr in 2022:2024) {
  for (riv in unique(stations$River)) {
    summer_data <- sal_data %>%
      filter(River == riv, year == yr, month >= 6, month <= 8)
    
    if (nrow(summer_data) > 0) {
      sal_quantiles_summer <- bind_rows(sal_quantiles_summer,
                                        data.frame(
                                          River = riv,
                                          Year = yr,
                                          season = "summer",
                                          sal_q10 = quantile(summer_data$salinity, 0.10),
                                          sal_q90 = quantile(summer_data$salinity, 0.90),
                                          sal_mean = mean(summer_data$salinity),
                                          sal_n = nrow(summer_data)))}}}

temp_quantiles_summer <- data.frame()
for (yr in 2022:2024) {
  for (riv in unique(stations$River)) {
    summer_data <- temp_data %>%
      filter(River == riv, year == yr, month >= 6, month <= 8)
    
    if (nrow(summer_data) > 0) {
      temp_quantiles_summer <- bind_rows(temp_quantiles_summer,
                                         data.frame(
                                           River = riv,
                                           Year = yr,
                                           season = "summer",
                                           temp_q10 = quantile(summer_data$temp, 0.10),
                                           temp_q90 = quantile(summer_data$temp, 0.90),
                                           temp_mean = mean(summer_data$temp),
                                           temp_n = nrow(summer_data)))}}}

summer_quantiles <- sal_quantiles_summer %>%
  left_join(temp_quantiles_summer, by = c("River", "Year", "season"))

sal_monthly <- sal_data %>%
  group_by(River, year, month) %>%
  summarise(mean_salinity = mean(salinity), n_obs = n(), .groups = "drop")

temp_monthly <- temp_data %>%
  group_by(River, year, month) %>%
  summarise(mean_temp = mean(temp), n_obs = n(), .groups = "drop")

env_monthly <- sal_monthly %>%
  full_join(temp_monthly, by = c("River", "year", "month")) %>%
  rename(Year = year, Month = month)

env_averages_all <- data.frame()
for (i in 1:nrow(oyster_collection_dates)) {
  river <- oyster_collection_dates$River[i]
  year <- oyster_collection_dates$Year[i]
  coll_date <- oyster_collection_dates$collection_date[i]
  coll_month <- month(coll_date)
  
  months_to_avg <- c(coll_month - 2, coll_month - 1, coll_month)
  
  env_subset <- env_monthly %>%
    filter(River == river,
           ((Year == year & Month %in% months_to_avg[months_to_avg > 0]) |
              (Year == year - 1 & Month %in% (months_to_avg[months_to_avg <= 0] + 12))))
  
  if (nrow(env_subset) > 0) {
    env_averages_all <- bind_rows(env_averages_all,
                                  data.frame(
                                    River = river,
                                    Year = year,
                                    collection_date = coll_date,
                                    mean_salinity_3mo = mean(env_subset$mean_salinity, na.rm = TRUE),
                                    mean_temp_3mo = mean(env_subset$mean_temp, na.rm = TRUE),
                                    n_months = nrow(env_subset)))}}

env_data_final <- env_quantiles_precollection %>%
  full_join(env_averages_all, by = c("River", "Year", "collection_date"))


# visualizations
river_colors <- c("Severn" = "#6baed6", "South" = "#0c2c84", "West" = "#c7e9b4")

temp_data_summer <- temp_data %>%
  filter(year %in% 2022:2024, month >= 6, month <= 9) %>%
  mutate(Year = factor(year))

sal_data_summer <- sal_data %>%
  filter(year %in% 2022:2024, month >= 6, month <= 9) %>%
  mutate(Year = factor(year))

# first plot let's do a box plot with all temp and sal data
# quartiles, median, outlier points
fig_temp_box <- ggplot(temp_data_summer, aes(x = Year, y = temp, fill = River)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  scale_fill_manual(values = river_colors) +
  labs(x = "Year", y = "Temperature (°C)", title = "") +
  theme_classic(base_size = 16) +
  theme(legend.position = "bottom")

fig_sal_box <- ggplot(sal_data_summer, aes(x = Year, y = salinity, fill = River)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  scale_fill_manual(values = river_colors) +
  labs(x = "Year", y = "Salinity (ppt)", title = "") +
  theme_classic(base_size = 16) +
  theme(legend.position = "bottom")

fig_temporal_box <- fig_temp_box + fig_sal_box
print(fig_temporal_box)

# temp line plot with the point at the mean and error bars showing the q10 and the q90
fig_temp_line <- ggplot(summer_quantiles, aes(x = Year, y = temp_mean, color = River, group = River)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = temp_q10, ymax = temp_q90), width = 0.15, linewidth = 1) +
  scale_color_manual(values = river_colors) +
  scale_x_continuous(breaks = c(2022, 2023, 2024)) +
  labs(x = "Year", y = "Temperature (°C)", title = "") +
  theme_classic(base_size = 16) +
  theme(legend.position = "bottom")

fig_sal_line <- ggplot(summer_quantiles, aes(x = Year, y = sal_mean, color = River, group = River)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = sal_q10, ymax = sal_q90), width = 0.15, linewidth = 1) +
  scale_color_manual(values = river_colors) +
  scale_x_continuous(breaks = c(2022, 2023, 2024)) +
  labs(x = "Year", y = "Salinity (ppt)", title = "") +
  theme_classic(base_size = 16) +
  theme(legend.position = "bottom")

fig_temporal_line <- fig_temp_line + fig_sal_line
print(fig_temporal_line)

# violin plots with ALL temp and sal data
# overlaied point is the summer mean
fig_temp_violin <- ggplot(temp_data_summer, aes(x = Year, y = temp, fill = River)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_point(data = summer_quantiles, aes(x = factor(Year), y = temp_mean, color = River),
             size = 3, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = river_colors) +
  scale_color_manual(values = river_colors) +
  labs(x = "Year", y = "Temperature (°C)", title = "") +
  theme_classic(base_size = 16) +
  theme(legend.position = "bottom")

fig_sal_violin <- ggplot(sal_data_summer, aes(x = Year, y = salinity, fill = River)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_point(data = summer_quantiles, aes(x = factor(Year), y = sal_mean, color = River),
             size = 3, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = river_colors) +
  scale_color_manual(values = river_colors) +
  labs(x = "Year", y = "Salinity (ppt)", title = "") +
  theme_classic(base_size = 16) +
  theme(legend.position = "bottom")

fig_temporal_violin <- fig_temp_violin + fig_sal_violin
print(fig_temporal_violin)

write.csv(env_data_final, "/Users/madelineeppley/Desktop/MarineGEO_oysterdisease/env_data_final.csv")
