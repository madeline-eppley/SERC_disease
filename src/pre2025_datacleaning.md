# preserving an old version of the analysis from pre-combined file (2022, 2023, and 2024 were separate tabs in the August 2025 version of the datasheet). 

```R
#### Madeline G Eppley (MGE) data processing for marineGEO oyster disease data
## data source: Lohan Lab @ SERC
## goals: standardize data columns across 4 years of datasets, merge into 1 df, plot heatmaps of diseases
## exports: data for GLM analysis

### code chunk 1 ###
library(readxl)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(vegan)

library(broom)
library(nnet)
library(viridis)

#rawfile <- "/Users/madelineeppley/Desktop/MarineGeo_OysterData_SevernSouthRhodeWest_2025.08.25.xlsx"
rawfile <- "/Users/madelineeppley/Desktop/MarineGeo_OysterData_SevernSouthRhodeWest_2025.11.19.xlsx"

# read in data but just do 2022-2024 since there is no bioeroder data for 2021
#data_2022 <- read_excel(rawfile, sheet = "RawData_2022")
#data_2023 <- read_excel(rawfile, sheet = "Data2023")
#data_2024 <- read_excel(rawfile, sheet = "Data2024")
data <- read_excel(rawfile, sheet = "Data") # this is now all years of data together, 485 x 65

colnames(data) # going to rename a few of the data columns to keep my code functionality that I wrote pre-2025 dataset
data <- data %>%
  mutate(
    Score_Rectum = DermoScore_Rectum,
    Count_Rectum = DermoCellCount_Rectum,
    Score_Mantle = DermoScore_Mantle,
    Count_Mantle = DermoCellCount_Mantle)

# 2022 data has some samples from the Rhode river, but we don't have these for 2023, 2024, or 2025
# let's check by making subsets
data_2022 <- data %>% filter(Year == 2022)
data_2023 <- data %>% filter(Year == 2023)
data_2024 <- data %>% filter(Year == 2024)
data_2025 <- data %>% filter(Year == 2025)

sum(grepl("Rhode", data_2022$River)) #30 data points
sum(grepl("Rhode", data_2023$River)) #0 data points
sum(grepl("Rhode", data_2024$River)) #0 data points
sum(grepl("Rhode", data_2025$River)) #0 data points

# now let's take a quick look at how many data points we have per year.
counts_2022 <- as.data.frame(table(data_2022$River))
colnames(counts_2022) <- c("River", "Count")

counts_2023 <- as.data.frame(table(data_2023$River))
colnames(counts_2023) <- c("River", "Count")

counts_2024 <- as.data.frame(table(data_2024$River))
colnames(counts_2024) <- c("River", "Count")

counts_2025 <- as.data.frame(table(data_2025$River))
colnames(counts_2025) <- c("River", "Count")

# let's make one composite figure
p22 <- ggplot(counts_2022, aes(x = River, y = Count)) +
  geom_bar(stat = "identity", fill = "#6baed6") +
  labs(title = "2022",
       x = "River",
       y = "Individuals") +
  theme_minimal()

p23 <- ggplot(counts_2023, aes(x = River, y = Count)) +
  geom_bar(stat = "identity", fill = "#0c2c84") +
  labs(title = "2023",
       x = "River",
       y = "Individuals") +
  theme_minimal()

p24 <- ggplot(counts_2024, aes(x = River, y = Count)) +
  geom_bar(stat = "identity", fill = "#c7e9b4") +
  labs(title = "2024",
       x = "River",
       y = "Individuals") +
  theme_minimal()

p25 <- ggplot(counts_2025, aes(x = River, y = Count)) +
  geom_bar(stat = "identity", fill = "mediumblue") +
  labs(title = "2025",
       x = "River",
       y = "Individuals") +
  theme_minimal()

# arrange plot
(p22 | p23) / (p24 | p25)

### code chunk 2 ###
data_2022 <- data_2022 %>% 
  filter(River != "Rhode")

# do some renaming for 2022 since the col names don't quite match - use 2023/2024 col names
data_2022 <- data_2022 %>%
  mutate(
    Score_Rectum = DermoScore_Rectum,
    Count_Rectum = DermoCellCount_Rectum,
    Score_Mantle = DermoScore_Mantle,
    Count_Mantle = DermoCellCount_Mantle,
    Initials_BioEroder_A = Initials_BE_A,
    Initials_BioEroder_B = Initials_BE_B,
    # this dataset is weird because it only has A and B but not L/R valve so we will just duplicate everything across valves
    Cliona_LeftValve_A = Cliona_A,
    Cliona_RightValve_A = Cliona_A,
    `Cliona_LeftValve_B_%` = Cliona_B,
    `Cliona_RightValve_B_%` = Cliona_B,
    Cliona_Avg_A = Cliona_A,
    Cliona_Avg_B = Cliona_B,
    `Polydora_LeftValve_B_%` = Polydora_LeftValve_B,
    `Polydora_RightValve_B_%` = Polydora_RightValve_B
  ) %>%
  select(-DermoScore_Rectum, -DermoCellCount_Rectum, -DermoScore_Mantle, -DermoCellCount_Mantle,
         -Initials_BE_A, -Initials_BE_B, -Cliona_A, -Cliona_B,
         -Polydora_LeftValve_B, -Polydora_RightValve_B) # remove the old cols now to not get confused

# 2024 needs to have the first column renamed
data_2024 <- data_2024 %>%
  rename(Specimen_ID = free)

# i just standardized to the 2023/2024 data so just re-store 2023 here
data_2023 <- data_2023

# now just need to update a few columns in 2025
data_2025 <- data_2025 %>%
  mutate(
    Score_Rectum = DermoScore_Rectum,
    Count_Rectum = DermoCellCount_Rectum,
    Score_Mantle = DermoScore_Mantle,
    Count_Mantle = DermoCellCount_Mantle
  ) %>%
  select(-DermoScore_Rectum, -DermoCellCount_Rectum, 
         -DermoScore_Mantle, -DermoCellCount_Mantle)

# check
colnames(data_2022)
colnames(data_2023)
colnames(data_2024)
colnames(data_2025)
# just a few small differences in notes but we will filter those out anyway

### code chunk 3 ###
## noticed here that our 2023 data did not produce bioeroder averages, so we need to make those
clean_num <- function(x) {
  as.numeric(str_extract(x, "\\d+"))
} # first make a function that gets rid of the characters ">" etc. 

data_2023 <- data_2023 %>%
  mutate(
    Cliona_LeftValve_A_num  = clean_num(Cliona_LeftValve_A),
    Cliona_RightValve_A_num = as.numeric(Cliona_RightValve_A),
    Cliona_Avg_A = rowMeans(
      cbind(Cliona_LeftValve_A_num, Cliona_RightValve_A_num),
      na.rm = TRUE),
    Cliona_LeftValve_B_num  = clean_num(`Cliona_LeftValve_B_%`),
    Cliona_RightValve_B_num = clean_num(`Cliona_RightValve_B_%`),
    Cliona_Avg_B = rowMeans(
      cbind(Cliona_LeftValve_B_num, Cliona_RightValve_B_num),
      na.rm = TRUE),
    Polydora_LeftValve_A_num  = clean_num(Polydora_LeftValve_A),
    Polydora_RightValve_A_num = clean_num(Polydora_RightValve_A),
    Polydora_Average_A = rowMeans(
      cbind(Polydora_LeftValve_A_num, Polydora_RightValve_A_num),
      na.rm = TRUE),
    Polydora_LeftValve_B_num  = clean_num(`Polydora_LeftValve_B_%`),
    Polydora_RightValve_B_num = `Polydora_RightValve_B_%`,
    Polydora_Average_B = rowMeans(
      cbind(Polydora_LeftValve_B_num, Polydora_RightValve_B_num),
      na.rm = TRUE))

# what are the common cols that we want to keep?
# ID, river, year, month, dermo (rectum and mantle), polydora, cliona, oyster length
keep_cols <- c("Specimen_ID", "River", "Year", "Year_Collected", "Month_Collected", "Day_Collected", 
               "Score_Rectum", "Count_Rectum", "Score_Mantle", "Count_Mantle",
               "Polydora_LeftValve_A", "Polydora_RightValve_A", "Polydora_Average_A",
               "Polydora_LeftValve_B_%", "Polydora_RightValve_B_%", "Polydora_Average_B",
               "Cliona_LeftValve_A", "Cliona_RightValve_A", "Cliona_Avg_A",
               "Cliona_LeftValve_B_%", "Cliona_RightValve_B_%", "Cliona_Avg_B",
               "Length_mm", "Blister_Count_A", "Blister_Count_B")

# now we can subset and merge into one df
data_2022_subset <- data_2022[, intersect(names(data_2022), keep_cols)]
data_2023_subset <- data_2023[, intersect(names(data_2023), keep_cols)]
data_2024_subset <- data_2024[, intersect(names(data_2024), keep_cols)]
colnames(data_2022_subset)
dim(data_2022_subset) # 90 inds x 25 cols
colnames(data_2023_subset)
dim(data_2023_subset) # 88 inds x 25 cols
colnames(data_2024_subset)
dim(data_2024_subset)  # 100 inds x 25 cols

# ok all is looking good to go, now let's see if we can rbind
all_data <- rbind.data.frame(data_2022_subset, data_2023_subset, data_2024_subset)
dim(all_data) # nice! 278 x 25 is all of our inds

# let's check the distribution of some data points
table(all_data$River) # we have 89 Severn, 99 South, 90 West
table(all_data$Year) # we havd 90 2022, 88 2023, and 100 2024 samples 
range(all_data$Length_mm) # not a number yet, we need to do some data cleaning

#write.csv(all_data, "/Users/madelineeppley/Desktop/marineGEO_oysterdisease/all_data.csv")

#### read data back in here now that we have a clean exported copy
all_data <- read.csv("/Users/madelineeppley/Desktop/marineGEO_oysterdisease/all_data.csv")

### data cleaning ####
## lots we need to do here! 
all_data$River <- trimws(all_data$River)
all_data$Specimen_ID <- trimws(all_data$Specimen_ID)
all_data$Length_mm <- as.numeric(all_data$Length_mm)
range(all_data$Length_mm) # not working here we have NA/NA

all_data$Score_Rectum <- as.numeric(as.character(all_data$Score_Rectum))
all_data$Score_Mantle <- as.numeric(as.character(all_data$Score_Mantle))
quantile(all_data$Score_Rectum) # just starting to look at values here
quantile(all_data$Score_Mantle) # mantle less intensely infected than rectum

# average the A and B observer values and combine into new col
all_data$Polydora_outside <- (all_data$Polydora_Average_A + all_data$Polydora_Average_B) / 2
all_data$Cliona <- (all_data$Cliona_Avg_A + all_data$Cliona_Avg_B) / 2
all_data$Polydora_blister <- (all_data$Blister_Count_A + all_data$Blister_Count_B) / 2

# now get rid of the old columns to simplify df
all_data$Polydora_LeftValve_A <- NULL
all_data$Polydora_RightValve_A <- NULL
all_data$Polydora_Average_A <- NULL
all_data$Polydora_Average_B <- NULL
all_data$Polydora_LeftValve_B_. <- NULL
all_data$Polydora_RightValve_B_. <- NULL
all_data$Cliona_LeftValve_A <- NULL
all_data$Cliona_RightValve_A <- NULL
all_data$Cliona_LeftValve_B_. <- NULL
all_data$Cliona_RightValve_B_. <- NULL
all_data$Blister_Count_A <- NULL
all_data$Blister_Count_B <- NULL
all_data$Count_Rectum <- NULL ## so for now I'm going to take this out, not sure how to analyze raw counts
all_data$Count_Mantle <- NULL ## same for here for raw counts remove for now
all_data$Cliona_Avg_A <- NULL
all_data$Cliona_Avg_B <- NULL
all_data$Year_Collected <- NULL
all_data$X <- NULL

#### plot heatmap
disease_by_river_year <- aggregate(
  cbind(Score_Rectum, Score_Mantle, Polydora_outside, Cliona, Polydora_blister) ~ River + Year,
  data = all_data,
  FUN = mean,
  na.rm = TRUE)

str(disease_by_river_year)

all_data$dermo_infected <- ifelse(all_data$Score_Rectum > 0, 1, 0)
all_data$polydora_infected <- ifelse(all_data$Polydora_outside > 0, 1, 0)
all_data$cliona_infected <- ifelse(all_data$Cliona > 0, 1, 0)
all_data$n_infections <- all_data$dermo_infected + all_data$polydora_infected + all_data$cliona_infected

# we have one individual at this point that has infection NAs, remove
all_data <- all_data[!(all_data$Specimen_ID == "MG_0335"), ]

# make infection groups for single infection vs coinfection
all_data$infection_group <- ifelse(all_data$n_infections == 0, "Uninfected",
                                   ifelse(all_data$n_infections == 1, "Single", "Co-infected"))
all_data$infection_group <- factor(all_data$infection_group, levels = c("Uninfected", "Single", "Co-infected"))

# summary info
table(all_data$infection_group)
table(all_data$River, all_data$infection_group)
table(all_data$Year, all_data$infection_group)

# visualize infections
infection_counts <- table(all_data$Year, all_data$River, all_data$infection_group)
layout(matrix(c(1, 2, 3, 4), nrow = 1), widths = c(1, 1, 1, 0.4))
par(mar = c(5, 4, 4, 1))
barplot(t(infection_counts[, "Severn", ]), beside = FALSE, col = c("#eff3ff", "#6baed6", "#08519c"),
        main = "Severn", xlab = "Year", ylab = "Count of Oysters")
barplot(t(infection_counts[, "South", ]), beside = FALSE, col = c("#eff3ff", "#6baed6", "#08519c"),
        main = "South", xlab = "Year", ylab = "Count of Oysters")
barplot(t(infection_counts[, "West", ]), beside = FALSE, col = c("#eff3ff", "#6baed6", "#08519c"),
        main = "West", xlab = "Year", ylab = "Count of Oysters")
par(mar = c(0, 0, 0, 0))
plot.new()
legend("left", legend = c("Uninfected", "Single", "Co-infected"), 
       fill = c("#eff3ff", "#6baed6", "#08519c"), bty = "n", cex = 1.2)
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))


#pca analysis
# let's include dermo rectum, dermo mantle, polydora outside, cliona, and polydora blister
pca_matrix <- cbind(all_data$Score_Rectum, all_data$Score_Mantle, 
                    all_data$Polydora_outside, all_data$Cliona, all_data$Polydora_blister)
colnames(pca_matrix) <- c("Dermo_Rectum", "Dermo_Mantle", "Polydora", "Cliona", "Blister")

pca_result <- prcomp(pca_matrix, scale. = TRUE, center = TRUE)
summary(pca_result)

#dev.off()
# pca by river
river_colors <- c("Severn" = "#0c2c84", "South" = "#c7e9b4", "West" = "#41b6c4")
plot(pca_result$x[, 1], pca_result$x[, 2], col = river_colors[all_data$River], 
     pch = ifelse(all_data$Year == 2022, 16, 17), cex = 1.5,
     xlab = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
     ylab = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)"),
     main = "PCA")
legend("topleft", legend = c("Severn", "South", "West", "2022", "2024"), 
       col = c(river_colors, "black", "black"), pch = c(16, 16, 16, 16, 17), cex = 0.8)

```
