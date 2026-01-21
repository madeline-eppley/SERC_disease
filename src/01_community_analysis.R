## parasite community analysis R script
# M.G. Eppley, upload date 01/21/2026

#### Madeline G Eppley (MGE) data processing for marineGEO oyster disease data
## data source: MarineGEO, CDE/Lohan Lab @ SERC
## goals: standardize data columns across 4 years of datasets, merge into 1 df, plot heatmaps of diseases

### code chunk 1 ###
library(car)
library(readxl)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(vegan)
library(broom)
library(dplyr)

setwd("~/Desktop/marineGEO_oysterdisease")

# read in data
rawfile <- "/Users/madelineeppley/Desktop/MarineGeo_OysterData_SevernSouthRhodeWest_2025.11.19.xlsx"
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
# first we want to remove the Rhode River data and the 2021 data.
data <- data %>% 
  filter(River != "Rhode")
data <- data %>% 
  filter(Year != "2021")

# first make a function that we can use to get rid of any symbols in the columns. we need everything to be a clean number.
clean_num <- function(x) {
  as.numeric(str_extract(as.character(x), "\\d+\\.?\\d*"))}

# ok now we have an interesting case in the 2022 cliona data where we have Avg_A and Avg_B values, but not Left Valve/Right Valve
# i'm going to fill in the values for L/R valves with the Avg so that the next function step will work. 
data <- data %>%
  mutate(
    Cliona_LeftValve_A = ifelse(Year == 2022, Cliona_Avg_A, Cliona_LeftValve_A),
    Cliona_RightValve_A = ifelse(Year == 2022, Cliona_Avg_A, Cliona_RightValve_A),
    `Cliona_LeftValve_B_%` = ifelse(Year == 2022, Cliona_Avg_B, `Cliona_LeftValve_B_%`),
    `Cliona_RightValve_B_%` = ifelse(Year == 2022, Cliona_Avg_B, `Cliona_RightValve_B_%`))

# calculate bioeroder averages for all years - should work now that we have our boring sponge valves filled
data <- data %>%
  mutate( # first do the boring sponge averages
    Cliona_LeftValve_A_num = clean_num(Cliona_LeftValve_A),
    Cliona_RightValve_A_num = clean_num(Cliona_RightValve_A),
    Cliona_Avg_A = rowMeans(
      cbind(Cliona_LeftValve_A_num, Cliona_RightValve_A_num),
      na.rm = TRUE),
    Cliona_LeftValve_B_num = clean_num(`Cliona_LeftValve_B_%`),
    Cliona_RightValve_B_num = clean_num(`Cliona_RightValve_B_%`),
    Cliona_Avg_B = rowMeans(
      cbind(Cliona_LeftValve_B_num, Cliona_RightValve_B_num),
      na.rm = TRUE),
    # now the polydora averages
    Polydora_LeftValve_A_num = clean_num(Polydora_LeftValve_A),
    Polydora_RightValve_A_num = clean_num(Polydora_RightValve_A),
    Polydora_Average_A = rowMeans(
      cbind(Polydora_LeftValve_A_num, Polydora_RightValve_A_num),
      na.rm = TRUE),
    Polydora_LeftValve_B_num = clean_num(`Polydora_LeftValve_B_%`),
    Polydora_RightValve_B_num = clean_num(`Polydora_RightValve_B_%`),
    Polydora_Average_B = rowMeans(
      cbind(Polydora_LeftValve_B_num, Polydora_RightValve_B_num),
      na.rm = TRUE))

# keep only the new columns that i want
keep_cols <- c("Specimen_ID", "River", "Year", "Year_Collected", "Month_Collected", "Day_Collected",
               "Score_Rectum", "Count_Rectum", "Score_Mantle", "Count_Mantle",
               "Polydora_LeftValve_A", "Polydora_RightValve_A", "Polydora_Average_A",
               "Polydora_LeftValve_B_%", "Polydora_RightValve_B_%", "Polydora_Average_B",
               "Cliona_LeftValve_A", "Cliona_RightValve_A", "Cliona_Avg_A",
               "Cliona_LeftValve_B_%", "Cliona_RightValve_B_%", "Cliona_Avg_B",
               "Length_mm", "Blister_Count_A", "Blister_Count_B")

all_data <- data[, intersect(names(data), keep_cols)]

# summary stats
dim(all_data) #370 x 25 nice, increased from 278 inds in the 2022-2024 only data
table(all_data$River) # just severn, south, west
table(all_data$Year) # 2022-2025 with fairly even distribution across years

### code chunk 3 ###
# trim any whitespaces around columns (there are some whitespaces, I checked lol)
all_data$River <- trimws(all_data$River)
all_data$Specimen_ID <- trimws(all_data$Specimen_ID)

# now get everything into numeric format
all_data$Length_mm <- as.numeric(as.character(all_data$Length_mm))
all_data$Score_Rectum <- as.numeric(as.character(all_data$Score_Rectum))
all_data$Score_Mantle <- as.numeric(as.character(all_data$Score_Mantle))
all_data$Blister_Count_A <- as.numeric(as.character(all_data$Blister_Count_A))
all_data$Blister_Count_B <- as.numeric(as.character(all_data$Blister_Count_B))

# a little sanity check with the ranges of our new numeric variables
range(all_data$Length_mm, na.rm = TRUE)
quantile(all_data$Score_Rectum, na.rm = TRUE) # just taking a quick peek at the data, looks like rectum is more infected than mantle
quantile(all_data$Score_Mantle, na.rm = TRUE)

# now average the A and B observer values
all_data$Polydora_outside <- (all_data$Polydora_Average_A + all_data$Polydora_Average_B) / 2
all_data$Cliona <- (all_data$Cliona_Avg_A + all_data$Cliona_Avg_B) / 2
all_data$Polydora_blister <- (all_data$Blister_Count_A + all_data$Blister_Count_B) / 2

# remove the old cols
all_data <- all_data %>%
  select(-c(Polydora_LeftValve_A, Polydora_RightValve_A, Polydora_Average_A, Polydora_Average_B,
            `Polydora_LeftValve_B_%`, `Polydora_RightValve_B_%`,
            Cliona_LeftValve_A, Cliona_RightValve_A, `Cliona_LeftValve_B_%`, `Cliona_RightValve_B_%`,
            Blister_Count_A, Blister_Count_B, Count_Rectum, Count_Mantle,
            Cliona_Avg_A, Cliona_Avg_B, Year_Collected))

# summary checks again
summary(all_data)
dim(all_data)

# ok now what i'm noticing is that we have a handful of individuals (anywhere from 1-7 in each col) that still have NAs
# i think the best way to deal with this data is to remove these individuals
# katrina should confirm that we want to do this - i don't think they'll be helpful data points down the line

# Create infection status variables
all_data$dermo_infected <- ifelse(all_data$Score_Rectum > 0, 1, 0)
all_data$polydora_infected <- ifelse(all_data$Polydora_outside > 0, 1, 0)
all_data$cliona_infected <- ifelse(all_data$Cliona > 0, 1, 0)
all_data$n_infections <- all_data$dermo_infected + all_data$polydora_infected + all_data$cliona_infected

# Check for individuals with NA infections and investigate
na_infections <- all_data[is.na(all_data$n_infections), ]
nrow(na_infections) # only 8 currently - let's just remove those!

# remove those inds
all_data <- all_data[!is.na(all_data$n_infections), ]
dim(all_data) # we just removed 8, this is now down to 362, looks good to me!

# ok now make a variable for our co-infected/uninfected/single infections
# i think katrina and I should also talk about this, how do we want to document each possible combination of co-infections?
all_data$infection_group <- ifelse(all_data$n_infections == 0, "Uninfected",
                                   ifelse(all_data$n_infections == 1, "Single", "Co-infected"))
all_data$infection_group <- factor(all_data$infection_group, levels = c("Uninfected", "Single", "Co-infected"))

# ok again more summary info
table(all_data$infection_group) # wow only 1 individual that was completely uninfected! 
table(all_data$River, all_data$infection_group)
table(all_data$Year, all_data$infection_group)

# export the shiny new clean data!
#write.csv(all_data, "/Users/madelineeppley/Desktop/marineGEO_oysterdisease/marineGEO_cleanv20260113.csv", row.names = FALSE)
#all_data <- read.csv("/Users/madelineeppley/Desktop/marineGEO_oysterdisease/marineGEO_cleanv20260113.csv")
str(all_data)

### code chunk 4 ###
# visualize infections!!!
infection_counts <- table(all_data$Year, all_data$River, all_data$infection_group)
layout(matrix(c(1, 2, 3, 4), nrow = 1), widths = c(1, 1, 1, 0.4))
par(mar = c(5, 4, 4, 1))
barplot(t(infection_counts[, "Severn", ]), beside = FALSE, col = c("#08519c","#6baed6","#eff3ff"),
        main = "Severn", xlab = "Year", ylab = "Individuals")
barplot(t(infection_counts[, "South", ]), beside = FALSE, col = c("#08519c","#6baed6","#eff3ff"),
        main = "South", xlab = "Year", ylab = "Individuals")
barplot(t(infection_counts[, "West", ]), beside = FALSE, col = c("#08519c","#6baed6","#eff3ff"),
        main = "West", xlab = "Year", ylab = "Individuals")
par(mar = c(0, 0, 0, 0))
plot.new()
legend("left", legend = c("Co-infected", "Single", "Uninfected"), 
       fill = c("#08519c","#6baed6","#eff3ff"), bty = "n", cex = 1)
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))


#pca analysis
# let's include dermo rectum, dermo mantle, polydora outside, cliona, and polydora blister
pca_matrix <- cbind(all_data$Score_Rectum, all_data$Score_Mantle, 
                    all_data$Polydora_outside, all_data$Cliona, all_data$Polydora_blister)
colnames(pca_matrix) <- c("Dermo_Rectum", "Dermo_Mantle", "Polydora", "Cliona", "Blister")

pca_result <- prcomp(pca_matrix, scale. = TRUE, center = TRUE) # scaled and centered
summary(pca_result)

#dev.off()
# pca by river
river_colors <- c("Severn" = "#0c2c84", "South" = "#c7e9b4", "West" = "#41b6c4")
year_shapes <- c("2022" = 16, "2023" = 17, "2024" = 15, "2025" = 18)
point_shapes <- year_shapes[as.character(all_data$Year)]
plot(pca_result$x[, 1], pca_result$x[, 2], 
     col = river_colors[all_data$River], 
     pch = point_shapes, cex = 1.5,
     xlab = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
     ylab = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)"),
     main = "")
legend("topleft", 
       legend = c("Severn", "South", "West", "2022", "2023", "2024", "2025"),
       y.intersp = 0.75,
       col = c(river_colors, "black", "black", "black", "black"), 
       pch = c(16, 16, 16, 16, 17, 15, 18), cex = 1)

# add in ellipses for rivers
xlim_range <- range(pca_result$x[, 1])
ylim_range <- range(pca_result$x[, 2])
x_padding <- diff(xlim_range) * 0.15
y_padding <- diff(ylim_range) * 0.15

plot(pca_result$x[, 1], pca_result$x[, 2], 
     col = river_colors[all_data$River], 
     pch = point_shapes, cex = 1.5,
     xlab = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
     ylab = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)"),
     main = "PCA - River", 
     xlim = c(xlim_range[1] - x_padding, xlim_range[2] + x_padding),
     ylim = c(ylim_range[1] - y_padding, ylim_range[2] + y_padding))
for(river in names(river_colors)) {
  river_data <- pca_result$x[all_data$River == river, 1:2]
  ordiellipse(river_data, groups = rep(1, nrow(river_data)), 
              kind = "sd", conf = 0.95, 
              col = river_colors[river], lwd = 2, draw = "polygon",
              alpha = 50, border = river_colors[river])}
legend("topleft", 
       legend = c("Severn", "South", "West", "2022", "2023", "2024", "2025"), 
       col = c(river_colors, "black", "black", "black", "black"),
       y.intersp = 0.75,
       pch = c(16, 16, 16, 16, 17, 15, 18), cex = 1)


# now add in ellipses per year
year_colors <- c("2022" = "#6baed6", "2023" = "#0c2c84", "2024" = "#c7e9b4", "2025" = "#08519c")

plot(pca_result$x[, 1], pca_result$x[, 2], 
     col = year_colors[as.character(all_data$Year)], 
     pch = 16, cex = 1.5,
     xlab = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
     ylab = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)"),
     main = "PCA - Year", 
xlim = c(xlim_range[1] - x_padding, xlim_range[2] + x_padding),
ylim = c(ylim_range[1] - y_padding, ylim_range[2] + y_padding))
for(year in names(year_colors)) {
  year_data <- pca_result$x[as.character(all_data$Year) == year, 1:2]
  if(nrow(year_data) > 2) {  # Need at least 3 points for ellipse
    ordiellipse(year_data, groups = rep(1, nrow(year_data)), 
                kind = "sd", conf = 0.95, 
                col = year_colors[year], lwd = 2, draw = "polygon",
                alpha = 50, border = year_colors[year])}}
legend("topleft", 
       legend = c("2022", "2023", "2024", "2025"),
       y.intersp = 0.75,
       col = year_colors, 
       pch = 16, cex = 1)

### understanding the data more ####
### now we want to achieve 3 things - rectum vs mantle intensity, composition of the co-infections, and the average length of the oysters in the study. 

## create categories for co-infection data
all_data <- all_data %>%
  mutate(
    infection_combo = case_when(
      dermo_infected == 0 & polydora_infected == 0 & cliona_infected == 0 ~ "None",
      dermo_infected == 1 & polydora_infected == 0 & cliona_infected == 0 ~ "Dermo only",
      dermo_infected == 0 & polydora_infected == 1 & cliona_infected == 0 ~ "Polydora only",
      dermo_infected == 0 & polydora_infected == 0 & cliona_infected == 1 ~ "Cliona only",
      dermo_infected == 1 & polydora_infected == 1 & cliona_infected == 0 ~ "Dermo + Polydora",
      dermo_infected == 1 & polydora_infected == 0 & cliona_infected == 1 ~ "Dermo + Cliona",
      dermo_infected == 0 & polydora_infected == 1 & cliona_infected == 1 ~ "Polydora + Cliona",
      dermo_infected == 1 & polydora_infected == 1 & cliona_infected == 1 ~ "Dermo + Polydora + Cliona"))
# now check the distributions
table(all_data$infection_combo)

# let's make a coinfection summary by river and year
coinf_summary <- all_data %>%
  group_by(River, Year, infection_combo) %>%
  summarise(n = n(), .groups = "drop")

# and now also proportions of each type of coinfection per river and per year
coinf_prop <- coinf_summary %>%
  group_by(River, Year) %>%
  mutate(prop = n / sum(n))

# plot!
ggplot(coinf_prop, aes(x = factor(Year), y = prop, fill = infection_combo)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ River) +
  labs(
    x = "Year",
    y = "Proportion of samples",
    fill = "Infections",
    title = "Co-infection combinations by river and year") +
  theme_minimal()

## improve the plot
coinf_colors <- c(
  "None"                      = "lightyellow",
  "Polydora only"             = "#c7e9c0",
  "Cliona only"               = "#6baed6",
  #"Dermo only"                = "#c6dbef",
  "Dermo + Polydora"          = "#08519c",
  #"Dermo + Cliona"            = "#74c476",
  "Polydora + Cliona"         = "#74c476",
  "Dermo + Polydora + Cliona" = "#9ecae1")

ggplot(coinf_prop, aes(x = factor(Year), y = prop, fill = infection_combo)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ River) +
  scale_fill_manual(values = coinf_colors, drop = FALSE) +
  labs(
    x = "Year",
    y = "Proportion of samples",
    fill = "Infections",
    title = "Co-infection groups") +
  theme_minimal() +
  theme(
    legend.position = "bottom")


### now the rectum vs mantle comparisons
dermo_long <- all_data %>%
  select(River, Year, Score_Rectum, Score_Mantle) %>%
  pivot_longer(
    cols = c(Score_Rectum, Score_Mantle),
    names_to = "Tissue",
    values_to = "Dermo_Score") %>%
  mutate(
    Tissue = recode(Tissue,
                    Score_Rectum = "Rectum",
                    Score_Mantle = "Mantle"))

tissue_colors <- c(
  "Rectum" = "#6baed6",
  "Mantle" = "#31a354")

ggplot(dermo_long, aes(x = Tissue, y = Dermo_Score, fill = Tissue)) +
  geom_boxplot(outlier.alpha = 0.4, width = 0.6) +
  facet_grid(River ~ Year) +
  scale_fill_manual(values = tissue_colors) +
  labs(
    x = "",
    y = "Dermo score",
    title = "Dermo score in rectum vs mantle") +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 9))

# box plots for each year and river
ggplot(dermo_long, aes(x = Tissue, y = Dermo_Score, fill = Tissue)) +
  geom_boxplot(outlier.alpha = 0.4) +
  facet_grid(River ~ Year) +
  labs(
    x = "",
    y = "Dermo score",
    title = "Dermo Score in Rectum vs Mantle") +
  theme_minimal() +
  theme(legend.position = "none")

# i kind of want to see this visualized as the difference between avg rectum and mantle?
# i think by subtracting we can get a value that will characterize more infected mantle or rectum.
all_data <- all_data %>%
  mutate(
    dermo_diff = Score_Rectum - Score_Mantle)

summary(all_data$dermo_diff)
dermo_diff_summary <- all_data %>%
  group_by(River, Year) %>%
  summarise(
    mean_diff = mean(dermo_diff, na.rm = TRUE),
    sd_diff   = sd(dermo_diff, na.rm = TRUE),
    n         = n(),
    .groups = "drop")

# here, values above 0 = rectum score is heavier than mantle score
## it looks like we have heavier rectum scores than mantle scores, so these are early infections.
# scores near 0 would be an even infection of rectum & mantle
# scores below 0 would be a heavier infection of mantle than rectum.
river_colors <- c(
  "Severn" = "#6baed6",
  "South"  = "#0c2c84",
  "West"   = "#c7e9b4")

ggplot(dermo_diff_summary,
       aes(x = factor(Year), y = mean_diff,
           color = River, group = River)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = river_colors) +
  labs(
    x = "Year",
    y = "Mean Dermo (rectum − mantle) score",
    title = "Dermo (rectum - mantle) score") +
  theme_minimal()

## ok last step - now let's make a plot of the oyster sizes
# use the same river colors nad plot year as each individual panel
ggplot(all_data, aes(x = Length_mm, fill = River, color = River)) +
  geom_histogram(
    bins = 30,
    alpha = 0.55,
    position = "identity") +
  facet_wrap(~ Year) +
  scale_fill_manual(values = river_colors) +
  scale_color_manual(values = river_colors) +
  labs(
    x = "Shell length (mm)",
    y = "Count",
    title = "Oyster shell length distributions") +
  theme_minimal() +
  theme(
    legend.position = "bottom")

# maybe also a box plot that has average shell length
ggplot(all_data, aes(x = River, y = Length_mm, fill = River)) +
  geom_boxplot(alpha = 0.75, width = 0.6, outlier.alpha = 0.4) +
  facet_wrap(~ Year) +
  scale_fill_manual(values = river_colors) +
  labs(
    x = "",
    y = "Shell length (mm)",
    title = "Average oyster shell length") +
  theme_minimal() +
  theme(
    legend.position = "none")

### code chunk 5 - model prep ###
# ultimate goal is to understand if the individuals that are co-infected with dermo & polydora are more intensely infected
# subset data to our co-infected dermo and polydora individuals
coinfected_data <- all_data %>%
  filter(dermo_infected == 1 & polydora_infected == 1)
nrow(coinfected_data) # 276 co-infected with dermo and polydora

# quick look at the intensity metrics
summary(coinfected_data$Score_Rectum)
summary(coinfected_data$Polydora_outside)
summary(coinfected_data$Polydora_blister)

# check for correlations between the variables
cor.test(coinfected_data$Score_Rectum, coinfected_data$Polydora_outside) # no correlation
cor.test(coinfected_data$Score_Rectum, coinfected_data$Polydora_blister) # weak positive
cor.test(coinfected_data$Polydora_outside, coinfected_data$Polydora_blister) # moderate positive r=0.38 (makes sense, these are both polydora)

# check the variance inflation factor for the polydora variables
model_both <- lm(Score_Rectum ~ Polydora_outside + Polydora_blister + River + Year, 
                 data = coinfected_data)
vif(model_both) # these look fine probably, ~1-1.2

# pull out our main variables of interest
score_rec <- coinfected_data$Score_Rectum
poly_blister <- coinfected_data$Polydora_blister
