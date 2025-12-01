## SERC disease analysis workflow
- I received two disease files from Calli Wise (SERC staff in CDE lab) that contained Dermo and bioeroder data (+ some extras) for the MarineGEO project
- Three years of data (2022, 2023, and 2024) were together in the first file. After the completion of the 2025 data collection (yay!) I was sent a follow-up file with the 2025 data.
- This file contains analysis steps in R that I followed to analyze the data, as well as my thoughts/workflow, which are also heavily annotated in the R scripts
- High-quality PDFs of figures are often generated (alongside png uploads in this .md file) and are stored in the /fig folder of this repo.

## Analysis goals
After our first meeting (2025/10/7), Katrina and I decided to focus on answering these intitial 3 questions with the dataset. Here's the analysis plan: 

| Research Question | Statistical Method | Response Variable | Predictor Variables | Visualization |
|-------------------|-------------------|-------------------|---------------------|---------------|
| How does the relationship between Dermo + bioeroder vary across rivers and over time? | GLM | Bioeroder presence or absence | Dermo intensity * river * year | Interaction plot? |
| How many oysters are co-infected and how does co-infection vary across rivers and over time? | Multinomial logistic regression (?, check) | Infection groups (co-infected, single infection, uninfected) | River, year | Venn diagram with all co-infections and stacked bar plot for each river with year on x |
| Do different rivers have different disease prevalences/intensities? | PCA | Matrix of all infections subset by river (includes dermo intensity, bioeroder presence, and co-infection status) | river | PCA with individual oysters = points colored by river and/or heatmaps of infection prevalence by river |

### Step 1: Figure out what's present/absent in the data
The first thing to do is figure out which rivers are represented in our data, and how many data points are available for each year. This will help me make decisions about what data to keep. 

```R
#### Madeline G Eppley (MGE) data processing for marineGEO oyster disease data
## data source: Lohan Lab @ SERC
## goals: standardize data columns across 4 years of datasets, merge into 1 df, plot heatmaps of diseases
## exports: data for GLM analysis

library(readxl)
library(broom)
library(vegan)
library(nnet)
library(ggplot2)
library(viridis)
library(tidyverse)
library(patchwork)

rawfile <- "/Users/madelineeppley/Desktop/MarineGeo_OysterData_SevernSouthRhodeWest_2025.08.25.xlsx"
rawfile2025 <- "/Users/madelineeppley/Desktop/MarineGeo_OysterData_SevernSouthRhodeWest_2025.11.19.xlsx"

# read in data but just do 2022-2024 since there is no bioeroder data for 2021
data_2022 <- read_excel(rawfile, sheet = "RawData_2022")
data_2023 <- read_excel(rawfile, sheet = "Data2023")
data_2024 <- read_excel(rawfile, sheet = "Data2024")
data_2025 <- read_excel(rawfile2025, sheet = "Data")

# 2022 data has some samples from the Rhode river, but we don't have these for 2023 or 2024
# it looks like there is some Rhode data from 2025, let's check
sum(grepl("Rhode", data_2025$River)) #30 data points
sum(grepl("Rhode", data_2022$River)) #30 data points
sum(grepl("Rhode", data_2023$River)) #0 data points
sum(grepl("Rhode", data_2024$River)) #0 data points

# let's compare that to the number of datapoints we have for other rivers in other years?
sum(grepl("Severn", data_2025$River)) #151
sum(grepl("Severn", data_2024$River)) #31
sum(grepl("Severn", data_2023$River)) #28 
sum(grepl("Severn", data_2022$River)) #30 

# we may have a scaling issue with the number of data points per year .. but for now let's get rid of the Rhode data and we will think about the scaling later
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
```

### Data points per year, per river
The figure generated by the above code script is below. There are clearly have unequal data points per year:

<img width="853" height="853" alt="image" src="https://github.com/user-attachments/assets/6d4c7139-b161-4e0b-bb08-3957ea6fa2e8" />

> [!NOTE]
> For now, I'm not worried about subsetting (2025 data) and will use all available individuals, but this is something to flag to discuss with Katrina. However, I am going to get rid of the Rhode data from 2022 and 2025, since we don't have Rhode data for 2023 and 2024. We should also discuss.


### Standardizing data columns
The data format is pretty messy, I have columns with slightly different names year-to-year, and we also have a good amount of missingness throughout. I should tackle this step-by-step:
- My first decision is to standardize the data to the 2023 format since it is the most complete.
- The 2022 dataframe will need the most changes. Let's make a list of column names in 2022 to change, then use mutate() for efficiency
- Another I notice is that the 2025 data has squashes, which we don't have from any other year.
- I'm also going to want to combine some columns to make new variables for co-infections (Y/N or sum), infected with dermo only, infected with polydora only, and uninfected. we should add these once all column names are standardized.

Here's how I handled some of the initial data wrangling to get all of the columns to have the same names:

```R
### code chunk 2 ###
data_2022 <- data_2022 %>% 
  filter(River != "Rhode")
data_2025 <- data_2025 %>% 
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
```

### Producing bioeroder averages
I noticed that the 2023 data did not average the bioeroder data, so I will need to calculate those averages. During my column filtering above I also saw that the 2022 data has A and B data collectors but not L/R valves, so I duplicated the single value across valves. This will probably get averaged/fixed later anyway. Also one weird R syntax thing that I just learned is that if a variable ends with % you need to encase the variable name in ``. There are also some crazy symbols (e.g., >, ~) in the dataset, so I need to get rid of those with clean_num as as.numeric. 






