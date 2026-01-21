## R code for ordinal regression modeling
### M. Eppley 1/21/2026

# subset to co-infected individuals (dermo + polydora)
coinfected_data <- all_data %>%
  filter(dermo_infected == 1 & polydora_infected == 1)

# prepare variables for ordinal models
coinfected_data <- coinfected_data %>%
  mutate(
    Score_Rectum_ord = factor(Score_Rectum, levels = c(0.5, 1, 2, 3, 4, 5), ordered = TRUE),
    River = factor(River),
    Year = factor(Year),
    DermoIntensity = as.factor(Score_Rectum))

# check sample sizes
table(coinfected_data$Score_Rectum)
table(coinfected_data$River)
table(coinfected_data$Year)

# check correlations between polydora metrics
cor.test(coinfected_data$Score_Rectum, coinfected_data$Polydora_outside)
cor.test(coinfected_data$Score_Rectum, coinfected_data$Polydora_blister)
cor.test(coinfected_data$Polydora_outside, coinfected_data$Polydora_blister)

## compare candidate models
model_formulas <- list(
  M1 = Score_Rectum_ord ~ Polydora_blister + Year + River,
  M2 = Score_Rectum_ord ~ Polydora_outside + Year + River,
  M3 = Score_Rectum_ord ~ Polydora_blister + River,
  M4 = Score_Rectum_ord ~ Polydora_blister + Year,
  M5 = Score_Rectum_ord ~ Polydora_blister,
  M6 = Score_Rectum_ord ~ Polydora_blister * Year + River,
  M7 = Score_Rectum_ord ~ Polydora_blister * River + Year,
  M8 = Score_Rectum_ord ~ Polydora_blister + Polydora_outside + Year + River,
  M9 = Score_Rectum_ord ~ Polydora_blister * Polydora_outside + Year + River)

model_list <- list()
for (name in names(model_formulas)) {
  fmla <- model_formulas[[name]]
  model <- try(clm(fmla, data = coinfected_data), silent = TRUE)
  model_list[[name]] <- model}

# compare AICs
model_aics <- sapply(model_list, function(m) {
  if (inherits(m, "try-error")) return(NA)
  return(AIC(m))})
sort(model_aics)

# best model: M2 (shell cover + Year + River)
model_M2 <- clm(Score_Rectum_ord ~ Polydora_outside + Year + River, data = coinfected_data)
summary(model_M2)

# test proportional odds assumption
nominal_test(model_M2)

# partial proportional odds model allowing polydora to vary across thresholds
model_final <- clm(
  Score_Rectum_ord ~ Polydora_outside + Year + River,
  nominal = ~ Polydora_outside,
  data = coinfected_data)
summary(model_final)
AIC(model_final)

# model comparison table
model_comparison_table <- data.frame(
  Model = c("Shell cover + Year + River", "Both metrics + Year + River", "Blister count + Year + River"),
  AIC = c(AIC(model_M2), AIC(model_list$M8), AIC(model_list$M1)),
  Polydora_p = c(
    summary(model_M2)$coefficients["Polydora_outside", "Pr(>|z|)"],
    summary(model_list$M8)$coefficients["Polydora_outside", "Pr(>|z|)"],
    summary(model_list$M1)$coefficients["Polydora_blister", "Pr(>|z|)"])) %>%
  mutate(
    Delta_AIC = AIC - min(AIC),
    AIC = round(AIC, 1),
    Delta_AIC = round(Delta_AIC, 1),
    Polydora_p = round(Polydora_p, 3)) %>%
  arrange(AIC)

print(model_comparison_table)

## model viz
# mean polydora metrics by dermo intensity
means_by_dermo <- coinfected_data %>%
  group_by(Score_Rectum) %>%
  summarise(
    n = n(),
    mean_outside = mean(Polydora_outside, na.rm = TRUE),
    se_outside = sd(Polydora_outside, na.rm = TRUE) / sqrt(n()),
    mean_blister = mean(Polydora_blister, na.rm = TRUE),
    se_blister = sd(Polydora_blister, na.rm = TRUE) / sqrt(n()))

# pairwise comparisons with kruskal-Wallis and wilcoxon tests
kruskal.test(Polydora_outside ~ Score_Rectum, data = coinfected_data)
kruskal.test(Polydora_blister ~ Score_Rectum, data = coinfected_data)

pairwise_outside <- pairwise.wilcox.test(coinfected_data$Polydora_outside,
                                         coinfected_data$DermoIntensity,
                                         p.adjust.method = "BH")
pairwise_blister <- pairwise.wilcox.test(coinfected_data$Polydora_blister,
                                         coinfected_data$DermoIntensity,
                                         p.adjust.method = "BH")

# compact letter display for sig groups
library(multcompView)
letters_outside <- multcompLetters(pairwise_outside$p.value)
letters_blister <- multcompLetters(pairwise_blister$p.value)

cld_outside <- data.frame(
  DermoIntensity = c("0.5", "1", "2", "3", "4", "5"),
  letter = c("a", "a", "ab", "b", "ab", "ab"))
cld_outside <- left_join(cld_outside,
                         means_by_dermo %>% mutate(DermoIntensity = as.character(Score_Rectum)),
                         by = "DermoIntensity")

cld_blister <- data.frame(
  DermoIntensity = c("0.5", "1", "2", "3", "4", "5"),
  letter = c("a", "a", "a", "b", "ab", "ab"))
cld_blister <- left_join(cld_blister,
                         means_by_dermo %>% mutate(DermoIntensity = as.character(Score_Rectum)),
                         by = "DermoIntensity")

# figure with polydora shell cover vs dermo intensity
fig1_left <- ggplot(means_by_dermo, aes(x = Score_Rectum, y = mean_outside)) +
  geom_line(color = "black", linewidth = 1.2) +
  geom_point(size = 4, color = "black") +
  geom_errorbar(aes(ymin = mean_outside - se_outside, ymax = mean_outside + se_outside),
                width = 0.2, color = "black", linewidth = 0.8) +
  geom_text(data = cld_outside,
            aes(x = as.numeric(DermoIntensity), y = mean_outside + se_outside + 3, label = letter),
            size = 5, fontface = "bold") +
  scale_x_continuous(breaks = c(0.5, 1, 2, 3, 4, 5)) +
  labs(x = "Dermo Intensity", y = "Polydora Shell Cover (%)") +
  theme_classic(base_size = 16)

# figure with polydora blister count vs dermo intensity
fig1_right <- ggplot(means_by_dermo, aes(x = Score_Rectum, y = mean_blister)) +
  geom_line(color = "black", linewidth = 1.2) +
  geom_point(size = 4, color = "black") +
  geom_errorbar(aes(ymin = mean_blister - se_blister, ymax = mean_blister + se_blister),
                width = 0.2, color = "black", linewidth = 0.8) +
  geom_text(data = cld_blister,
            aes(x = as.numeric(DermoIntensity), y = mean_blister + se_blister + 1, label = letter),
            size = 5, fontface = "bold") +
  scale_x_continuous(breaks = c(0.5, 1, 2, 3, 4, 5)) +
  labs(x = "Dermo Intensity", y = "Polydora Blister Count") +
  theme_classic(base_size = 16)

# combined panel
fig1 <- fig1_left + fig1_right
print(fig1)

# figure with odds ratios from best model
coef_data <- data.frame(
  Variable = c("Polydora Shell Cover", "Year 2024", "Year 2025", "River South"),
  OR = exp(c(
    coef(model_M2)["Polydora_outside"],
    coef(model_M2)["Year2024"],
    coef(model_M2)["Year2025"],
    coef(model_M2)["RiverSouth"])),
  Lower = exp(c(
    coef(model_M2)["Polydora_outside"] - 1.96 * summary(model_M2)$coefficients["Polydora_outside", "Std. Error"],
    coef(model_M2)["Year2024"] - 1.96 * summary(model_M2)$coefficients["Year2024", "Std. Error"],
    coef(model_M2)["Year2025"] - 1.96 * summary(model_M2)$coefficients["Year2025", "Std. Error"],
    coef(model_M2)["RiverSouth"] - 1.96 * summary(model_M2)$coefficients["RiverSouth", "Std. Error"])),
  Upper = exp(c(
    coef(model_M2)["Polydora_outside"] + 1.96 * summary(model_M2)$coefficients["Polydora_outside", "Std. Error"],
    coef(model_M2)["Year2024"] + 1.96 * summary(model_M2)$coefficients["Year2024", "Std. Error"],
    coef(model_M2)["Year2025"] + 1.96 * summary(model_M2)$coefficients["Year2025", "Std. Error"],
    coef(model_M2)["RiverSouth"] + 1.96 * summary(model_M2)$coefficients["RiverSouth", "Std. Error"])),
  P_value = c(
    summary(model_M2)$coefficients["Polydora_outside", "Pr(>|z|)"],
    summary(model_M2)$coefficients["Year2024", "Pr(>|z|)"],
    summary(model_M2)$coefficients["Year2025", "Pr(>|z|)"],
    summary(model_M2)$coefficients["RiverSouth", "Pr(>|z|)"]))

coef_data$Sig <- ifelse(coef_data$P_value < 0.001, "***",
                        ifelse(coef_data$P_value < 0.01, "**",
                               ifelse(coef_data$P_value < 0.05, "*", "ns")))

fig2 <- ggplot(coef_data, aes(x = reorder(Variable, OR), y = OR)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, linewidth = 1) +
  geom_text(aes(label = Sig), vjust = -1.5, size = 6) +
  coord_flip() +
  labs(x = "", y = "Odds Ratio (95% CI)") +
  theme_classic(base_size = 16) +
  theme(axis.text.y = element_text(size = 14))

print(fig2)

# fig with temporal trends in dermo intensity
mean_dermo_year <- coinfected_data %>%
  group_by(Year, River) %>%
  summarise(
    mean_dermo = mean(Score_Rectum, na.rm = TRUE),
    se_dermo = sd(Score_Rectum, na.rm = TRUE) / sqrt(n()),
    .groups = "drop")

fig3 <- ggplot(mean_dermo_year, aes(x = Year, y = mean_dermo, color = River, group = River)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = mean_dermo - se_dermo, ymax = mean_dermo + se_dermo),
                width = 0.15, linewidth = 1) +
  scale_color_manual(values = river_colors) +
  annotate("text", x = 2.8, y = 3,
           label = expression(paste("2024: ", beta, " = -1.40, p < 0.001")),
           size = 5, fontface = "bold") +
  labs(x = "Year", y = "Mean Dermo Intensity") +
  theme_classic(base_size = 16) +
  theme(legend.position = c(0.85, 0.85),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white", color = "black"))

print(fig3)
