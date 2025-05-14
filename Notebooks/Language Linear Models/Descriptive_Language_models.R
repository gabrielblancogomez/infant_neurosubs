# Descriptive statistics of language 

# ============================================================================ #
# 1. SETUP AND DATA PREPARATION
# ============================================================================ #
# Load required packages
library(tidyverse)       # For data manipulation and visualization
library(lme4)            # For mixed-effects models
library(lmerTest)        # For p-values in mixed models
library(performance)     # For model diagnostics
library(emmeans)         # For estimated marginal means
library(viridis)         # For color-blind friendly palettes
library(psych)           # For descriptive statistics
library(car)             # For additional diagnostics
library(naniar)          # For missing data visualization
library(gt)              # For publication-ready tables
library(gtsummary)       # For summary tables
library(flextable)       # Alternative for tables
library(ggridges)        # For density ridge plots
library(cowplot)         # For arranging multiple plots
library(MuMIn)  # For r.squaredGLMM()
library(mice) # For imputation 
library(effsize) # For effect size calculations


# Load the data
lang_data <- read.csv("C:/Users/gabot/OneDrive - McGill University/Desktop/Github_repos/infant_neurosubs/Datasets/neurosub_solutions/neurosubs_df_2025_04_24_5D.csv", header = TRUE, sep = ",")


# Check data structure
print("Data Structure Summary:")
print(paste("Number of participants:", nrow(lang_data)))
print(paste("Variables in dataset:", ncol(lang_data)))


# Get column names for receptive language measures at different timepoints
receptive_cols <- grep("receptive_", names(lang_data), value = TRUE)

# Check for missing data in receptive language measures
missing_summary <- sapply(lang_data[receptive_cols], function(x) sum(is.na(x)))
missing_percent <- round(missing_summary / nrow(lang_data) * 100, 1)
missing_df <- data.frame(
  Timepoint = receptive_cols,
  Missing_Count = missing_summary,
  Missing_Percent = missing_percent
)
print("Missing data by timepoint:")
print(missing_df)

# Check normality of receptive language measures
print("Normality check by timepoint (Shapiro-Wilk test):")
norm_tests <- sapply(lang_data[receptive_cols], function(x) {
  if(sum(!is.na(x)) > 3) {
    shapiro.test(x[!is.na(x)])$p.value
  } else {
    NA
  }
})
norm_df <- data.frame(
  Timepoint = receptive_cols,
  p_value = round(norm_tests, 3),
  Normal = ifelse(norm_tests > 0.05, "Yes", "No")
)
print(norm_df)

# Check for missing values in nonverbal IQ
print("Missing values in nonverbal IQ at 6 months:")
print(sum(is.na(lang_data$nonverbal_iq_6)))
print(paste("Percent missing:", round(sum(is.na(lang_data$nonverbal_iq_6)) / nrow(lang_data) * 100, 1), "%"))

# ============================================================================ #
# 2. ENHANCED MISSINGNESS VISUALIZATION
# ============================================================================ #

# Create a missingness pattern visualization across all timepoints
missing_vis <- lang_data %>%
  select(subject, lpa_class, receptive_cols) %>%
  pivot_longer(
    cols = all_of(receptive_cols),
    names_to = "timepoint",
    values_to = "value"
  ) %>%
  mutate(
    timepoint = gsub("receptive_", "", timepoint),
    is_missing = is.na(value)
  )

# Plot 1: Overall missingness patterns
plot_missing_overall <- ggplot(missing_vis, aes(x = timepoint, fill = is_missing)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("FALSE" = "#3498db", "TRUE" = "#e74c3c"),
                    labels = c("FALSE" = "Present", "TRUE" = "Missing")) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    title = "Data Missingness by Timepoint",
    subtitle = "Proportion of missing values at each assessment",
    x = "Age (months)",
    y = "Proportion",
    fill = "Data Status"
  )

# Plot 2: Missingness by LPA class
plot_missing_by_class <- ggplot(missing_vis, aes(x = timepoint, fill = is_missing)) +
  geom_bar(position = "fill") +
  facet_wrap(~lpa_class, ncol = 1) +
  scale_fill_manual(values = c("FALSE" = "#3498db", "TRUE" = "#e74c3c"),
                    labels = c("FALSE" = "Present", "TRUE" = "Missing")) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    title = "Data Missingness by LPA Class and Timepoint",
    x = "Age (months)",
    y = "Proportion",
    fill = "Data Status"
  )

# Plot 3: Subject-level missingness pattern (UpSet plot approach)
# Create a binary matrix of missingness (1=missing, 0=present)
miss_matrix <- lang_data %>%
  select(subject, receptive_cols) %>%
  mutate(across(all_of(receptive_cols), ~ifelse(is.na(.), 1, 0))) %>%
  select(-subject)

# Count combinations of missingness
miss_patterns <- miss_matrix %>%
  group_by_all() %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(desc(count))

# Convert to long form for visualization
miss_patterns_long <- miss_patterns %>%
  mutate(pattern_id = row_number()) %>%
  pivot_longer(
    cols = all_of(receptive_cols),
    names_to = "timepoint",
    values_to = "is_missing"
  ) %>%
  mutate(
    timepoint = gsub("receptive_", "", timepoint)
  )

# Create missingness combinations plot
plot_miss_patterns <- ggplot(miss_patterns_long, 
                             aes(x = timepoint, y = factor(pattern_id), fill = factor(is_missing))) +
  geom_tile(color = "white") +
  geom_text(aes(label = count), data = miss_patterns %>% mutate(pattern_id = row_number()), 
            x = "receptive_6mo", hjust = -0.5, size = 3) +
  scale_fill_manual(values = c("0" = "#3498db", "1" = "#e74c3c"),
                    labels = c("0" = "Present", "1" = "Missing")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  labs(
    title = "Common Patterns of Missing Data",
    subtitle = "Each row represents a distinct pattern of missing values",
    x = "Timepoint (months)",
    y = "Pattern",
    fill = "Data Status"
  )

# Plot 4: Individual subject missingness trajectories (subsample for clarity)
set.seed(123)  # For reproducibility
sampled_subjects <- sample(unique(lang_data$subject), 
                           min(30, length(unique(lang_data$subject))))

# Create trajectory plot for sampled subjects
subject_miss <- missing_vis %>%
  filter(subject %in% sampled_subjects)

plot_subject_miss <- ggplot(subject_miss, 
                            aes(x = timepoint, y = subject, fill = is_missing)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("FALSE" = "#3498db", "TRUE" = "#e74c3c"),
                    labels = c("FALSE" = "Present", "TRUE" = "Missing")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Individual Subject Data Availability",
    subtitle = "Sample of 30 subjects showing data presence/absence patterns",
    x = "Timepoint (months)",
    y = "Individual Subjects",
    fill = "Data Status"
  )
plot_missing_overall
plot_missing_by_class
 

# Combine missingness plots into a single figure
missing_plots <- plot_grid(
  plot_missing_overall,
  plot_missing_by_class,
  ncol = 2,
  labels = "AUTO"
)

missing_plots
# Save the combined plot


ggsave("./missingness_patterns.pdf", missing_plots, width = 12, height = 10)

# =========================================================================== #
# DESCRIPTIVE ANALYSIS OF LANGUAGE VARS
# =========================================================================== #  


to_long <- function(data, cols, value_name) {
  long <- data %>%
    pivot_longer(
      cols = all_of(cols),
      names_to = "timepoint_var",
      values_to = value_name
    ) %>%
    mutate(
      timepoint = as.numeric(gsub(".*_", "", gsub("mo", "", timepoint_var))),
      timepoint_centered = timepoint - 6,
      timepoint_factor = factor(timepoint),
      hc_class = factor(hc_class),
      site = factor(site),
      sex = factor(sex),
      group_type = factor(group_type),
      subject = factor(subject)
    )
  return(long)
}

# ============================================================================
# EXPRESSIVE LANGUAGE ANALYSIS
# ============================================================================
long_expressive <- to_long(lang_data, expressive_cols, "expressive")

desc_expressive <- long_expressive %>%
  group_by(timepoint) %>%
  summarise(
    n = sum(!is.na(expressive)),
    mean = mean(expressive, na.rm = TRUE),
    sd = sd(expressive, na.rm = TRUE),
    min = min(expressive, na.rm = TRUE),
    max = max(expressive, na.rm = TRUE),
    variance = var(expressive, na.rm = TRUE),
    .groups = "drop"
  )

# HTML Table (gt)
gt_expressive <- desc_expressive %>%
  gt() %>%
  tab_header(
    title = "Descriptive Statistics of Expressive Language",
    subtitle = "Overall by Timepoint (Months)"
  ) %>%
  fmt_number(columns = c(mean, sd, variance), decimals = 2) %>%
  cols_label(
    timepoint = "Timepoint (months)",
    n = "N", mean = "Mean", sd = "SD",
    min = "Min", max = "Max", variance = "Variance"
  )


gtsave(gt_expressive, "C:/Users/gabot/OneDrive - McGill University/Desktop/Github_repos/infant_neurosubs/Tables/Supplementary/expressive_language_descrptive.html")

# Word Table
ft_expressive <- flextable(desc_expressive)
docx_expressive <- read_docx() %>%
  body_add_par("Descriptive Statistics of Expressive Language", style = "heading 1") %>%
  body_add_flextable(ft_expressive)
print(docx_expressive, target = "C:/Users/gabot/OneDrive - McGill University/Desktop/Github_repos/infant_neurosubs/Tables/Supplementary/expressive_language_descrptive.docx")

# ============================================================================
# RECEPTIVE LANGUAGE ANALYSIS
# ============================================================================
long_receptive <- to_long(lang_data, receptive_cols, "receptive")

desc_receptive <- long_receptive %>%
  group_by(timepoint) %>%
  summarise(
    n = sum(!is.na(receptive)),
    mean = mean(receptive, na.rm = TRUE),
    sd = sd(receptive, na.rm = TRUE),
    min = min(receptive, na.rm = TRUE),
    max = max(receptive, na.rm = TRUE),
    variance = var(receptive, na.rm = TRUE),
    .groups = "drop"
  )

# HTML Table (gt)
gt_receptive <- desc_receptive %>%
  gt() %>%
  tab_header(
    title = "Descriptive Statistics of Receptive Language",
    subtitle = "Overall by Timepoint (Months)"
  ) %>%
  fmt_number(columns = c(mean, sd, variance), decimals = 2) %>%
  cols_label(
    timepoint = "Timepoint (months)",
    n = "N", mean = "Mean", sd = "SD",
    min = "Min", max = "Max", variance = "Variance"
  )

gtsave(gt_receptive,"C:/Users/gabot/OneDrive - McGill University/Desktop/Github_repos/infant_neurosubs/Tables/Supplementary/receptive_language_descrptive.html")

# Word Table
ft_receptive <- flextable(desc_receptive)
docx_receptive <- read_docx() %>%
  body_add_par("Descriptive Statistics of Receptive Language", style = "heading 1") %>%
  body_add_flextable(ft_receptive)
print(docx_receptive, target = "C:/Users/gabot/OneDrive - McGill University/Desktop/Github_repos/infant_neurosubs/Tables/Supplementary/receptive_language_descrptive.docx")


# ============================================================================ #

# SEX DIFFERNCES

# ============================================================================


# Asses if there are differences in sex between the different EEG variables

eeg_variables <- c( "front_gamma_6","auditory_con_6" ,"lang_comp_con_6","speech_con_6_left", 
                   "gamma_lat_6")

# Conduct students t-test or Mann Whitney U test depending on normality

sex_results_df <- data.frame(
  Variable = character(),
  n_female=numeric(),
  n_male=numeric(),
  Test = character(),
  Statistic = numeric(),
  DF = numeric(),
  P_value = numeric(),
  Effect_Size = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each EEG variable
for (var in eeg_variables) {
  
  # Check for missing data
  var_data <- lang_data[[var]]
  sex <- lang_data$sex
  valid_idx <- !is.na(var_data) & !is.na(sex)
  var_data <- var_data[valid_idx]
  sex <- sex[valid_idx]
  # get n of males
  nmales= lang_data %>% filter( sex=="M") %>% nrow()
  nfemales= lang_data %>% filter( sex=="F") %>% nrow()
  
  # Skip if not enough data
  if (length(unique(sex)) < 2 || length(var_data) < 3) next
  
  # Shapiro-Wilk test for normality
  shapiro_p <- shapiro.test(var_data)$p.value
  is_normal <- shapiro_p > 0.05
  
  if (is_normal) {
    test_name <- "Independent t-test"
    test_result <- t.test(var_data ~ sex)
    effect_size <- effsize::cohen.d(var_data ~ sex)$estimate
    df_val <- as.numeric(test_result$parameter)
    stat_val <- as.numeric(test_result$statistic)
    p_val <- test_result$p.value
    n_female=nfemales
    n_male=nmales
  } else {
    test_name <- "Wilcoxon rank-sum test"
    test_result <- wilcox.test(var_data ~ sex)
    effect_size <- effsize::cliff.delta(var_data ~ sex)$estimate
    df_val <- NA
    stat_val <- as.numeric(test_result$statistic)
    p_val <- test_result$p.value
    n_female=nfemales
    n_male=nmales
  }
  
  # Append results to dataframe
  sex_results_df <- rbind(
    sex_results_df,
    data.frame(
      Variable = var,
      Test = test_name,
      n_female=nfemales,
      n_male=nmales,
      Statistic = round(stat_val, 3),
      DF = ifelse(is.na(df_val), NA, round(df_val, 1)),
      P_value = signif(p_val, 3),
      Effect_Size = round(effect_size, 3),
      stringsAsFactors = FALSE

    )
  )
}

# Create a flextable for Word
sex_results_ft <- flextable(sex_results_df)
sex_results_ft <- autofit(sex_results_ft)
sex_results_ft
 
# Save as Word document
doc <- read_docx() %>%
  body_add_par("EEG Variable Comparisons by Sex", style = "heading 1") %>%
  body_add_flextable(sex_results_ft)

print(doc, target = "C:/Users/gabot/OneDrive - McGill University/Desktop/Github_repos/infant_neurosubs/Tables/Supplementary/eeg_comparisons_by_sex.docx")


# ============================================================================ #
# 4. VISUALIZATION OF RAW DATA AND TRAJECTORIES
# ============================================================================ #

# Plot 1: Raw data with group means
p1 <- ggplot(long_data, aes(x = timepoint, y = receptive, color = lpa_class)) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun = mean, geom = "line", aes(group = lpa_class), size = 1.2) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.5) +
  scale_color_viridis_d(end = 0.8, option = "D") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    text = element_text(size = 12)
  ) +
  labs(
    title = "Language Receptive Trajectories by LPA Class",
    x = "Age (months)",
    y = "Receptive Language Score",
    color = "LPA Class"
  )

# Plot 2: Individual trajectories with smoothed group trends
p2 <- ggplot(long_data, aes(x = timepoint, y = receptive, group = subject, color = lpa_class)) +
  geom_line(alpha = 0.15) +
  geom_smooth(aes(group = lpa_class), method = "loess", se = TRUE, span = 1) +
  scale_color_viridis_d(end = 0.8, option = "D") +
  facet_wrap(~lpa_class) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    text = element_text(size = 12)
  ) +
  labs(
    title = "Individual Trajectories with LPA Class Trends",
    x = "Age (months)",
    y = "Receptive Language Score"
  )

# Plot 3: Density ridges by timepoint and LPA class
p3 <- ggplot(long_data, aes(x = receptive, y = timepoint_factor, fill = lpa_class)) +
  geom_density_ridges(alpha = 0.7, scale = 0.9, rel_min_height = 0.01) +
  scale_fill_viridis_d(end = 0.8, option = "D") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    text = element_text(size = 12)
  ) +
  labs(
    title = "Distribution of Receptive Language by Age and LPA Class",
    x = "Receptive Language Score",
    y = "Age (months)",
    fill = "LPA Class"
  )
p1
p2
p3
# Combine plots
trajectory_plots <- plot_grid(
  p1, p3,
  ncol = 1,
  labels = "AUTO"
)
trajectory_plots
# Save the combined plot
ggsave("trajectory_visualizations.pdf", trajectory_plots, width = 10, height = 15)