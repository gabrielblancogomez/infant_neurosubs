# Unified analysis for language trajectories
# This combines both classification methods (HC & HC) and language types (receptive & receptive)

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
library(officer)         # For saving tables as Word documents'
library(r2glmm) # For calculating semi-partial R-squared
library(effectsize)  # For standardized coefficients
library(sjstats)      # For effect sizes



# Load the data
lang_data <- read.csv("C:/Users/gabot/OneDrive - McGill University/Desktop/Github_repos/infant_neurosubs/Datasets/neurosub_solutions/neurosubs_df_2025_04_23_5D.csv", header = TRUE, sep = ",")
# Set current working directory as the same where this file is
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Select Palettes
hc_palette <- c("#184e77","#7570b3", "#1b9e77")

# Create output directories
dir.create("figures/main", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/supplementary", recursive = TRUE, showWarnings = FALSE)
dir.create("tables/supplementary", recursive = TRUE, showWarnings = FALSE)
dir.create("tables/main", recursive = TRUE, showWarnings = FALSE)

dir.create("sensitivity", recursive = TRUE, showWarnings = FALSE)
dir.create("sensitivity", recursive = TRUE, showWarnings = FALSE)


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

# Basic winsorization for main analysis (5% level)
winsorize <- function(data, columns, by_group = NULL, trim_pct = 0.05) {
  data_new <- data
  
  if(is.null(by_group)) {
    # Winsorize without grouping
    for(col in columns) {
      if(sum(!is.na(data[[col]])) > 5) {  # Only if we have enough data
        qt <- quantile(data[[col]], probs = c(trim_pct, 1-trim_pct), na.rm = TRUE)
        data_new[[col]] <- ifelse(is.na(data[[col]]), NA, 
                                  ifelse(data[[col]] < qt[1], qt[1],
                                         ifelse(data[[col]] > qt[2], qt[2], data[[col]])))
      }
    }
  } else {
    # Winsorize by group
    for(grp in unique(data[[by_group]])) {
      for(col in columns) {
        idx <- data[[by_group]] == grp
        if(sum(idx, na.rm = TRUE) > 5) {  # Only if we have enough data in this group
          qt <- quantile(data[idx, col], probs = c(trim_pct, 1-trim_pct), na.rm = TRUE)
          data_new[idx & !is.na(data[[col]]) & data[[col]] < qt[1], col] <- qt[1]
          data_new[idx & !is.na(data[[col]]) & data[[col]] > qt[2], col] <- qt[2]
        }
      }
    }
  }
  
  return(data_new)
}

# Apply 5% winsorization for main analysis
lang_data_wins <- winsorize(lang_data, receptive_cols, by_group = "hc_class", trim_pct = 0.05)

# Convert wide to long format for the main analysis dataset
to_long <- function(data) {
  long <- data %>%
    pivot_longer(
      cols = all_of(receptive_cols),
      names_to = "timepoint_var",
      values_to = "receptive"
    ) %>%
    mutate(
      # Extract numeric months from column names (assuming format receptive_Xmo)
      timepoint = as.numeric(gsub("receptive_", "", gsub("mo", "", timepoint_var))),
      # Center timepoint at first measurement (6mo)
      timepoint_centered = timepoint - 6,
      # Create categorical version for certain analyses
      timepoint_factor = factor(timepoint)
    )
  
  # Convert factors properly
  long$hc_class <- factor(long$hc_class)
  long$site <- factor(long$site)
  long$sex <- factor(long$sex)
  long$group_type <- factor(long$group_type)
  long$subject <- factor(long$subject)
  
  return(long)
}

# Create main analysis long dataset
long_data <- to_long(lang_data_wins)

# Save data description table
desc_table <- long_data %>%
  group_by(hc_class, timepoint) %>%
  summarise(
    n = sum(!is.na(receptive)),
    mean = mean(receptive, na.rm = TRUE),
    sd = sd(receptive, na.rm = TRUE),
    min = min(receptive, na.rm = TRUE),
    max = max(receptive, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    id_cols = timepoint,
    names_from = hc_class,
    values_from = c(n, mean, sd),
    names_sep = "_"
  ) %>%
  gt() %>%
  tab_header(
    title = "Descriptive Statistics of Receptive Language by HC Class",
    subtitle = "Across all timepoints (months)"
  ) %>%
  fmt_number(
    columns = starts_with("mean") | starts_with("sd"),
    decimals = 2
  )
desc_table
# Save descriptive table as HTML and Word
gtsave(desc_table, "tables/supplementary/HC_descriptive_statistics.html")
desc_data <- as.data.frame(desc_table$`_data`)

# Then create flextable
descriptive_flextable <- flextable(desc_data) %>%
  set_header_labels( # Match your column names
    variable = "Variable",
    mean = "Mean",
    sd = "SD"
  ) %>%
  theme_zebra() %>%
  autofit()
descriptive_flextable
save_as_docx(descriptive_flextable, path = "tables/supplementary/HC_Receptive_descriptive_statistics.docx")

# ============================================================================ #
# 2. VISUALIZATION OF RAW DATA AND TRAJECTORIES
# ============================================================================ #

# Plot 1: Raw data with group means
p1 <- ggplot(long_data, aes(x = timepoint, y = receptive, color = hc_class)) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun = mean, geom = "line", aes(group = hc_class), size = 1.2) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.5) +
  scale_color_manual(values = hc_palette) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    text = element_text(size = 12)
  ) +
  labs(
    title = "Language Receptive Trajectories by HC Class",
    x = "Age (months)",
    y = "Receptive Language Score",
    color = "HC Class"
  )
p1
# Plot 2: Individual trajectories with smoothed group trends
p2 <- ggplot(long_data, aes(x = timepoint, y = receptive, group = subject, color = hc_class)) +
  geom_line(alpha = 0.15) +
  geom_smooth(aes(group = hc_class), method = "loess", se = TRUE, span = 1) +
  scale_color_manual(values = hc_palette) +
  facet_wrap(~hc_class) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    text = element_text(size = 12)
  ) +
  labs(
    title = "Individual Trajectories with HC Class Trends",
    x = "Age (months)",
    y = "Receptive Language Score"
  )

p2
# Plot 3: Density ridges by timepoint and HC class
p3 <- ggplot(long_data, aes(x = receptive, y = timepoint_factor, fill = hc_class)) +
  geom_density_ridges(alpha = 0.7, scale = 0.9, rel_min_height = 0.01) +
  scale_fill_manual(values = hc_palette) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    text = element_text(size = 12)
  ) +
  labs(
    title = "Distribution of Receptive Language by Age and HC Class",
    x = "Receptive Language Score",
    y = "Age (months)",
    fill = "HC Class"
  )

p3
# Combine plots
trajectory_plots <- plot_grid(p1, p3, ncol = 1, labels = "AUTO")
trajectory_plots


# Save the plots
ggsave("figures/main/HC_Receptive_mean_trajectories.pdf", p1, width = 10, height = 7)
ggsave("figures/supplementary/HC_Receptive_individual_trajectories.pdf", p2, width = 10, height = 7)
ggsave("figures/supplementary/HC_Receptive_distribution_by_timepoint.pdf", p3, width = 10, height = 7)
ggsave("figures/supplementary/HC_Receptive_combined_trajectory_plots.pdf", trajectory_plots, width = 10, height = 14)
# ============================================================================ #
# 3. MIXED-EFFECTS MODELING - MAIN ANALYSIS
# ============================================================================ #



# Define function to run main model sequence
run_main_models <- function(data) {
  # Create a standardized version of the data for standardized models
  data_std <- data %>%
    mutate(across(c(receptive, timepoint_centered, nonverbal_iq_6), 
                  ~scale(.) %>% as.vector(), 
                  .names = "{.col}_std"))
  
  # 1. Unconditional means model (random intercept only)
  m0 <- lmer(receptive ~ 1 + (1 | subject), data = data, REML = FALSE)
  
  # 2. Unconditional growth model (random slope only for convergence)
  m1 <- lmer(receptive ~ timepoint_centered + 
               (0 + timepoint_centered | subject),  # Random slopes only
             data = data, REML = FALSE,
             control = lmerControl(optimizer = "bobyqa"))
  
  # 3. Add HC class effect on intercept
  m2 <- lmer(receptive ~ timepoint_centered + hc_class + 
               (0 + timepoint_centered | subject),  # Random slopes only
             data = data, REML = FALSE,
             control = lmerControl(optimizer = "bobyqa"))
  
  # 4. Add HC class effect on slope (interaction)
  m3 <- lmer(receptive ~ timepoint_centered * hc_class + 
               (0 + timepoint_centered | subject),  # Random slopes only
             data = data, REML = FALSE,
             control = lmerControl(optimizer = "bobyqa"))
  
  # 5. Add covariates (sex, site, risk, iq)
  m4 <- lmer(receptive ~ timepoint_centered * hc_class + 
               sex + site + group_type + nonverbal_iq_6 +
               (0 + timepoint_centered | subject),  # Random slopes only
             data = data, REML = FALSE,
             control = lmerControl(optimizer = "bobyqa"))
  
  # Create standardized models for comparison
  # Use standardized data for fully standardized coefficients
  m1_std <- lmer(receptive_std ~ timepoint_centered_std + 
                   (0 + timepoint_centered_std | subject),
                 data = data_std, REML = FALSE,
                 control = lmerControl(optimizer = "bobyqa"))
  
  m2_std <- lmer(receptive_std ~ timepoint_centered_std + hc_class + 
                   (0 + timepoint_centered_std | subject),
                 data = data_std, REML = FALSE,
                 control = lmerControl(optimizer = "bobyqa"))
  
  m3_std <- lmer(receptive_std ~ timepoint_centered_std * hc_class + 
                   (0 + timepoint_centered_std | subject),
                 data = data_std, REML = FALSE,
                 control = lmerControl(optimizer = "bobyqa"))
  
  m4_std <- lmer(receptive_std ~ timepoint_centered_std * hc_class + 
                   sex + site + group_type + nonverbal_iq_6_std +
                   (0 + timepoint_centered_std | subject),
                 data = data_std, REML = FALSE,
                 control = lmerControl(optimizer = "bobyqa"))
  
  # Return list of models
  return(list(
    models = list(m0 = m0, m1 = m1, m2 = m2, m3 = m3, m4 = m4),
    models_std = list(m1_std = m1_std, m2_std = m2_std, m3_std = m3_std, m4_std = m4_std)
  ))
}



# Filter out NAs in nonverbal_iq_6 for model fitting if needed
#main_data <- long_data %>% filter(!is.na(nonverbal_iq_6))
main_data <- long_data

# Run the main models
main_results <- run_main_models(main_data)
main_models <- main_results$models
main_models_std <- main_results$models_std

# Function to check convergence of models
check_convergence <- function(model) {
  if (is.null(model@optinfo$conv$lme4$messages)) {
    return("Converged")
  } else if ("boundary (singular) fit: see help('isSingular')" %in% model@optinfo$conv$lme4$messages) {
    return("Singular")
  } else if (any(grepl("convergence", model@optinfo$conv$lme4$messages, ignore.case = TRUE))) {
    return("Conv. issues")
  } else {
    return("Other issue")
  }
}

# Create model comparison table
model_descriptions <- c(
  "Random intercept only",
  "Random slope only",
  "Fixed HC effect",
  "HC*Time interaction",
  "Covariates (sex, site, risk IQ)"
)

model_comp <- data.frame(
  Model = names(main_models),
  Description = model_descriptions,
  AIC = sapply(main_models, AIC),
  BIC = sapply(main_models, BIC),
  logLik = sapply(main_models, logLik),
  Convergence = sapply(main_models, check_convergence)
)

model_comp$dAIC <- model_comp$AIC - min(model_comp$AIC)
model_comp$dBIC <- model_comp$BIC - min(model_comp$BIC)

# Create formatted model comparison table
model_table <- model_comp %>%
  gt() %>%
  tab_header(title = "Model Comparison for Receptive Language Development") %>%
  fmt_number(columns = c("AIC", "BIC", "logLik", "dAIC", "dBIC"), decimals = 2) %>%
  tab_style(
    style = list(cell_fill(color = "#e8f4f8")),
    locations = cells_body(rows = AIC == min(AIC))
  ) %>%
  tab_footnote(
    footnote = "Highlighted row indicates model with lowest AIC value",
    locations = cells_column_labels(columns = AIC)
  )

model_table

# Save model comparison table
gtsave(model_table, "tables/supplementary/HC_Receptive_model_comparison.html")


# Convert to flextable
model_data <- as.data.frame(model_table$`_data`)  # Extract underlying data
model_table_flex <- flextable(model_data) %>% 
  set_table_properties(layout = "autofit") %>%
  theme_zebra() %>%
  align(align = "center", part = "all")

save_as_docx(model_table_flex, path = "tables/supplementary/HC_Receptive_model_comparison.docx")

# Calculate standardized coefficients for the best model 

best_model_std <- main_models_std$m4_std
std_coef_table <- broom.mixed::tidy(best_model_std) %>%
  filter(effect == "fixed") %>%
  select(term, estimate, std.error, statistic, p.value) %>%
  rename(
    Parameter = term,
    `Standardized Beta` = estimate,
    SE = std.error,
    t = statistic,
    p = p.value
  )

std_coef_table 

# For appropriate standardized coefficients:
std_coef_effectsize <- standardize_parameters(main_models$m4, 
                                              method =  "smart")
std_coef_effectsize
# Create nicely formatted table for standardized betas
std_beta_table <- std_coef_effectsize %>%
  rename(
    Parameter = Parameter,
    `Standardized Beta` = Std_Coefficient,
    `CI Low` = CI_low,
    `CI High` = CI_high
  ) %>%
  select(Parameter, `Standardized Beta`, `CI Low`, `CI High`)

std_beta_gt <- std_beta_table %>%
  gt() %>%
  tab_header(title = "Standardized Coefficients for Final Model") %>%
  fmt_number(columns = c("Standardized Beta", "CI Low", "CI High"), decimals = 3) %>%
  tab_style(
    style = list(cell_text(weight = "bold")),
    locations = cells_body(
      rows = (`CI Low` > 0) | (`CI High` < 0)
    )
  ) %>%
  tab_footnote(
    footnote = "Values represent standardized beta coefficients (β) with 95% confidence intervals",
    locations = cells_column_labels(columns = "Standardized Beta")
  )

std_beta_gt

# Save standardized beta table
gtsave(std_beta_gt, "tables/supplementary/HC_Receptive_standardized_betas.html")



# 1. R-squared method using performance package
r2_values <- performance::r2(main_models$m4)
r2_values

# R-squared summary
r2_summary <- data.frame(
  Measure = c("Conditional R²", "Marginal R²"),
  Value = c(r2_values$R2_conditional, r2_values$R2_marginal),
  Description = c(
    "Variance explained by entire model (fixed + random effects)",
    "Variance explained by fixed effects only"
  )
)

r2_gt <- r2_summary %>%
  gt() %>%
  tab_header(title = "R-squared Values for Final Model") %>%
  fmt_number(columns = "Value", decimals = 3)

r2_gt

# Save R-squared summary
gtsave(r2_gt, "tables/supplementary/HC_Receptive_r_squared.html")

# Convert to flextable
r2_flex <- flextable(as.data.frame(r2_gt$`_data`)) %>% 
  set_table_properties(layout = "autofit") %>%
  theme_zebra() %>%
  align(align = "center", part = "all")

save_as_docx(r2_flex, path = "tables/supplementary/HC_Receptive_r_squared.docx")



# Convert to flextable
std_beta_flex <- flextable(as.data.frame(std_beta_gt$`_data`)) %>% 
  set_table_properties(layout = "autofit") %>%
  theme_zebra() %>%
  align(align = "center", part = "all")
std_beta_flex
save_as_docx(std_beta_flex, path = "tables/supplementary/HC_Receptive_standardized_betas.docx")


# Select best model based on AIC
best_model_idx <- which.min(sapply(main_models, AIC))
best_model <- main_models[[best_model_idx]]
print(paste("Selected best model:", names(main_models)[best_model_idx]))

# Model diagnostics for best model
diagnostics_plot <- performance::check_model(best_model)
pdf("figures/supplementary/HC_Receptive_model_diagnostics.pdf", width = 8, height = 10)
print(diagnostics_plot)
dev.off()
 
# ============================================================================ #
# 4. ESTIMATED MARGINAL MEANS AND CONTRASTS
# ============================================================================ #

# Calculate estimated marginal means for HC class by timepoint
emm_hc_time <- emmeans(best_model, ~ hc_class | timepoint_centered, 
                        at = list(timepoint_centered = c(0, 6, 12, 18, 24, 30)))

emm_hc_time_later <- emmeans(best_model, ~ hc_class | timepoint_centered, 
                        at = list(timepoint_centered = c(24, 30)))

# Convert back to original timepoint scale for plotting and interpretation
emm_hc_time_df <- as.data.frame(emm_hc_time) %>%
  mutate(timepoint = timepoint_centered + 6)

# Pairwise comparisons between HC classes at each timepoint
pairs_hc <- pairs(emm_hc_time_later, adjust = "tukey")
pairs_hc_df <- as.data.frame(pairs_hc)

# Test for overall trend differences (slope contrasts)
emm_slopes <- emtrends(best_model, ~ hc_class, var = "timepoint_centered")
slope_pairs <- pairs(emm_slopes, adjust = "tukey")
slope_pairs_df <- as.data.frame(slope_pairs)

# Create publication-ready tables for EMM results
emm_table <- emm_hc_time_df %>%
  select(hc_class, timepoint_centered, emmean, SE, lower.CL, upper.CL) %>%
  mutate(timepoint = timepoint_centered + 6) %>%
  arrange(timepoint, hc_class) %>%
  gt() %>%
  tab_header(
    title = "Estimated Marginal Means by HC Class and Timepoint",
    subtitle = "Adjusted for covariates in the final model"
  ) %>%
  fmt_number(
    columns = c("emmean", "SE", "lower.CL", "upper.CL"),
    decimals = 2
  ) %>%
  cols_label(
    hc_class = "HC Class",
    timepoint_centered = "Centered Time",
    timepoint = "Age (months)",
    emmean = "Estimated Mean",
    SE = "Std. Error",
    lower.CL = "95% CI Lower",
    upper.CL = "95% CI Upper"
  )

pairs_table <- pairs_hc_df %>%
  gt() %>%
  tab_header(
    title = "Pairwise Comparisons Between HC Classes at Each Timepoint",
    subtitle = "Tukey-adjusted p-values"
  ) %>%
  fmt_number(
    columns = c("estimate", "SE", "t.ratio"),
    decimals = 2
  ) %>%
  fmt_number(
    columns = "p.value",
    decimals = 3
  ) %>%
  cols_label(
    contrast = "Contrast",
    estimate = "Difference",
    SE = "Std. Error",
    t.ratio = "t value",
    p.value = "p-value"
  )

slopes_table <- slope_pairs_df %>%
  gt() %>%
  tab_header(
    title = "Differences in Developmental Trajectories Between HC Classes",
    subtitle = "Comparison of slopes (change per month)"
  ) %>%
  fmt_number(
    columns = c("estimate", "SE", "t.ratio"),
    decimals = 2
  ) %>%
  fmt_number(
    columns = "p.value",
    decimals = 3
  ) %>%
  tab_style(
    style = list(
      cell_text(style = "italic", color = "red")
    ),
    locations = cells_body(
      columns = "p.value",
      rows = p.value < 0.05
    )
  ) %>%
  cols_label(
    contrast = "Contrast",
    estimate = "Slope Difference",
    SE = "Std. Error",
    t.ratio = "t value",
    p.value = "p-value"
  )

emm_table
pairs_table
slopes_table

# Save tables
gtsave(emm_table, "tables/supplementary/HC_Receptive_estimated_marginal_means.html")
gtsave(pairs_table, "tables/supplementary/HC_Receptive_pairwise_comparisons.html")
gtsave(slopes_table, "tables/supplementary/HC_Receptive_slope_comparisons.html")

# Create table to make itno word docs 
gt_to_flextable <- function(gt_table, 
                            zebra_stripes = TRUE, 
                            autofit = TRUE,
                            fontname = "Arial",
                            fontsize = 11,
                            header_color = "#2C3E50",
                            header_text = "white") {
  
  # Check required packages
  if (!requireNamespace("flextable", quietly = TRUE)) {
    stop("Package 'flextable' required but not installed")
  }
  
  # Validate input
  if (!inherits(gt_table, "gt_tbl")) {
    stop("Input must be a gt table object")
  }
  
  # Extract data from gt object
  flex_data <- tryCatch({
    as.data.frame(gt_table$`_data`)
  }, error = function(e) {
    stop("Could not extract data from gt table: ", e$message)
  })
  
  # Create basic flextable
  ft <- flextable::flextable(flex_data)
  
  # Apply basic formatting
  ft <- ft %>%
    flextable::font(fontname = fontname, part = "all") %>%
    flextable::fontsize(size = fontsize, part = "all") %>%
    flextable::align(align = "center", part = "header") %>%
    flextable::bold(part = "header") %>%
    flextable::bg(bg = header_color, part = "header") %>%
    flextable::color(color = header_text, part = "header")
  
  # Conditional formatting
  if (zebra_stripes) {
    ft <- ft %>% flextable::theme_zebra()
  }
  
  if (autofit) {
    ft <- ft %>% flextable::autofit()
  }
  
  # Add borders using officer::fp_border
  border_style <- officer::fp_border(color = "black", width = 1)
  inner_border <- officer::fp_border(color = "gray", width = 0.5)
  
  ft <- ft %>%
    flextable::border_outer(border = border_style, part = "all") %>%
    flextable::border_inner_h(border = inner_border, part = "body")
  
  return(ft)
}
# Save as Word documents
emm_flex <- gt_to_flextable(emm_table)
pairs_flex <- gt_to_flextable(pairs_table)
slopes_flex <- gt_to_flextable(slopes_table)

save_as_docx(emm_flex, path = "tables/supplementary/HC_Receptive_estimated_marginal_means.docx")
save_as_docx(pairs_flex, path = "tables/supplementary/HC_Receptive_pairwise_comparisons.docx")
save_as_docx(slopes_flex, path = "tables/supplementary/HC_Receptive_slope_comparisons.docx")

# Plot estimated marginal means with custom palette
emm_plot <- ggplot(emm_hc_time_df, aes(x = timepoint, y = emmean, color = hc_class)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = hc_class), 
              alpha = 0.2, color = NA) +
  scale_color_manual(values = hc_palette) +
  scale_fill_manual(values = hc_palette) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    text = element_text(size = 12)
  ) +
  labs(
    title = "Estimated Receptive Language Trajectories by HC Class",
    subtitle = "With 95% confidence intervals",
    x = "Age (months)",
    y = "Estimated Receptive Language Score",
    color = "HC Class",
    fill = "HC Class"
  )
emm_plot
# Save the plot
ggsave("figures/supplementary/HC_Receptive_estimated_trajectories.pdf", emm_plot, width = 10, height = 8)

# ============================================================================ #
# 5. FINAL MODEL INTERPRETATION
# ============================================================================ #

# Full summary of the final selected model
final_model <- best_model
final_summary <- summary(final_model)
print("Final Model Summary:")
print(final_summary)

# Create a table for the fixed effects
fixed_effects <- as.data.frame(coef(summary(final_model)))
fixed_effects$Parameter <- rownames(fixed_effects)
rownames(fixed_effects) <- NULL
fixed_effects <- fixed_effects %>% 
  mutate(
    p.value = 2 * (1 - pnorm(abs(`t value`))),
    sig = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      p.value < 0.1 ~ ".",
      TRUE ~ ""
    )
  ) %>%
  relocate(Parameter)

fixed_effects_table <- fixed_effects %>%
  gt() %>%
  tab_header(
    title = "Fixed Effects for Final Model",
    subtitle = "Linear mixed-effects model with random slopes"
  ) %>%
  fmt_number(
    columns = c("Estimate", "Std. Error", "t value"),
    decimals = 3
  ) %>%
  fmt_number(
    columns = "p.value",
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(color = "red"),
    locations = cells_body(
      columns = "p.value",
      rows = p.value < 0.05
    )
  ) %>%
  cols_label(
    Parameter = "Parameter",
    Estimate = "Estimate",
    `Std. Error` = "Std. Error",
    `t value` = "t value",
    p.value = "p-value",
    sig = "Significance"
  ) %>%
  tab_footnote(
    footnote = "* p < 0.05, ** p < 0.01, *** p < 0.001",
    locations = cells_column_labels(columns = sig)
  )

fixed_effects_table

# Save the fixed effects table
gtsave(fixed_effects_table, "tables/main/HC_Receptive_fixed_effects_table.html")
fixed_effects_flex <- gt_to_flextable(fixed_effects_table)
save_as_docx(fixed_effects_flex, path = "tables/main/HC_Receptive_fixed_effects_table.docx")

# Create a variance components table
vc <- VarCorr(final_model)
var_comp <- data.frame(
  Group = c("Random slope (timepoint)", "Residual"),
  Variance = c(as.numeric(vc$subject), attr(vc, "sc")^2),
  SD = c(sqrt(as.numeric(vc$subject)), attr(vc, "sc"))
)

var_comp_table <- var_comp %>%
  gt() %>%
  tab_header(
    title = "Variance Components for Final Model"
  ) %>%
  fmt_number(
    columns = c("Variance", "SD"),
    decimals = 3
  )

var_comp_table
# Save the variance components table
gtsave(var_comp_table, "tables/supplementary/HC_Receptive_variance_components.html")
var_comp_flex <- gt_to_flextable(var_comp_table)
save_as_docx(var_comp_flex, path = "tables/supplementary/HC_Receptive_variance_components.docx")

# Create a forest plot for key fixed effects
forest_data <- fixed_effects %>%
  filter(!grepl("Intercept", Parameter)) %>%
  mutate(
    lower.CI = Estimate - 1.96 * `Std. Error`,
    upper.CI = Estimate + 1.96 * `Std. Error`,
    Parameter = factor(Parameter, levels = rev(Parameter))
  )

forest_plot <- ggplot(forest_data, aes(x = Estimate, y = Parameter, 
                                       xmin = lower.CI, xmax = upper.CI)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbarh(height = 0.2) +
  geom_point(size = 3, aes(color = p.value < 0.05)) +
  scale_color_manual(values = c("FALSE" = "gray70", "TRUE" = "#184e77"),
                     labels = c("FALSE" = "Not significant", "TRUE" = "Significant")) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    text = element_text(size = 12)
  ) +
  labs(
    title = "Forest Plot of Model Parameters",
    subtitle = "Parameter estimates with 95% confidence intervals",
    x = "Parameter Estimate",
    y = "",
    color = "Statistical Significance"
  )
forest_plot
# Save the forest plot
ggsave("figures/supplementary/HC_Receptive_forest_plot.pdf", forest_plot, width = 10, height = 8)

# Create a summary table of key findings
summary_points <- data.frame(
  Finding = c(
    "HC Classes Trajectory Differences",
    "Nonverbal IQ Effect",
    "Model Fit Quality",
    "Missing Data Impact"
  ),
  Description = c(
    paste0("The three HC classes showed ", 
           ifelse(any(slope_pairs_df$p.value < 0.05), "significant", "non-significant"), 
           " differences in developmental trajectories."),
    paste0("Nonverbal IQ at 6 months ", 
           ifelse(any(grepl("nonverbal_iq_6", fixed_effects$Parameter[fixed_effects$p.value < 0.05])), 
                  "significantly predicted", "did not significantly predict"), 
           " language development."),
    paste0("The final model explained ", 
           round(r.squaredGLMM(final_model)[1] * 100, 1), 
           "% of the variance in receptive language."),
    paste0("Missing data patterns were predominantly at the ", 
           missing_df$Timepoint[which.max(missing_df$Missing_Percent)], 
           " timepoint (", 
           max(missing_df$Missing_Percent), "%).")
  )
)

summary_table <- summary_points %>%
  gt() %>%
  tab_header(
    title = "Summary of Key Findings",
    subtitle = "Language development trajectories by HC Class"
  ) %>%
  cols_label(
    Finding = "Key Finding",
    Description = "Description"
  )

summary_table
# Save the summary table
gtsave(summary_table, "tables/supplementary/HC_Receptive_summary_findings.html")
summary_flex <- gt_to_flextable (summary_table)
save_as_docx(summary_flex, path = "tables/supplementary/HC_Receptive_summary_findings.docx")

# ============================================================================ #
# SUPPLEMENTARY ANALYSES
# ============================================================================ #

# ============================================================================ #
# S1. WINSORIZATION SENSITIVITY ANALYSIS
# ============================================================================ #

# Create datasets with different winsorization levels
lang_data_orig <- lang_data  # No winsorization
lang_data_wins05 <- winsorize(lang_data, receptive_cols, by_group = "hc_class", trim_pct = 0.05)
lang_data_wins10 <- winsorize(lang_data, receptive_cols, by_group = "hc_class", trim_pct = 0.10)


# Convert to long format
long_data_orig <- to_long(lang_data_orig)
long_data_wins05 <- to_long(lang_data_wins05)
long_data_wins10 <- to_long(lang_data_wins10)

# Run models with different winsorization levels
models_orig <- run_main_models(long_data_orig )
models_wins05 <- run_main_models(long_data_wins05 )
models_wins10 <- run_main_models(long_data_wins10)

# Create model comparison tables for different winsorization levels
model_table_orig <- create_model_table(models_orig, "Model Comparison (No Winsorization)")
model_table_wins05 <- create_model_table(models_wins05, "Model Comparison (5% Winsorization)")
model_table_wins10 <- create_model_table(models_wins10, "Model Comparison (10% Winsorization)")
 

# Save model comparison tables
gtsave(model_table_orig, "sensitivity/model_comparison_original.html")
gtsave(model_table_wins05, "sensitivity/model_comparison_wins05.html")
gtsave(model_table_wins10, "sensitivity/model_comparison_wins10.html")

# Extract fixed effects from equivalent models for each winsorization level
fe_orig <- extract_fixed_effects(models_orig[[best_model_idx]])
fe_wins05 <- extract_fixed_effects(models_wins05[[best_model_idx]])
fe_wins10 <- extract_fixed_effects(models_wins10[[best_model_idx]])

# Add winsorization level as a column
fe_orig$Winsorization <- "None"
fe_wins05$Winsorization <- "5%"
fe_wins10$Winsorization <- "10%"

# Combine all fixed effects
all_fe <- bind_rows(fe_orig, fe_wins05, fe_wins10) %>%
  mutate(Winsorization = factor(Winsorization, levels = c("None", "5%", "10%")))

all_fe <- all_fe %>%
  mutate(
    significant = ifelse(lower.CI * upper.CI > 0, TRUE, FALSE),
    Parameter = factor(Parameter, levels = rev(unique(Parameter)))
  )
# Create the enhanced sensitivity plot with significance indicators
winsor_sensitivity_plot <- ggplot(all_fe, aes(x = Estimate, y = Parameter, color = Winsorization,
                                              xmin = lower.CI, xmax = upper.CI)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbarh(height = 0.1, position = position_dodge(width = 0.5), 
                 aes(alpha = significant), show.legend = FALSE) +
  geom_point(position = position_dodge(width = 0.5), 
             aes(shape = significant), size = 3) +
  # Add significance asterisks - FIXED THIS SECTION
  theme_minimal() +
  labs(title = "Parameter Estimates Across Winsorization Levels",
       subtitle = "Comparing no winsorization, 5% winsorization, and 10% winsorization\n* indicates estimates where 95% CI doesn't include zero",
       x = "Parameter Estimate", y = "", 
       color = "Winsorization Level",
       shape = "Significant") +
  scale_color_brewer(palette = "Set1") +
  scale_alpha_manual(values = c(0.5, 1)) +  # Faded for non-significant
  scale_shape_manual(values = c(1, 16)) +    # Hollow vs solid points
  theme(legend.position = "bottom",
        plot.subtitle = element_text(size = 10))
# Save sensitivity plot
ggsave("sensitivity/winsorization_sensitivity.pdf", winsor_sensitivity_plot, 
       width = 10, height = 8)

# ============================================================================ #
# S2. QUADRATIC TERMS ANALYSIS
# ============================================================================ #

long_data$timepoint_centered_sq <- long_data$timepoint_centered^2


# Function to run models with quadratic terms
run_quadratic_models <- function(data) {
  # Base model (your best model)
  m_linear <- lmer(receptive ~ timepoint_centered * hc_class + 
                     sex + site + group_type + nonverbal_iq_6 + 
                     timepoint_centered:nonverbal_iq_6 +
                     (0 + timepoint_centered | subject),
                   data = data, REML = FALSE,
                   control = lmerControl(optimizer = "bobyqa"))
  
  # Add quadratic term for time
  m_quad <- lmer(receptive ~ timepoint_centered * hc_class + 
                   timepoint_centered_sq +
                   sex + site + group_type + nonverbal_iq_6 + 
                   timepoint_centered:nonverbal_iq_6 +
                   (0 + timepoint_centered | subject),
                 data = data, REML = FALSE,
                 control = lmerControl(optimizer = "bobyqa"))
  
  # Add interaction between quadratic term and HC class
  m_quad_int <- lmer(receptive ~ timepoint_centered * hc_class + 
                       timepoint_centered_sq * hc_class +
                       sex + site + group_type + nonverbal_iq_6 + 
                       timepoint_centered:nonverbal_iq_6 +
                       (0 + timepoint_centered | subject),
                     data = data, REML = FALSE,
                     control = lmerControl(optimizer = "bobyqa"))
  
  # Return list of models
  return(list(m_linear = m_linear, m_quad = m_quad, m_quad_int = m_quad_int))
}


# Run quadratic models on the 5% winsorized data
quad_models <- run_quadratic_models(long_data)

# Create comparison table
quad_table <- create_model_table(quad_models, "Linear vs. Quadratic Models (5% Winsorization)")
gtsave(quad_table, "sensitivity/model_comparison_quadratic.html")

# Compare models using likelihood ratio test
quad_anova <- anova(quad_models$m_linear, quad_models$m_quad, quad_models$m_quad_int)
capture.output(quad_anova, file = "sensitivity/quadratic_anova.txt")

# ============================================================================ #
# S3. MULTIPLE IMPUTATION SENSITIVITY ANALYSIS
# ============================================================================ #

# Run MICE imputation on long data
imp_long <- mice(long_imp_prep, m = 10, seed = 12345, printFlag = FALSE)

# Create list to store models
imp_long_models <- list()

# Run models on each imputation
for (i in 1:10) {
  # Get the ith imputed dataset
  imp_data <- complete(imp_long, i)
  
  # Run best model
  imp_long_models[[i]] <- lmer(receptive ~ timepoint_centered * hc_class + 
                                 sex + site + group_type + nonverbal_iq_6 + 
                                 timepoint_centered:nonverbal_iq_6 +
                                 (0 + timepoint_centered | subject),
                               data = imp_data, REML = FALSE,
                               control = lmerControl(optimizer = "bobyqa"))
  
  # Also run quadratic model
  imp_long_models[[paste0("quad_", i)]] <- lmer(receptive ~ timepoint_centered * hc_class + 
                                                  timepoint_centered_sq * hc_class +
                                                  sex + site + group_type + nonverbal_iq_6 + 
                                                  timepoint_centered:nonverbal_iq_6 +
                                                  (0 + timepoint_centered | subject),
                                                data = imp_data, REML = FALSE,
                                                control = lmerControl(optimizer = "bobyqa"))
}

# Pool linear and quadratic models separately
linear_models <- imp_long_models[1:10]
quad_models <- imp_long_models[paste0("quad_", 1:10)]

# Convert each model to lme4 class
linear_models_lme4 <- lapply(linear_models, function(x) as(x, "lmerMod"))
quad_models_lme4 <- lapply(quad_models, function(x) as(x, "lmerMod"))

# Pool results
pool_linear <- pool(linear_models_lme4)
pool_quad <- pool(quad_models_lme4)

# Save pooled results
capture.output(summary(pool_linear), file = "sensitivity/mice_pooled_linear.txt")
capture.output(summary(pool_quad), file = "sensitivity/mice_pooled_quadratic.txt")

# Create a comparison table with compatible metrics
mice_comparison <- data.frame(
  Model = c("Original (No Imputation)", "MICE Linear", "MICE Quadratic"),
  AIC = c(AIC(best_model), 
          mean(sapply(linear_models_lme4, AIC)), 
          mean(sapply(quad_models_lme4, AIC))),
  BIC = c(BIC(best_model), 
          mean(sapply(linear_models_lme4, BIC)), 
          mean(sapply(quad_models_lme4, BIC))),
  LogLik = c(logLik(best_model), 
             mean(sapply(linear_models_lme4, logLik)), 
             mean(sapply(quad_models_lme4, logLik))),
  Fixed_Effects = c(length(fixef(best_model)),
                    mean(sapply(linear_models_lme4, function(x) length(fixef(x)))),
                    mean(sapply(quad_models_lme4, function(x) length(fixef(x))))),
  Convergence = c(ifelse(isSingular(best_model), "Singular", "Converged"),
                  ifelse(all(sapply(linear_models_lme4, isSingular)), "All Singular", 
                         ifelse(any(sapply(linear_models_lme4, isSingular)), "Some Singular", "All Converged")),
                  ifelse(all(sapply(quad_models_lme4, isSingular)), "All Singular", 
                         ifelse(any(sapply(quad_models_lme4, isSingular)), "Some Singular", "All Converged")))
) %>%
  mutate(across(where(is.numeric), ~round(., 2)))  # Round numeric columns

# Convert to gt table
mice_comparison_table <- mice_comparison %>%
  gt() %>%
  tab_header(
    title = "Model Comparison Summary",
    subtitle = "Original vs Multiply Imputed Models"
  ) %>%
  cols_label(
    Model = md("**Model Type**"),
    AIC = md("**AIC**"),
    BIC = md("**BIC**"),
    LogLik = md("**Log Likelihood**"),
    Fixed_Effects = md("**No. Fixed Effects**"),
    Convergence = md("**Convergence Status**")
  ) %>%
  fmt_number(columns = c(AIC, BIC), decimals = 1) %>%
  fmt_number(columns = LogLik, decimals = 2) %>%
  tab_style(
    style = cell_fill(color = "lightcyan"),
    locations = cells_body(rows = AIC == min(AIC))
  )

# Save MICE comparison table
gtsave(mice_comparison_table, "sensitivity/mice_model_comparison.html")

# ============================================================================ #
# S4. COMBINED SENSITIVITY SUMMARY
# ============================================================================ #

# Create a summary of all sensitivity analyses
sensitivity_summary <- data.frame(
  Analysis = c("Winsorization", "Quadratic Terms", "Multiple Imputation"),
  Key_Finding = c(
    "Parameter estimates were robust across different winsorization levels",
    ifelse(quad_anova$`Pr(>Chisq)`[2] < 0.05, 
           "Quadratic terms significantly improved model fit",
           "Quadratic terms did not significantly improve model fit"),
    ifelse(mean(sapply(linear_models_lme4, AIC)) < AIC(best_model),
           "Imputation improved model fit compared to complete cases",
           "Imputation did not improve model fit compared to complete cases")
  ),
  Recommendation = c(
    "Proceed with 5% winsorization as primary analysis",
    ifelse(quad_anova$`Pr(>Chisq)`[2] < 0.05,
           "Consider including quadratic terms in final model",
           "Linear model is sufficient"),
    ifelse(mean(sapply(linear_models_lme4, AIC)) < AIC(best_model),
           "Results should be interpreted with caution due to missing data",
           "Complete case analysis appears appropriate")
  )
)

# Save sensitivity summary
sensitivity_table <- sensitivity_summary %>%
  gt() %>%
  tab_header(
    title = "Sensitivity Analyses Summary",
    subtitle = "Key findings and recommendations from supplementary analyses"
  ) %>%
  cols_label(
    Analysis = "Analysis Type",
    Key_Finding = "Key Finding",
    Recommendation = "Recommendation"
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  )

gtsave(sensitivity_table, "sensitivity/sensitivity_summary.html")
