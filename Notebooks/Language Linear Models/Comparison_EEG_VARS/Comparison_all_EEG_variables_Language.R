# ============================================================================= #

# Comparison of EEG variables vs subtyping

# ============================================================================= #


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
library(officer)         # For saving tables as Word documents
library(e1071)  # For kurtosis and skewness


# Function to create better table 

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

# Load the data
lang_data <- read.csv("C:/Users/gabot/OneDrive - McGill University/Desktop/Github_repos/infant_neurosubs/Datasets/neurosub_solutions/neurosubs_df_2025_04_23_5D.csv", header = TRUE, sep = ",")

# Change the "" in outcome to nan
lang_data[lang_data == ""] <- NA

# Set current working directory as the same where this file is
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Select Palettes
hc_palette <- c("#184e77","#7570b3", "#1b9e77")

# Create output directories
dir.create("figures/main", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/supplementary", recursive = TRUE, showWarnings = FALSE)
dir.create("tables/supplementary", recursive = TRUE, showWarnings = FALSE)
dir.create("tables/main", recursive = TRUE, showWarnings = FALSE)



# Check data structure
print("Data Structure Summary:")
print(paste("Number of participants:", nrow(lang_data)))
print(paste("Variables in dataset:", ncol(lang_data)))

# Get column names for expressive language measures at different timepoints
expressive_cols <- grep("expressive_", names(lang_data), value = TRUE)
receptive_cols <- grep("receptive_", names(lang_data), value = TRUE)


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
lang_data_wins <- winsorize(lang_data, expressive_cols, by_group = "hc_class", trim_pct = 0.05)
lang_data_wins_recep <- winsorize(lang_data, receptive_cols, by_group = "hc_class", trim_pct = 0.05)

# Convert wide to long format for the main analysis dataset
to_long <- function(data) {
  long <- data %>%
    pivot_longer(
      cols = all_of(expressive_cols),
      names_to = "timepoint_var",
      values_to = "expressive"
    ) %>%
    mutate(
      # Extract numeric months from column names (assuming format expressive_Xmo)
      timepoint = as.numeric(gsub("expressive_", "", gsub("mo", "", timepoint_var))),
      # Center timepoint at first measurement (6mo)
      timepoint_centered = timepoint - 6,
      # Create categorical version for certain analyses
      timepoint_factor = factor(timepoint)
    )
  
  # Convert factors properly
  long$hc_class <- factor(long$hc_class)
  long$lpa_class <- factor(long$lpa_class)
  
  long$site <- factor(long$site)
  long$sex <- factor(long$sex)
  long$group_type <- factor(long$group_type)
  long$subject <- factor(long$subject)
  
  return(long)
}

# Create main analysis long dataset
long_data <- to_long(lang_data_wins)



# ============================================================================ #
# 3. MIXED-EFFECTS MODELING - MAIN ANALYSIS
# ============================================================================ #


# ============================================================================ #
# Enhanced EEG Features Comparison Analysis
# ============================================================================ #

# Define data to be used
main_data <- long_data

# EEG variables including clustering methods
eeg_variables <- c("front_gamma_6", "auditory_con_6", "lang_comp_con_6", 
                   "speech_con_6_left", "gamma_lat_6", "hc_class", "lpa_class")

# == EXPRESSIVE LANGUAGE MODEL ==

# Function to run individual feature models with p-value extraction
run_individual_feature_models <- function(data,language_var) {
  individual_models <- list()
  model_summaries <- list()
  
  for(var in eeg_variables) {
    # Handle factor variables differently
    if(var %in% c("hc_class", "lpa_class")) {
      # Model with categorical predictors
      formula <- as.formula(paste0(language_var, " ~ timepoint_centered * ", var, 
                                   " + sex + site + group_type + nonverbal_iq_6 + 
                                   (0 + timepoint_centered | subject)"))
      
      
      # 5. Add covariates (sex, site, risk, iq)
      m4 <- lmer(expressive ~ timepoint_centered * lpa_class + 
                   sex + site + group_type + nonverbal_iq_6+
                   (0 + timepoint_centered | subject),  # Random slopes only
                 data = data, REML = FALSE,
                 control = lmerControl(optimizer = "bobyqa"))
      print(var)
      print(all.equal(formula, formula(m4)))
      
      
      m_individual <- lmer(formula, data = data, REML = FALSE,
                           control = lmerControl(optimizer = "bobyqa"))
      
    } else {
      # Model with continuous predictors
      m_individual <- lmer(expressive ~ timepoint_centered * data[[var]] +
                             sex + site + group_type + nonverbal_iq_6 +
                             (0 + timepoint_centered | subject),
                           data = data, REML = FALSE,
                           control = lmerControl(optimizer = "bobyqa"))
    }
    
    individual_models[[var]] <- m_individual
    
    # Extract model summary
    model_sum <- summary(m_individual)
    model_summaries[[var]] <- model_sum
  }
  
  return(list(models = individual_models, summaries = model_summaries))
}

# Run individual feature models
model_results <- run_individual_feature_models(main_data, "expressive")
individual_models <- model_results$models

individual_models
model_summaries <- model_results$summaries

model_summaries
# Create enhanced comparison table with p-values

create_enhanced_comparison_table <- function(individual_models, model_summaries) {
  # Prepare comparison data
  comparison_data <- data.frame(
    Feature = names(individual_models),
    stringsAsFactors = FALSE
  )
  
  # Extract model fit metrics
  comparison_data$AIC <- sapply(individual_models, function(m) {
    tryCatch(AIC(m), error = function(e) NA_real_)
  })
  comparison_data$BIC <- sapply(individual_models, function(m) {
    tryCatch(BIC(m), error = function(e) NA_real_)
  })
  
  # Safely extract R-squared values
  r2_results <- tryCatch({
    sapply(individual_models, function(m) r.squaredGLMM(m))
  }, error = function(e) {
    warning("Error in r.squaredGLMM calculation")
    matrix(NA_real_, nrow = 2, ncol = length(individual_models))
  })
  
  comparison_data$R2_Marginal <- r2_results[1, ]
  comparison_data$R2_Conditional <- r2_results[2, ]
  
  # Extract p-values for main effects and interactions
  comparison_data$Main_p_value <- NA_real_
  comparison_data$Interaction_p_value <- NA_real_
  
  for(var in names(model_summaries)) {
    summ <- model_summaries[[var]]
    coef_table <- coef(summ)
    model <- individual_models[[var]]
    
    if(var %in% c("hc_class", "lpa_class")) {
      # For categorical variables, use previous approach
      main_effect_rows <- grep(paste0("^", var), rownames(coef_table), value = TRUE)
      interaction_rows <- grep(paste0("timepoint_centered:", var), rownames(coef_table), value = TRUE)
      
      # Get minimum p-value (most significant effect)
      if(length(main_effect_rows) > 0) {
        p_values <- sapply(main_effect_rows, function(row) {
          t_val <- coef_table[row, "t value"]
          2 * pt(abs(t_val), df = df.residual(model), lower.tail = FALSE)
        })
        comparison_data$Main_p_value[comparison_data$Feature == var] <- min(p_values, na.rm = TRUE)
      }
      
      if(length(interaction_rows) > 0) {
        p_values <- sapply(interaction_rows, function(row) {
          t_val <- coef_table[row, "t value"]
          2 * pt(abs(t_val), df = df.residual(model), lower.tail = FALSE)
        })
        comparison_data$Interaction_p_value[comparison_data$Feature == var] <- min(p_values, na.rm = TRUE)
      }
    } else {
      # For continuous variables, use exact literal string match
      
      # For main effect - look for literal "data[[var]]"
      main_effect_row <- which(rownames(coef_table) == "data[[var]]")
      if(length(main_effect_row) > 0) {
        t_val <- coef_table[main_effect_row, "t value"]
        p_val <- 2 * pt(abs(t_val), df = df.residual(model), lower.tail = FALSE)
        comparison_data$Main_p_value[comparison_data$Feature == var] <- p_val
      }
      
      # For interaction effect - look for literal "timepoint_centered:data[[var]]"
      interaction_row <- which(rownames(coef_table) == "timepoint_centered:data[[var]]")
      if(length(interaction_row) > 0) {
        t_val <- coef_table[interaction_row, "t value"]
        p_val <- 2 * pt(abs(t_val), df = df.residual(model), lower.tail = FALSE)
        comparison_data$Interaction_p_value[comparison_data$Feature == var] <- p_val
      }
    }
  }
  
  # Replace any NaN or Inf values with NA
  comparison_data$Main_p_value[is.nan(comparison_data$Main_p_value) | 
                                 is.infinite(comparison_data$Main_p_value)] <- NA_real_
  comparison_data$Interaction_p_value[is.nan(comparison_data$Interaction_p_value) | 
                                        is.infinite(comparison_data$Interaction_p_value)] <- NA_real_
  
  # Add significance indicators
  comparison_data$Main_Significant <- ifelse(!is.na(comparison_data$Main_p_value) & 
                                               comparison_data$Main_p_value < 0.05, "Yes", "No")
  comparison_data$Interaction_Significant <- ifelse(!is.na(comparison_data$Interaction_p_value) & 
                                                      comparison_data$Interaction_p_value < 0.05, "Yes", "No")
  
  # Format p-values
  comparison_data$Main_p_value_formatted <- NA_character_
  comparison_data$Interaction_p_value_formatted <- NA_character_
  
  # Only format valid p-values
  valid_main_p <- !is.na(comparison_data$Main_p_value)
  if(any(valid_main_p)) {
    comparison_data$Main_p_value_formatted[valid_main_p] <- 
      sprintf("%.3f", comparison_data$Main_p_value[valid_main_p])
  }
  
  valid_int_p <- !is.na(comparison_data$Interaction_p_value)
  if(any(valid_int_p)) {
    comparison_data$Interaction_p_value_formatted[valid_int_p] <- 
      sprintf("%.3f", comparison_data$Interaction_p_value[valid_int_p])
  }
  
  # Add asterisks for significance
  comparison_data$Main_p_value_formatted <- ifelse(
    !is.na(comparison_data$Main_p_value),
    paste0(
      comparison_data$Main_p_value_formatted,
      ifelse(comparison_data$Main_p_value < 0.001, "***",
             ifelse(comparison_data$Main_p_value < 0.01, "**",
                    ifelse(comparison_data$Main_p_value < 0.05, "*", "")))
    ),
    "NA"
  )
  
  comparison_data$Interaction_p_value_formatted <- ifelse(
    !is.na(comparison_data$Interaction_p_value),
    paste0(
      comparison_data$Interaction_p_value_formatted,
      ifelse(comparison_data$Interaction_p_value < 0.001, "***",
             ifelse(comparison_data$Interaction_p_value < 0.01, "**",
                    ifelse(comparison_data$Interaction_p_value < 0.05, "*", "")))
    ),
    "NA"
  )
  
  # Create formatted table
  formatted_table <- comparison_data %>%
    dplyr::select(Feature, AIC, BIC, R2_Marginal, R2_Conditional, 
           Main_p_value_formatted, Interaction_p_value_formatted) %>%
    gt() %>%
    tab_header(
      title = "Model Comparison: Subgrouping vs Individual EEG Features",
      subtitle = "Predictive Power for Expressive Language Trajectories"
    ) %>%
    fmt_number(
      columns = c("AIC", "BIC", "R2_Marginal", "R2_Conditional"),
      decimals = 3
    ) %>%
    cols_label(
      Feature = "Model/Feature",
      AIC = "AIC",
      BIC = "BIC", 
      R2_Marginal = "Marginal R²",
      R2_Conditional = "Conditional R²",
      Main_p_value_formatted = "Main Effect p-value",
      Interaction_p_value_formatted = "Interaction with Time p-value"
    ) %>%
    tab_style(
      style = cell_fill(color = "#e8f4f8"),
      locations = cells_body(
        rows = AIC == min(AIC, na.rm = TRUE)
      )
    ) %>%
    tab_style(
      style = cell_text(color = "red", weight = "bold"),
      locations = cells_body(
        columns = "Main_p_value_formatted",
        rows = grepl("[*]", Main_p_value_formatted)
      )
    ) %>%
    tab_style(
      style = cell_text(color = "red", weight = "bold"),
      locations = cells_body(
        columns = "Interaction_p_value_formatted",
        rows = grepl("[*]", Interaction_p_value_formatted)
      )
    ) %>%
    tab_footnote(
      footnote = "* p < 0.05, ** p < 0.01, *** p < 0.001",
      locations = cells_column_labels(columns = c("Main_p_value_formatted", "Interaction_p_value_formatted"))
    ) %>%
    tab_footnote(
      footnote = "NA values indicate the effect could not be reliably estimated",
      locations = cells_column_labels(columns = c("Main_p_value_formatted", "Interaction_p_value_formatted"))
    )
  
  return(list(
    data = comparison_data,
    table = formatted_table
  ))
}

# Generate enhanced comparison table
enhanced_comparison <- create_enhanced_comparison_table(individual_models, model_summaries)
enhanced_comparison_table <- enhanced_comparison$table
comparison_data <- enhanced_comparison$data
comparison_data

enhanced_comparison_table
comparison_data
# Save comparison table
gtsave(enhanced_comparison_table, "tables/main/Expressive_EEG_Feature_Comparison.html")
enhanced_comparison_flex <- gt_to_flextable(enhanced_comparison_table)
save_as_docx(enhanced_comparison_flex, path = "tables/main/Expressive_Enhanced_EEG_Feature_Comparison.docx")

# Create AIC comparison plot
aic_plot <- ggplot(comparison_data, aes(x = reorder(Feature, -AIC), y = AIC)) +
  geom_bar(stat = "identity", fill = "#184e77") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Expressive language Model Comparison by AIC",
    subtitle = "Lower values indicate better model fit",
    x = "EEG Feature",
    y = "Akaike Information Criterion (AIC)"
  )

aic_plot
# Save AIC plot
ggsave("figures/main/expressive_EEG_AIC_Comparison.pdf", aic_plot, width = 10, height = 6)

# Summary table of best predictors
best_predictors <- comparison_data %>%
  arrange(AIC) %>%
  dplyr::select(Feature, AIC, R2_Marginal, R2_Conditional, 
         Main_p_value, Interaction_p_value) %>%
  mutate(
    Main_Significant = Main_p_value < 0.05,
    Interaction_Significant = Interaction_p_value < 0.05,
    Overall_Ranking = row_number()
  )
best_predictors

# Create summary table
summary_table <- best_predictors %>%
  dplyr::select(Feature, Overall_Ranking, Main_Significant, Interaction_Significant) %>%
  gt() %>%
  tab_header(
    title = "Ranking of EEG Features as Predictors of Expressive Language Trajectories",
    subtitle = "Based on model fit and significance"
  ) %>%
  cols_label(
    Feature = "EEG Feature",
    Overall_Ranking = "Rank (by AIC)",
    Main_Significant = "Significant Main Effect",
    Interaction_Significant = "Significant Time Interaction"
  ) %>%
  tab_style(
    style = cell_fill(color = "#e8f4f8"),
    locations = cells_body(
      rows = Overall_Ranking <= 3
    )
  )

summary_table
# Save summary table
gtsave(summary_table, "tables/main/Expressive_EEG_Feature_Rankings.html")
summary_flex <- gt_to_flextable(summary_table)
save_as_docx(summary_flex, path = "tables/main/Expressive_EEG_Feature_Rankings.docx")




# == RECEPTIVE LANGUAGE MODEL ==


# Convert wide to long format for the main analysis dataset
to_long_receptive <- function(data) {
  long <- data %>%
    pivot_longer(
      cols = all_of(receptive_cols),
      names_to = "timepoint_var",
      values_to = "receptive"
    ) %>%
    mutate(
      # Extract numeric months from column names (assuming format expressive_Xmo)
      timepoint = as.numeric(gsub("receptive_", "", gsub("mo", "", timepoint_var))),
      # Center timepoint at first measurement (6mo)
      timepoint_centered = timepoint - 6,
      # Create categorical version for certain analyses
      timepoint_factor = factor(timepoint)
    )
  
  # Convert factors properly
  long$hc_class <- factor(long$hc_class)
  long$lpa_class <- factor(long$lpa_class)
  
  long$site <- factor(long$site)
  long$sex <- factor(long$sex)
  long$group_type <- factor(long$group_type)
  long$subject <- factor(long$subject)
  
  return(long)
}

# Create main analysis long dataset
long_data <- to_long_receptive(lang_data_wins_recep)

long_data
main_data <- long_data

#Function to run individual feature models with p-value extraction
run_individual_feature_models <- function(data,language_var) {
  individual_models <- list()
  model_summaries <- list()
  
  for(var in eeg_variables) {
    # Handle factor variables differently
    if(var %in% c("hc_class", "lpa_class")) {
      # Model with categorical predictors
      formula <- as.formula(paste0(language_var, " ~ timepoint_centered * ", var, 
                                   " + sex + site + group_type + nonverbal_iq_6 + 
                                   (0 + timepoint_centered | subject)"))
      
      m_individual <- lmer(formula, data = data, REML = FALSE,
                           control = lmerControl(optimizer = "bobyqa"))
    } else {
      # Model with continuous predictors
      m_individual <- lmer(receptive ~ timepoint_centered * data[[var]] +
                             sex + site + group_type + nonverbal_iq_6 +
                             (0 + timepoint_centered | subject),
                           data = data, REML = FALSE,
                           control = lmerControl(optimizer = "bobyqa"))
    }
    
    individual_models[[var]] <- m_individual
    
    # Extract model summary
    model_sum <- summary(m_individual)
    model_summaries[[var]] <- model_sum
  }
  
  return(list(models = individual_models, summaries = model_summaries))
}

main_data
# Run individual feature models
model_results <- run_individual_feature_models(main_data, "receptive")
individual_models <- model_results$models

individual_models
model_summaries <- model_results$summaries

model_summaries
# Create enhanced comparison table with p-values

create_enhanced_comparison_table <- function(individual_models, model_summaries) {
  # Prepare comparison data
  comparison_data <- data.frame(
    Feature = names(individual_models),
    stringsAsFactors = FALSE
  )
  
  # Extract model fit metrics
  comparison_data$AIC <- sapply(individual_models, function(m) {
    tryCatch(AIC(m), error = function(e) NA_real_)
  })
  comparison_data$BIC <- sapply(individual_models, function(m) {
    tryCatch(BIC(m), error = function(e) NA_real_)
  })
  
  # Safely extract R-squared values
  r2_results <- tryCatch({
    sapply(individual_models, function(m) r.squaredGLMM(m))
  }, error = function(e) {
    warning("Error in r.squaredGLMM calculation")
    matrix(NA_real_, nrow = 2, ncol = length(individual_models))
  })
  
  comparison_data$R2_Marginal <- r2_results[1, ]
  comparison_data$R2_Conditional <- r2_results[2, ]
  
  # Extract p-values for main effects and interactions
  comparison_data$Main_p_value <- NA_real_
  comparison_data$Interaction_p_value <- NA_real_
  
  for(var in names(model_summaries)) {
    summ <- model_summaries[[var]]
    coef_table <- coef(summ)
    model <- individual_models[[var]]
    
    if(var %in% c("hc_class", "lpa_class")) {
      # For categorical variables, use previous approach
      main_effect_rows <- grep(paste0("^", var), rownames(coef_table), value = TRUE)
      interaction_rows <- grep(paste0("timepoint_centered:", var), rownames(coef_table), value = TRUE)
      
      # Get minimum p-value (most significant effect)
      if(length(main_effect_rows) > 0) {
        p_values <- sapply(main_effect_rows, function(row) {
          t_val <- coef_table[row, "t value"]
          2 * pt(abs(t_val), df = df.residual(model), lower.tail = FALSE)
        })
        comparison_data$Main_p_value[comparison_data$Feature == var] <- min(p_values, na.rm = TRUE)
      }
      
      if(length(interaction_rows) > 0) {
        p_values <- sapply(interaction_rows, function(row) {
          t_val <- coef_table[row, "t value"]
          2 * pt(abs(t_val), df = df.residual(model), lower.tail = FALSE)
        })
        comparison_data$Interaction_p_value[comparison_data$Feature == var] <- min(p_values, na.rm = TRUE)
      }
    } else {
      # For continuous variables, use exact literal string match
      
      # For main effect - look for literal "data[[var]]"
      main_effect_row <- which(rownames(coef_table) == "data[[var]]")
      if(length(main_effect_row) > 0) {
        t_val <- coef_table[main_effect_row, "t value"]
        p_val <- 2 * pt(abs(t_val), df = df.residual(model), lower.tail = FALSE)
        comparison_data$Main_p_value[comparison_data$Feature == var] <- p_val
      }
      
      # For interaction effect - look for literal "timepoint_centered:data[[var]]"
      interaction_row <- which(rownames(coef_table) == "timepoint_centered:data[[var]]")
      if(length(interaction_row) > 0) {
        t_val <- coef_table[interaction_row, "t value"]
        p_val <- 2 * pt(abs(t_val), df = df.residual(model), lower.tail = FALSE)
        comparison_data$Interaction_p_value[comparison_data$Feature == var] <- p_val
      }
    }
  }
  
  # Replace any NaN or Inf values with NA
  comparison_data$Main_p_value[is.nan(comparison_data$Main_p_value) | 
                                 is.infinite(comparison_data$Main_p_value)] <- NA_real_
  comparison_data$Interaction_p_value[is.nan(comparison_data$Interaction_p_value) | 
                                        is.infinite(comparison_data$Interaction_p_value)] <- NA_real_
  
  # Add significance indicators
  comparison_data$Main_Significant <- ifelse(!is.na(comparison_data$Main_p_value) & 
                                               comparison_data$Main_p_value < 0.05, "Yes", "No")
  comparison_data$Interaction_Significant <- ifelse(!is.na(comparison_data$Interaction_p_value) & 
                                                      comparison_data$Interaction_p_value < 0.05, "Yes", "No")
  
  # Format p-values
  comparison_data$Main_p_value_formatted <- NA_character_
  comparison_data$Interaction_p_value_formatted <- NA_character_
  
  # Only format valid p-values
  valid_main_p <- !is.na(comparison_data$Main_p_value)
  if(any(valid_main_p)) {
    comparison_data$Main_p_value_formatted[valid_main_p] <- 
      sprintf("%.3f", comparison_data$Main_p_value[valid_main_p])
  }
  
  valid_int_p <- !is.na(comparison_data$Interaction_p_value)
  if(any(valid_int_p)) {
    comparison_data$Interaction_p_value_formatted[valid_int_p] <- 
      sprintf("%.3f", comparison_data$Interaction_p_value[valid_int_p])
  }
  
  # Add asterisks for significance
  comparison_data$Main_p_value_formatted <- ifelse(
    !is.na(comparison_data$Main_p_value),
    paste0(
      comparison_data$Main_p_value_formatted,
      ifelse(comparison_data$Main_p_value < 0.001, "***",
             ifelse(comparison_data$Main_p_value < 0.01, "**",
                    ifelse(comparison_data$Main_p_value < 0.05, "*", "")))
    ),
    "NA"
  )
  
  comparison_data$Interaction_p_value_formatted <- ifelse(
    !is.na(comparison_data$Interaction_p_value),
    paste0(
      comparison_data$Interaction_p_value_formatted,
      ifelse(comparison_data$Interaction_p_value < 0.001, "***",
             ifelse(comparison_data$Interaction_p_value < 0.01, "**",
                    ifelse(comparison_data$Interaction_p_value < 0.05, "*", "")))
    ),
    "NA"
  )
  
  # Create formatted table
  formatted_table <- comparison_data %>%
    dplyr::select(Feature, AIC, BIC, R2_Marginal, R2_Conditional, 
                  Main_p_value_formatted, Interaction_p_value_formatted) %>%
    gt() %>%
    tab_header(
      title = "Model Comparison: Subgrouping vs Individual EEG Features",
      subtitle = "Predictive Power for Receptive Language Trajectories"
    ) %>%
    fmt_number(
      columns = c("AIC", "BIC", "R2_Marginal", "R2_Conditional"),
      decimals = 3
    ) %>%
    cols_label(
      Feature = "Model/Feature",
      AIC = "AIC",
      BIC = "BIC", 
      R2_Marginal = "Marginal R²",
      R2_Conditional = "Conditional R²",
      Main_p_value_formatted = "Main Effect p-value",
      Interaction_p_value_formatted = "Interaction with Time p-value"
    ) %>%
    tab_style(
      style = cell_fill(color = "#e8f4f8"),
      locations = cells_body(
        rows = AIC == min(AIC, na.rm = TRUE)
      )
    ) %>%
    tab_style(
      style = cell_text(color = "red", weight = "bold"),
      locations = cells_body(
        columns = "Main_p_value_formatted",
        rows = grepl("[*]", Main_p_value_formatted)
      )
    ) %>%
    tab_style(
      style = cell_text(color = "red", weight = "bold"),
      locations = cells_body(
        columns = "Interaction_p_value_formatted",
        rows = grepl("[*]", Interaction_p_value_formatted)
      )
    ) %>%
    tab_footnote(
      footnote = "* p < 0.05, ** p < 0.01, *** p < 0.001",
      locations = cells_column_labels(columns = c("Main_p_value_formatted", "Interaction_p_value_formatted"))
    ) %>%
    tab_footnote(
      footnote = "NA values indicate the effect could not be reliably estimated",
      locations = cells_column_labels(columns = c("Main_p_value_formatted", "Interaction_p_value_formatted"))
    )
  
  return(list(
    data = comparison_data,
    table = formatted_table
  ))
}

# Generate enhanced comparison table
enhanced_comparison <- create_enhanced_comparison_table(individual_models, model_summaries)
enhanced_comparison_table <- enhanced_comparison$table
comparison_data <- enhanced_comparison$data
comparison_data

enhanced_comparison_table
comparison_data
# Save comparison table
gtsave(enhanced_comparison_table, "tables/main/Receptive_EEG_Feature_Comparison.html")
enhanced_comparison_flex <- gt_to_flextable(enhanced_comparison_table)
save_as_docx(enhanced_comparison_flex, path = "tables/main/Receptive_Enhanced_EEG_Feature_Comparison.docx")

# Create AIC comparison plot
aic_plot <- ggplot(comparison_data, aes(x = reorder(Feature, -AIC), y = AIC)) +
  geom_bar(stat = "identity", fill = "#184e77") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Receptive language Model Comparison by AIC",
    subtitle = "Lower values indicate better model fit",
    x = "EEG Feature",
    y = "Akaike Information Criterion (AIC)"
  )

aic_plot
# Save AIC plot
ggsave("figures/main/receptive_EEG_AIC_Comparison.pdf", aic_plot, width = 10, height = 6)

# Summary table of best predictors
best_predictors <- comparison_data %>%
  arrange(AIC) %>%
  dplyr::select(Feature, AIC, R2_Marginal, R2_Conditional, 
                Main_p_value, Interaction_p_value) %>%
  mutate(
    Main_Significant = Main_p_value < 0.05,
    Interaction_Significant = Interaction_p_value < 0.05,
    Overall_Ranking = row_number()
  )
best_predictors

# Create summary table
summary_table <- best_predictors %>%
  dplyr::select(Feature, Overall_Ranking, Main_Significant, Interaction_Significant) %>%
  gt() %>%
  tab_header(
    title = "Ranking of EEG Features as Predictors of Receptive Language Trajectories",
    subtitle = "Based on model fit and significance"
  ) %>%
  cols_label(
    Feature = "EEG Feature",
    Overall_Ranking = "Rank (by AIC)",
    Main_Significant = "Significant Main Effect",
    Interaction_Significant = "Significant Time Interaction"
  ) %>%
  tab_style(
    style = cell_fill(color = "#e8f4f8"),
    locations = cells_body(
      rows = Overall_Ranking <= 3
    )
  )

summary_table
# Save summary table
gtsave(summary_table, "tables/main/Receptive_EEG_Feature_Rankings.html")
summary_flex <- gt_to_flextable(summary_table)
save_as_docx(summary_flex, path = "tables/main/Receptive_EEG_Feature_Rankings.docx")


model_summaries$lpa_class


# =========================================================================== #

# EEG Feature Network Connectivity Analysis

# =========================================================================== #
# Create correlation network
correlation_matrix <- cor(main_data[, eeg_features], use = "pairwise.complete.obs")

# Create igraph object
network <- graph_from_adjacency_matrix(
  abs(correlation_matrix) > 0.5,  # Threshold for connectivity
  mode = "undirected",
  weighted = TRUE
)

# Network metrics
network_metrics <- data.frame(
  Feature = eeg_features,
  Degree = degree(network),
  Betweenness = betweenness(network),
  Closeness = closeness(network)
)

# Visualize network
network_plot <- ggraph(network, layout = "fr") +
  geom_edge_link(aes(width = E(network)$weight), alpha = 0.5) +
  geom_node_point(size = 10, color = "#184e77") +
  geom_node_text(aes(label = name), repel = FALSE) +
  theme_minimal() +
  labs(
    title = "EEG Feature Connectivity Network",
    subtitle = "Nodes represent features, edges represent strong correlations"
  )
network_plot
# Save network plot
ggsave("figures/eeg_network_connectivity.pdf", network_plot, 
       width = 10, height = 8)





