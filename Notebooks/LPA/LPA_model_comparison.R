

#New function to compare across all variable sets
compare_all_models <- function(base_dir = "./", variable_sets = c("5D", "3D_con", "2D"), var_interest="Class") {
  all_results <- list()
  combined_fit <- data.frame()
  
  for (vs in variable_sets) {
    cat("\n===== Analyzing", vs, "models =====\n")
    results <- analyze_lpa_models(vs, base_dir, var_interest)
    all_results[[vs]] <- results
    
    # Add variable set identifier and combine fit stats
    fit_stats <- results$fit_comparison
    fit_stats$Variable_Set <- vs
    combined_fit <- bind_rows(combined_fit, fit_stats)
  }
  
  # Create comparison directory
  comp_dir <- file.path("C:/Users/gabot/OneDrive - McGill University/Desktop/Github_repos/infant_neurosubs/Tables/Mplus/", "Model_comparisons")
  if (!dir.exists(comp_dir)) dir.create(comp_dir)
  
  figure_dir=file.path("C:/Users/gabot/OneDrive - McGill University/Desktop/Github_repos/infant_neurosubs/Figures/Mplus/", "Model_comparisons")
  if (!dir.exists(figure_dir)) dir.create(figure_dir)
  
  # Save combined results
  write.csv(combined_fit, file.path(comp_dir, "Mplus_all_fit_statistics.csv"), row.names = FALSE)
  
  # Set the background to white for all my plots
  
  theme_set(
    theme_minimal() +
      theme(
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)
      )
  )
  
  
  # Create comparison plots
  p_combined <- ggplot(combined_fit, aes(x = Model, y = BIC, color = Variable_Set)) +
    geom_point(size = 3) +
    geom_line(aes(group = Variable_Set)) +
    labs(title = "BIC Comparison Across All Models",
         x = "Number of Classes",
         y = "BIC Value",
         color = "Variable Set")
  
  ggsave(file.path(figure_dir, "combined_bic_comparison.png"), p_combined, width = 10, height = 6)
  
  # Adjusted Bic 
  p2_combined <- ggplot(combined_fit, aes(x = Model, y = Adjusted_BIC, color = Variable_Set)) +
    geom_point(size = 3) +
    geom_line(aes(group = Variable_Set)) +
    labs(title = "Adjusted BIC Comparison Across All Models",
         x = "Number of Classes",
         y = "Adjusted-BIC Value",
         color = "Variable Set")
  
  ggsave(file.path(figure_dir, "combined_adj_bic_comparison.png"), p2_combined, width = 10, height = 6)
  
  
  # Adjusted AIC 
  p3_combined <- ggplot(combined_fit, aes(x = Model, y = AIC, color = Variable_Set)) +
    geom_point(size = 3) +
    geom_line(aes(group = Variable_Set)) +
    labs(title = "AIC Comparison Across All Models",
         x = "Number of Classes",
         y = "AIC Value",
         color = "Variable Set")
  
  ggsave(file.path(figure_dir, "combined_AIC_comparison.png"), p3_combined, width = 10, height = 6)
  
  
  
  # Entropy
  p4_combined <- ggplot(combined_fit, aes(x = Model, y = Entropy, color = Variable_Set)) +
    geom_point(size = 3) +
    geom_line(aes(group = Variable_Set)) +
    labs(title = "Entropy Comparison Across All Models",
         x = "Number of Classes",
         y = "Entropy Value",
         color = "Variable Set") 
  
  ggsave(file.path(figure_dir, "combined_entropy_comparison.png"), p4_combined, width = 10, height = 6)
  
  
  # Print and return results
  cat("\n===== Combined Model Comparison =====\n")
  print(combined_fit)
  print(p_combined)
  
  return(list(
    all_results = all_results,
    combined_fit = combined_fit,
    comparison_plot = p_combined
  ))
}


analyze_lpa_models <- function(variable_set, base_dir = "./", var_interest="class") {
  # Load required packages
  if (!require("MplusAutomation")) install.packages("MplusAutomation")
  if (!require("tidyverse")) install.packages("tidyverse")
  if (!require("gridExtra")) install.packages("gridExtra")
  library(MplusAutomation)
  library(tidyverse)
  library(gridExtra)
  
  # Helper functions
  check_replication <- function(model) {
    if (!is.null(model$replicated)) {
      return(ifelse(model$replicated, "Yes", "No"))
    } else {
      return("NA")
    }
  }
  
  calculate_small_classes <- function(class_counts, threshold = 0.05) {
    proportions <- class_counts$proportion
    return(sum(proportions < threshold))
  }
  
  format_class_counts <- function(class_counts, model_name) {
    data.frame(
      Model = model_name,
      Class = paste("Class", 1:nrow(class_counts)),
      Size = class_counts$count,
      Proportion = class_counts$proportion,
      stringsAsFactors = FALSE
    )
  }
  
  select_best_model <- function(fit_comp) {
    fit_comp %>%
      filter(Small_Classes == 0) %>%
      arrange(BIC, -Entropy) %>%
      slice(1)
  }
  
  # Create reports directory
  report_dir <- file.path(base_dir, variable_set, "reports")
  if (!dir.exists(report_dir)) {
    dir.create(report_dir, recursive = TRUE)
  }
  
  # Find ALL .out files in the directory
  out_files <- list.files(
    path = file.path(base_dir, variable_set),
    pattern = "\\.out$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  if (length(out_files) == 0) {
    stop("No .out files found in directory: ", file.path(base_dir, variable_set))
  }
  
  # Read all models and extract class numbers from file names
  models <- list()
  for (file in out_files) {
    # Extract class number (2-5) from filename
    class_num <- str_extract(basename(file), "(?<=_)[2-5](?=cls|_)")
    if (!is.na(class_num)) {
      models[[paste0(class_num, "_class")]] <- readModels(file)
    } else {
      # Alternative pattern if first one didn't match
      class_num <- str_extract(basename(file), "[2-5](?=classes|_)")
      if (!is.na(class_num)) {
        models[[paste0(class_num, "_class")]] <- readModels(file)
      }
    }
  }
  
  if (length(models) == 0) {
    stop("No valid models found in output files - could not extract class numbers")
  }
  
  # Extract fit statistics
  fit_comparison <- map_dfr(names(models), function(model_name) {
    model <- models[[model_name]]
    class_num <- str_extract(model_name, "[2-5]")
    
    data.frame(
      Variable_Set = variable_set,
      Model = paste0(class_num, "-Class"),
      AIC = model$summaries$AIC,
      BIC = model$summaries$BIC,
      Adjusted_BIC = model$summaries$aBIC,
      Entropy = model$summaries$Entropy,
      Replication_Status = check_replication(model),
      Small_Classes = calculate_small_classes(model$class_counts$modelEstimated),
      stringsAsFactors = FALSE
    )
  })
  
  # Extract class counts
  class_sizes <- map_dfr(names(models), function(model_name) {
    model <- models[[model_name]]
    class_num <- str_extract(model_name, "[2-5]")
    counts <- model$class_counts$modelEstimated
    data.frame(
      Model = paste0(class_num, "-Class"),
      Class = paste("Class", seq_len(nrow(counts))),
      Size = counts$count,
      Proportion = counts$proportion,
      stringsAsFactors = FALSE
    )
  })
  
  # Create wide format for class sizes
  if (nrow(class_sizes) > 0) {
    class_sizes_wide <- class_sizes %>%
      dplyr::select(Model, Class, Size) %>%
      tidyr::pivot_wider(names_from = Class, values_from = Size, names_prefix = "Class_")
    
    # Merge with fit comparison
    fit_comparison <- fit_comparison %>%
      left_join(class_sizes_wide, by = "Model")
  } else {
    warning("No class sizes found to merge")
  }
  
  
  # Print results to console
  cat("\n=== Fit Statistics Comparison ===\n")
  print(fit_comparison)
  
  cat("\n=== Class Sizes ===\n")
  print(class_sizes)
  
  # Generate and display plots
  p1 <- ggplot(fit_comparison, aes(x = Model)) +
    geom_point(aes(y = AIC, color = "AIC"), size = 3) +
    geom_point(aes(y = BIC, color = "BIC"), size = 3) +
    geom_point(aes(y = Adjusted_BIC, color = "Adjusted BIC"), size = 3) +
    labs(title = paste(variable_set, "Model Comparison: AIC, BIC, and Adjusted BIC"),
         x = "Model",
         y = "Value",
         color = "Fit Statistic") +
    theme_minimal()
  
  p2 <- ggplot(fit_comparison, aes(x = Model, y = Entropy)) +
    geom_line(size = 1, color = "red") +
    geom_point(size = 3, color = "red") +
    labs(title = paste(variable_set, "Model Comparison: Entropy"),
         x = "Model",
         y = "Entropy") +
    theme_minimal()
  
  # Display plots in R
  grid.arrange(p1, p2, ncol = 2)
  
  # Identify best model
  best_model <- select_best_model(fit_comparison)
  cat("\n=== Best Model ===\n")
  print(best_model)
  
  # Generate reports
  ## 1. Fit statistics table
  write.table(
    fit_comparison,
    file = file.path(report_dir, paste0(variable_set, "_fit_statistics.txt")),
    sep = "\t",
    row.names = FALSE
  )
  
  ## 2. Class sizes table
  write.table(
    class_sizes,
    file = file.path(report_dir, paste0(variable_set, "_class_sizes.txt")),
    sep = "\t",
    row.names = FALSE
  )
  
  ## 3. Plots
  ggsave(
    file.path(report_dir, paste0(variable_set, "_fit_statistics_plots.png")),
    arrangeGrob(p1, p2, ncol = 2),
    width = 12,
    height = 6
  )
  
  ## 4. Best model selection
  writeLines(
    c(paste("Best model for", variable_set, "variable set:"),
      capture.output(print(best_model))),
    file.path(report_dir, paste0(variable_set, "_best_model.txt"))
  )
  
  ## Save the variable of interest once again
  
  if (var_interest==variable_set) {
    
    var_interest_dir= "C:/Users/gabot/OneDrive - McGill University/Desktop/Github_repos/infant_neurosubs/Tables/Mplus"
    report_dir <- file.path(var_interest_dir)
    
    # Generate and display plots
    p1 <- ggplot(fit_comparison, aes(x = Model)) +
      geom_point(aes(y = AIC, color = "AIC"), size = 3) +
      geom_point(aes(y = BIC, color = "BIC"), size = 3) +
      geom_point(aes(y = Adjusted_BIC, color = "Adjusted BIC"), size = 3) +
      labs(title = paste("5 EEG Feature", "Model Comparison: AIC, BIC, and Adjusted BIC"),
           x = "Model",
           y = "Value",
           color = "Fit Statistic") +
      theme_minimal()
    
    p2 <- ggplot(fit_comparison, aes(x = Model, y = Entropy)) +
      geom_line(size = 1, color = "red") +
      geom_point(size = 3, color = "red") +
      labs(title = paste("5 EEG Features", "Model Comparison: Entropy"),
           x = "Model",
           y = "Entropy") +
      theme_minimal()
    
    
    
    ## 1.2 Fit statistics table
    write.csv(
      fit_comparison,
      file = file.path(report_dir, paste0(variable_set, "_fit_statistics.csv")),
      row.names = FALSE
    )
    
    ## 3.2 Plots
    ggsave(
      file.path("C:/Users/gabot/OneDrive - McGill University/Desktop/Github_repos/infant_neurosubs/Figures/Mplus", paste0(variable_set, "_LPA_fit_statistics_plots.png")),
      arrangeGrob(p1, p2, ncol = 2),
      width = 12,
      height = 6
    )
    
    ## 4.2 Best model selection
    writeLines(
      c(paste("Best model for", variable_set, "variable set:"),
        capture.output(print(best_model))),
      file.path(report_dir, paste0(variable_set, "_best_model.txt"))
    )
    
  }
  
  
  # Return all results
  return(list(
    fit_comparison = fit_comparison,
    class_sizes = class_sizes,
    best_model = best_model,
    report_dir = report_dir,
    plots = list(aic_bic_plot = p1, entropy_plot = p2)
  ))
}

results <- compare_all_models(base_dir = "C:/Users/gabot/OneDrive - McGill University/Desktop/Github_repos/infant_neurosubs/Datasets/mplus/May_2025/",
                              variable_sets = c("5D", "3D_con", "2D"), var_interest="5D")
