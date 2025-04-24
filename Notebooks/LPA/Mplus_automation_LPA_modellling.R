# Function to automate LPA analysis for classes 2-5 with multiple variable sets
run_comprehensive_lpa <- function(base_dir = "Datasets") {
  # Load required packages
  if (!require("MplusAutomation")) install.packages("MplusAutomation")
  library(MplusAutomation)
  
  # Define variable sets with shortened names (8 chars or less)
  variable_sets <- list(
    "5D" = c("aud_con6", "lang_c6", "spch_l6", "frnt_g6", "gma_lat6"),
    "3D_con" = c("aud_con6", "lang_c6", "spch_l6"),
    "2D" = c("gma_lat6", "spch_l6")
  )
  
  # Mapping of shortened names to original names for reference
  var_mapping <- list(
    "aud_con6" = "auditory_con_6",
    "lang_c6" = "lang_comp_con_6",
    "spch_l6" = "speech_con_6_left",
    "frnt_g6" = "front_gamma_6",
    "gma_lat6" = "gamma_lat_6"
  )
  
  # Generate timestamp for filenames
  timestamp <- format(Sys.Date(), "%Y%m%d") # Shorter format
  
  # Function to create input files for different class solutions
  create_lpa_input <- function(n_classes, vars, set_name, dat_path) {
    # Get just the filename without path
    dat_file_name <- basename(dat_path)
    
    # Create the variable sections
    names_section <- paste(vars, collapse = " ")
    usevars_section <- paste(vars, collapse = " ")
    
    # Create the complete template
    template <- sprintf(
      'TITLE: LPA %s %dcls %s

DATA:
  FILE IS "%s";

VARIABLE:
  NAMES ARE %s;
  
  USEVARIABLES ARE %s;
  MISSING ARE ALL (999);
  CLASSES = L(%d); 

ANALYSIS:
  TYPE IS MIXTURE;
  STARTS = 1500 100;
  STITERATIONS = 50;

MODEL:

OUTPUT:  
  SAMPSTAT;
  TECH11;

SAVEDATA:
  FILE IS "mplus_%s_%dcls_%s.csv";
  FORMAT IS FREE;
  RECORDLENGTH = 1000;
  SAVE = CPROB;
',
      substr(set_name, 1, 5), # Shorten set name
      n_classes,
      substr(timestamp, 3, 8), # Short date
      dat_file_name, # Just the filename
      names_section,
      usevars_section,
      n_classes,
      substr(set_name, 1, 5), # Shorten set name
      n_classes,
      substr(timestamp, 3, 8) # Short date
    )
    
    # Create input file in same directory as DAT file
    input_filename <- file.path(
      dirname(dat_path),
      sprintf("LPA_%s_%dcls_%s.inp", 
              substr(set_name, 1, 5), 
              n_classes, 
              substr(timestamp, 3, 8))
    )
    
    # Write file
    writeLines(template, input_filename)
    
    return(input_filename)
  }
  
  # Run analyses for each variable set
  all_results <- list()
  for (set_name in names(variable_sets)) {
    message("\nRunning analyses for variable set: ", set_name)
    set_results <- list()
    
    # Find the DAT file in the specific folder
    dat_path <- file.path(base_dir, set_name)
    dat_files <- list.files(path = dat_path, pattern = "\\.DAT$", ignore.case = TRUE, full.names = TRUE)
    
    if (length(dat_files) == 0) {
      warning("No .DAT file found in directory: ", dat_path)
      next
    }
    
    if (length(dat_files) > 1) {
      warning("Multiple .DAT files found in ", dat_path, ". Using the first one: ", dat_files[1])
    }
    
    dat_file <- dat_files[1]
    
    for (n in 2:5) {
      message("Running ", n, "-class solution...")
      
      # Create input file
      inp_file <- create_lpa_input(n, variable_sets[[set_name]], set_name, dat_file)
      
      # Run Mplus from the same directory
      wd_original <- getwd()
      setwd(dirname(dat_file))
      
      tryCatch({
        run <- runModels(basename(inp_file), showOutput = TRUE, logFile = NULL)
        set_results[[paste0(n, "_classes")]] <- readModels(basename(inp_file))
        message("Successfully completed ", n, "-class solution for ", set_name)
      }, error = function(e) {
        warning("Failed to run ", n, "-class solution for ", set_name, ": ", e$message)
      }, finally = {
        setwd(wd_original)
      })
    }
    
    all_results[[set_name]] <- set_results
  }
  
  # Return all results with variable mapping
  return(list(results = all_results, variable_mapping = var_mapping))
}



# Run the comprehensive analysis
lpa_results <- run_comprehensive_lpa(base_dir = "C:/Users/gabot/OneDrive - McGill University/Desktop/Github_repos/infant_neurosubs/Datasets/mplus/May_2025")

