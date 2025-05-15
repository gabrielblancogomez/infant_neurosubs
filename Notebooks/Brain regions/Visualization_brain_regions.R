# Install and load required packages
if (!require("ggseg")) install.packages("ggseg")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("tibble")) install.packages("tibble")

library(ggseg)
library(ggplot2)
library(dplyr)
library(tibble)
# Define the color palette

brain_labels(dk)

color_palette <- c('#26c6da', '#1a759f', '#0a3d62')
ggplot() +
  geom_brain(atlas = dk, side = "lateral")


 # Create a function to plot brain regions
plot_brain_regions <- function(feature_name, regions_df) {
  # Plot
  p <- ggplot(regions_df) +
    geom_brain(atlas = dk,
               side = "lateral",
               aes(fill = p)) +
    labs(title = feature_name, fontsize=30) +
    scale_fill_gradient(guide = "none") +  # This removes the legend
    
    # Remove grid lines
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    
    # Remove axis labels
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    
    # Remove axis ticks
    theme(axis.ticks = element_blank()) +
    # Remove axis text
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank()) +
    
    # increase the font size of the title 
    theme(text = element_text(size = 20)) +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(p)
}

# Define regions for each feature (without hemisphere suffixes)
frontal_gamma_regions <- c("rostral middle frontal", "superior frontal")
gamma_lateralization_regions <- c("pars opercularis", "pars triangularis",
                                  "pars orbitalis", "rostral middle frontal",
                                  "lateral orbitofrontal", "caudal middle frontal",
                                  "frontal pole")
lang_comp_regions <- c("middle temporal", "inferior parietal",
                       "supramarginal", "pars opercularis",
                       "pars triangularis", "pars orbitalis")
speech_regions <- c("pars opercularis", "pars triangularis",
                    "pars orbitalis", "precentral",
                    "superior temporal")
auditory_regions <- c("superior temporal", "transverse temporal")

# Create data frames for each feature with proper hemi specification
frontal_gamma_data <- tibble(
  region = rep(frontal_gamma_regions, each = 2),
  p = 0.5,
  hemi = rep(c("left", "right"), times = length(frontal_gamma_regions)),
  feature = "Frontal Gamma Power", 
  color_palette= "#0a3d62"
)

gamma_lateralization_data <- tibble(
  region = rep(gamma_lateralization_regions, each = 2),
  p = 0.5,
  hemi = rep(c("left", "right"), times = length(gamma_lateralization_regions)),
  feature = "Gamma Lateralization"
)

lang_comp_data <- tibble(
  region = rep(lang_comp_regions, each = 2),
  p = 0.5,
  hemi = rep(c("left", "right"), times = length(lang_comp_regions)),
  feature = "Connectivity Language"
)

# Speech regions are left-hemisphere only
speech_data <- tibble(
  region = speech_regions,
  p = 0.5,
  hemi = "left",
  feature = "Connectivity Speech"
)

auditory_data <- tibble(
  region = rep(auditory_regions, each = 2),
  p = 0.5,
  hemi = rep(c("left", "right"), times = length(auditory_regions)),
  feature = "Connectivity Auditory"
)


# Create plots for each feature
frontal_gamma_plot <- plot_brain_regions("Frontal Gamma Power", frontal_gamma_data)
frontal_gamma_plot
gamma_lateralization_plot <- plot_brain_regions("Frontal Gamma Lateralization", gamma_lateralization_data)
lang_comp_plot <- plot_brain_regions("Connectivity Language Network", lang_comp_data)
speech_plot <- plot_brain_regions("Connectivity Speech Network", speech_data)
auditory_plot <- plot_brain_regions("Connectivity Auditory Network", auditory_data)
speech_plot

library(patchwork)

combined_plot <- (
  (frontal_gamma_plot + theme(plot.margin = margin(0,0,0,0))) | 
    (gamma_lateralization_plot + theme(plot.margin = margin(0,0,0,0)))
) / (
  (lang_comp_plot + theme(plot.margin = margin(0,0,0,0))) | 
    (speech_plot + theme(plot.margin = margin(0,0,0,0)))
) / (
  (auditory_plot + theme(plot.margin = margin(0,0,0,0))) | 
    plot_spacer()
) +
  plot_layout(heights = c(1, 1, 0.8)) &
  theme(plot.margin = margin(0,0,0,0))

# Display the combined plot
combined_plot

# Set current directory to where this code file lives
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Save the combined plot
ggsave("../../Figures/Supplementary/2025_04_23/brain_regions_combined.png", combined_plot, width = 12, height = 9, dpi = 300)
ggsave("../Visuals/dashboard_source/brain_regions_combined.png", combined_plot, width = 12, height = 9, dpi = 300)

# Alternatively, save each plot separately
ggsave("../../Figures/Supplementary/2025_04_23/frontal_gamma_power.png", frontal_gamma_plot, width = 6, height = 5, dpi = 300)
ggsave("../../Figures/Supplementary/2025_04_23/gamma_lateralization.png", gamma_lateralization_plot, width = 6, height = 5, dpi = 300)
ggsave("../../Figures/Supplementary/2025_04_23/connectivity_language.png", lang_comp_plot, width = 6, height = 5, dpi = 300)
ggsave("../../Figures/Supplementary/2025_04_23/connectivity_speech.png", speech_plot, width = 6, height = 5, dpi = 300)
ggsave("../../Figures/Supplementary/2025_04_23/connectivity_auditory.png", auditory_plot, width = 6, height = 5, dpi = 300)

