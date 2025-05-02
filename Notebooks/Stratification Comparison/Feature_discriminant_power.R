# Validation of which features are most informative 

# ============================================================================ #
# 1. SETUP AND DATA PREPARATION
# ============================================================================ #

# Load required packages
library(tidyverse)       # For data manipulation and visualization
library(clusterpval)     # For computing p-values for cluster differences
library(fastcluster)     # For efficient hierarchical clustering
library(RColorBrewer)    # For color palettes
library(ggdendro)        # For better dendrogram visualization
library(cowplot)    # For plot_grid()
library(webshot)    # For converting HTML to image
library(magick)     # For image handling
library(dendextend)  # For dendrogram manipulatio

# Set custom color palette (using your preferred colors)
hc_palette <- c("#184e77", "#7570b3", "#1b9e77")  # Blue, purple, green
options(ggplot2.discrete.colour = hc_palette)

# ============================================================================ #
# 2. DATA LOADING AND PREPROCESSING
# ============================================================================ #
# Select date
date= "2025_04_23"

# 
# Load the MICE imputed data
file_name= paste0("C:/Users/gabot/OneDrive - McGill University/Desktop/Github_repos/infant_neurosubs/Datasets/hierarchical/source/", 
                  date, "/combined/MICE/5Dcombined_data_", date,".csv")
cluster_data <- read.csv(file_name,  header = TRUE, 
  sep = ","
)

# Normalize the relevant features
cluster_data_normed <- scale(cluster_data[, c(
  "front_gamma_6",
  "auditory_con_6", 
  "lang_comp_con_6", 
  "speech_con_6_left", 
  "gamma_lat_6"
)])

# Create feature matrix
X <- as.matrix(cluster_data_normed)

# ============================================================================ #
# 3. HIERARCHICAL CLUSTERING
# ============================================================================ #

# Compute squared Euclidean distance matrix
dist_matrix <- dist(X, method = "euclidean")^2

# Perform hierarchical clustering (ward.D2 matches scipy's 'ward')
hcl <- hclust(dist_matrix, method = "ward.D2")

# ============================================================================ #
# 4. VISUALIZATION
# ============================================================================ #

# Basic dendrogram plot
plot(as.dendrogram(hcl), 
     leaflab = "none",
     main = "Hierarchical Clustering Dendrogram",
     xlab = "Observations",
     ylab = "Height")

# Add cut line for 3 clusters
abline(h = (hcl$height[nrow(X) - 3] + 0.1),  # Adjusted for better visibility
       lty = "dashed", 
       col = "darkgrey")

# Highlight clusters
rect_hier_clusters(hcl, 
                   k = 3, 
                   which = 1:3, 
                   border = hc_palette)

# ============================================================================ #
# 5. CLUSTER SIGNIFICANCE TESTING
# ============================================================================ #
# Initialize results dataframe
results <- data.frame(
  measure = colnames(X),
  statistic = numeric(ncol(X)),
  p.value = numeric(ncol(X)),
  corrected.p = numeric(ncol(X))
)

# Test if clusters 1 and 2 are significantly different

test_result_1_v_2 <- test_hier_clusters_exact(
  X,
  link = "ward.D",  # Required for exact test
  K = 3,
  k1 = 1,
  k2 = 2,
  hcl = hcl
)

# Print the test result
print(test_result_1_v_2)

# Test to see if clusters 1 and 3 are significantly different
test_result_1_v_3 <- test_hier_clusters_exact(
  X,
  link = "ward.D",  # Required for exact test
  K = 3,
  k1 = 1,
  k2 = 3,
  hcl = hcl
)

print(test_result_1_v_3)


# Test to see if clusters 2 and 3 are significantly different

test_result_3_v_2 <- test_hier_clusters_exact(
  X,
  link = "ward.D",  # Required for exact test
  K = 3,
  k1 = 2,
  k2 = 3,
  hcl = hcl
)
print (test_result_3_v_2)



# ============================================================================ #
# 6. ENHANCED VISUALIZATION 
# ============================================================================ #

  dendogram_file_name <-  paste0("C:/Users/gabot/OneDrive - McGill University/Desktop/Github_repos/infant_neurosubs/Figures/Supplementary/",
  date, "/Hierarchical_dendogram_" , date, ".png")


    # Open PNG device
  png(dendogram_file_name, 
      width = 8,         # 10 inches wide
      height = 6,         # 7 inches tall
      units = "in", 
      res = 600,          # High resolution (600 DPI)
      type = "cairo")     # Better anti-aliasing
  
  # Set graphical parameters
  op <- par(mar = c(5, 5, 4, 2) + 0.1,  # Larger margins
            lwd = 1.5,
            cex.main = 1.5,
            cex.lab = 1.3)
  
  # Create and customize dendrogram
  dend <- as.dendrogram(hcl) %>%
    set("branches_k_color", k = 3, value = hc_palette) %>%
    set("branches_lwd", 2) %>%
    set("labels_col", "transparent")
  
  # Plot dendrogram
  plot(dend,
       xlab = "Participants",
       ylab = "Dissimilarity (Height)",
       ylim = c(0, max(hcl$height) * 1.1),
       axes = FALSE)
  
  # Add custom axis
  axis(2, at = pretty(hcl$height), lwd = 2, cex.axis = 1.2, font.axis = 2)
  
  # Add cut line
  abline(h = (hcl$height[nrow(X) - 3] + 0.1), 
         lty = "dashed", col = "darkgrey", lwd = 3)
  

  # Close device
  dev.off()
  
    # Reset parameters
  par(op)
  
  
  
 # -----------------------------------------------------
  # Panel A: Dendrogram 
  # -----------------------------------------------------

## EEG Class difference 
class_diff_file_name <- paste0("C:/Users/gabot/OneDrive - McGill University/Desktop/Github_repos/infant_neurosubs/Figures/Supplementary/",
               date, "/EEG_Class_Differentiation_Power.png")
 
p1 <- ggdraw() + 
  draw_image(class_diff_file_name) #+
  labs(title = "Cluster Differentiation Power") +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5))

p1
# -----------------------------------------------------
# Panel B: Load from PNG file
# -----------------------------------------------------
p2_file<-dendogram_file_name
 
p2 <- ggdraw() + 
  draw_image(p2_file)  + labs(title = "Hierarchical Clustering Dendogram") +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.4))

p2
# -----------------------------------------------------
# Panel C: Sankey Diagram (from HTML)
# -----------------------------------------------------


sankey_png <- paste0("C:/Users/gabot/OneDrive - McGill University/Desktop/Github_repos/infant_neurosubs/Figures/Supplementary/",
               date, "/sankey_diagram_" , date, ".png")

p3 <- ggdraw() + 
  draw_image(sankey_png) +
  labs(title = "Sankey Diagram of Overlap") +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.4))

p3
# -----------------------------------------------------
# Combine all panels horizontally
# -----------------------------------------------------
final_plot <- plot_grid(
  p2, p3,
  nrow = 1,
  labels = "AUTO")

final_plot
# Save output
output_file <- paste0("C:/Users/gabot/OneDrive - McGill University/Desktop/Github_repos/infant_neurosubs/Figures/main/",
        "Figure_Sankey.png")
ggsave(output_file, final_plot, width = 10, height = 4, dpi = 300)


