library(nlme)
library(psych)
library(rstatix)
library(tidyr)
library(psych)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggpubr)
library(factoextra)
library(NbClust)
library(gridExtra) 
library(clValid)



# Functions
fviz_nbclust_2 <- function (x, FUNcluster = NULL, method = c("silhouette", "wss", 
                                                           "gap_stat"), diss = NULL, k.max = 10, nboot = 100, verbose = interactive(), 
                          barfill = "#184e77", barcolor = "#184e77", linecolor = "#184e77", 
                          print.summary = TRUE, ...) 
{
  set.seed(123)
  if (k.max < 2) 
    stop("k.max must bet > = 2")
  method = match.arg(method)
  if (!inherits(x, c("data.frame", "matrix")) & !("Best.nc" %in% 
                                                  names(x))) 
    stop("x should be an object of class matrix/data.frame or ", 
         "an object created by the function NbClust() [NbClust package].")
  if (inherits(x, "list") & "Best.nc" %in% names(x)) {
    best_nc <- x$Best.nc
    if (any(class(best_nc) == "numeric") ) 
      print(best_nc)
    else if (any(class(best_nc) == "matrix") )
      .viz_NbClust(x, print.summary, barfill, barcolor)
  }
  else if (is.null(FUNcluster)) 
    stop("The argument FUNcluster is required. ", "Possible values are kmeans, pam, hcut, clara, ...")
  else if (!is.function(FUNcluster)) {
    stop("The argument FUNcluster should be a function. ", 
         "Check if you're not overriding the specified function name somewhere.")
  }
  else if (method %in% c("silhouette", "wss")) {
    if (is.data.frame(x)) 
      x <- as.matrix(x)
    if (is.null(diss)) 
      diss <- stats::dist(x)
    v <- rep(0, k.max)
    if (method == "silhouette") {
      for (i in 2:k.max) {
        clust <- FUNcluster(x, i, ...)
        v[i] <- .get_ave_sil_width(diss, clust$cluster)
      }
    }
    else if (method == "wss") {
      for (i in 1:k.max) {
        clust <- FUNcluster(x, i, ...)
        v[i] <- .get_withinSS(diss, clust$cluster)
      }
    }
    df <- data.frame(clusters = as.factor(1:k.max), y = v, 
                     stringsAsFactors = TRUE)
    ylab <- "Total Within Sum of Square"
    if (method == "silhouette") 
      ylab <- "Average silhouette width"
    p <- ggpubr::ggline(df, x = "clusters", y = "y", group = 1, 
                        color = linecolor, ylab = ylab, xlab = "Number of clusters k", 
                        main = "Optimal number of clusters")
    if (method == "silhouette") 
      p <- p + geom_vline(xintercept = which.max(v), linetype = 2, 
                          color = linecolor)
    return(p)
  }
  else if (method == "gap_stat") {
    extra_args <- list(...)
    gap_stat <- cluster::clusGap(x, FUNcluster, K.max = k.max, 
                                 B = nboot, verbose = verbose, ...)
    if (!is.null(extra_args$maxSE)) 
      maxSE <- extra_args$maxSE
    else maxSE <- list(method = "firstSEmax", SE.factor = 1)
    p <- fviz_gap_stat(gap_stat, linecolor = linecolor, 
                       maxSE = maxSE)
    return(p)
  }
}

.viz_NbClust <- function (x, print.summary = TRUE, barfill = "#184e77", 
                          barcolor = "#184e77") 
{
  best_nc <- x$Best.nc
  if (any(class(best_nc) == "numeric") )
    print(best_nc)
  else if (any(class(best_nc) == "matrix") ) {
    best_nc <- as.data.frame(t(best_nc), stringsAsFactors = TRUE)
    best_nc$Number_clusters <- as.factor(best_nc$Number_clusters)
    if (print.summary) {
      ss <- summary(best_nc$Number_clusters)
      cat("Among all indices: \n===================\n")
      for (i in 1:length(ss)) {
        cat("*", ss[i], "proposed ", names(ss)[i], 
            "as the best number of clusters\n")
      }
      cat("\nConclusion\n=========================\n")
      cat("* According to the majority rule, the best number of clusters is ", 
          names(which.max(ss)), ".\n\n")
    }
    df <- data.frame(Number_clusters = names(ss), freq = ss, 
                     stringsAsFactors = TRUE)
    p <- ggpubr::ggbarplot(df, x = "Number_clusters", 
                           y = "freq", fill = barfill, color = barcolor) + 
      labs(x = "Number of clusters k", y = "Frequency among all indices", 
           title = paste0("Optimal number of clusters - k = ", 
                          names(which.max(ss))))
    return(p)
  }
}
# assign them to the factoextra namespace
environment(fviz_nbclust) <- asNamespace("factoextra")
assignInNamespace("fviz_nbclust",fviz_nbclust,"factoextra")
environment(.viz_NbClust) <- asNamespace("factoextra")
assignInNamespace(".viz_NbClust",.viz_NbClust,"factoextra")


#path="C:/Users/gabot/OneDrive - McGill University/Desktop/Neurosubtyping_code/Datasets/clustering/dimensionality_datasets/listwise_deletion"
path="C:/Users/gabot/OneDrive - McGill University/Desktop/Github_repos/infant_neurosubs/Datasets/hierarchical/source/2025_04_23/combined/MICE/"
# Set directory
setwd(path)

date <- "2025_04_23"


## List all files inside this path 
csv_files <- list.files(path = path, pattern = "\\.csv$", full.names = TRUE)

# Create analysis df 

analysis_clusters <- data.frame(
  id = character(),
  NB_bBest = numeric(),
  clvalid_con=numeric(),
  clvalid_dunn=numeric(),
  clvalid_sil=numeric(),
  
  clvalid_sil_score=numeric(),
  clValid_dunn_score=numeric(),
  clvalid_conn_score=numeric(),
  
  
  nb_sil_score = numeric(),
  nb_kl_score= numeric(),
  nb_ch_score= numeric(),
  nb_hartigan_score=numeric(),
  
  total_measures=numeric(),
  
  stringsAsFactors = FALSE
)





for (file in csv_files) {
  com_df <- read.csv(file)
  print(basename(file))
  cat("Reading:", basename(file), "\n")
  
  # Normalize the data
  com_df <- com_df %>% 
    mutate(across(everything(), ~ (.-min(.))/(max(.)-min(.))))
  
  # Calculate best number of clusters based on NB clust
  
  nb <- NbClust(com_df, distance = "euclidean", min.nc = 2,
                max.nc = 5, method = "ward.D2")
  
  # Extract the NBClust best choice 
  indices<-t(nb$Best.nc)[,1]
  nb_choice <- as.numeric(names(table(indices))[which.max(table(indices))])
  nb_choice
  
  
  # Calculate best number of clusters based on ClValid
  
  clvalid_results <- clValid(com_df, nClust = 2:5, 
                    clMethods = c("hierarchical"), validation = "internal")
  
  #Extract the CLValid best choice

  clvalid_df=optimalScores(clvalid_results)
  clvalid_df  

  # Summary
  #summary(intern)
  
  t(nb$Best.nc)
  

  # Create label for best Siljouette score
  
  sil=paste0("Highest silhouette is:",t(nb$Best.nc)[13,2])
  
  
  pdf_file <- paste0("C:/Users/gabot/OneDrive - McGill University/Desktop/Github_repos/infant_neurosubs/Figures/NBCLUST/NBCLUST_" , basename(file),date, ".pdf")
  
  # Create PDF file
  pdf(file = pdf_file, width = 12, height = 6)  # Adjust width and height as needed
  
  p1 <- fviz_nbclust_2(nb)

  p2 <- # Elbow method
    fviz_nbclust(com_df, kmeans, method = "wss") +
    geom_vline(xintercept = 4, linetype = 2)+
    labs(subtitle = "Elbow method")
  
  p3 <-# Silhouette method
    fviz_nbclust(com_df, kmeans, method = "silhouette")+
    labs(subtitle = "Silhouette method")+
    annotate("text", x = Inf, y = Inf, label = basename(file),
             hjust = 2.9, vjust =1, size = 4, color = "red")+
    annotate("text", x = Inf, y = Inf, label = sil,
             hjust = 1.1, vjust =1, size = 4, color = "black")
  
  # Arrange grid
  grid.arrange(p1, p2, p3, nrow = 1)
  dev.off()
  
  # Save P1 as png
  png_file <- paste0("C:/Users/gabot/OneDrive - McGill University/Desktop/Github_repos/infant_neurosubs/Figures/NBCLUST/NBCLUST_" , basename(file),date, ".png")
  
  ggsave(p1, filename = png_file, width = 6, height = 8) # Adjust width and height as needed
  

  
  new_df = data.frame(
  id = basename(file),
  NB_bBest = nb_choice,
  clvalid_con= clvalid_df[1,3],
  clvalid_dunn= clvalid_df[2,3],
  clvalid_sil= clvalid_df[3,3],
  
  
  clvalid_sil_score=clvalid_df[3,1],
  clValid_dunn_score=clvalid_df[2,1],
  clvalid_conn_score=clvalid_df[1,1],
  
  nb_sil_score =   t(nb$Best.nc)[13],
  nb_kl_score=   t(nb$Best.nc)[1],
  nb_ch_score=   t(nb$Best.nc)[2],
  nb_hartigan_score=  t(nb$Best.nc)[3],
  total_measures= sum(clvalid_df[3,1]+clvalid_df[2,1]+clvalid_df[1,1]))
  
  

  
  analysis_clusters <- rbind(analysis_clusters, new_df)
  
  
  analysis_clusters
  
  # Check elbow method
  
}

write.csv(analysis_clusters,"C:/Users/gabot/OneDrive - McGill University/Desktop/Github_repos/infant_neurosubs/Datasets/NBCLUST/nbclust_results.csv")
nb

