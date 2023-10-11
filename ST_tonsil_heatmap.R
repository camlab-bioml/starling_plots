rm(list=ls())

library(tidyr)
library(dplyr)
library(ggplot2)
library(glue)
library(stringr)
library(ggpubr)
library(ggsci)

library(rmote)
start_rmote(port=4339)

library(circlize)
library(ComplexHeatmap)

## we fixed parameters based on benchmarking starling results
## we then run starling on every cells from the 16 tonsil images 5 times
## The paper reports the best starling model via negative log likelihood scores
##    - the script generates 5 cluster centroid plots
## we map clusters to cell types for further analysis

## minmax a vector
minmax <- function(x, na.rm = TRUE) {
  return((x- min(x, na.rm = TRUE)) /(max(x, na.rm = TRUE)-min(x, na.rm = TRUE)))
}

col_fun <- colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
lgd <- Legend(col_fun = col_fun, title = "Exp.")

setwd("/home/campbell/yulee/project/result/tonsil/_2bm/dc/")
ct_labels <- c('Epithelial', 'B', 'Macrophage', 'B', 'CD4 T', 'B', 'CD4 T', #'IFNg+ B'
               'Mix', 'Mix', 'CD4 T', 'CD4 T', 'Stromal', 'Epithelial', 'Monocytes',
               'B', 'Neutrophils', 'Mix', 'B', 'Endothelial', 'B', 'B',
               'B', 'B', 'Unknown', 'Unknown', 'Epithelial')

file_names <- list.files()
f1 <- file_names[grep("250000", file_names)]
f2 <- f1[grep("n50", f1)]
f3 <- f2[grep("t3", f2)]
f4 <- f3[grep("nm2", f3)]
f5 <- f4[grep("lv1_", f4)]
f6 <- f5[grep("cs1", f5)]
f7 <- f6[grep("rr1", f6)]
f8 <- f7[grep("init_c.csv", f7)]
f9 <- f8[grep("KM", f8)]

for ( ii in 1:length(f9) ) {
  
  mat1 <- read.csv(f9[ii])[,-1]
  mat2 <- read.csv(paste0(strsplit(f9[ii], "_init_c.csv")[[1]], "_st_c.csv"))[,-1]
  
  samp_lb <- read.csv(paste0(strsplit(f9[ii], "_init_c.csv")[[1]], "_st_l.csv"))[,-1]
  
  nc <- dim(mat1)[1]
  
  count <- NULL
  for (c in 0:(nc-1)) {
    
    not_in_n_are_in  <- samp_lb[,"INIT"] != c & samp_lb[,"ST"] == c ## werent in & are in
    were_in_n_are_in <- samp_lb[,"INIT"] == c & samp_lb[,"ST"] == c ## were in & are in
    were_in_n_not_in <- samp_lb[,"INIT"] == c & samp_lb[,"ST"] != c ## were in & not in
    
    count <- rbind(count, c(sum(not_in_n_are_in), sum(were_in_n_are_in), sum(were_in_n_not_in), 
                            sum(samp_lb[,"INIT"] == c),
                            sum(samp_lb[,"ST"] == c)))
  }
  
  print(ii)
  dmat <- apply(rbind(mat1, mat2[,-25]), 2, minmax)
  
  i=0
  ma1 <- as.matrix(dmat[(i*nc+1):((i+1)*nc),])
  rownames(ma1) <- ct_labels
  haa1 = rowAnnotation(" " = anno_barplot(as.integer(count[,4]), axis = TRUE))
  htt1 <- Heatmap(ma1, right_annotation = haa1, column_title = "Initialization (Kmeans)", column_title_side = 'top', #column_names_rot = 0,
                  cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = FALSE)

  i=1
  ma2 <- as.matrix(dmat[(i*nc+1):((i+1)*nc),])
  rownames(ma2) <- ct_labels
  haa2 = rowAnnotation(" " = anno_barplot(as.integer(count[,5]), axis = TRUE))
  htt2 <- Heatmap(ma2, right_annotation = haa2, column_title = "STARLING", column_title_side = 'top',
                  cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = FALSE)

  pdf(paste0('/home/campbell/yulee/project/result/plots/', strsplit(f9[ii], "_init_c.csv")[[1]], "_tonsil_heatmap_init.pdf"))
  draw(htt1, annotation_legend_list = lgd, row_title = "Cluster", column_title = "Marker", column_title_side = 'bottom')
  dev.off()
  
  pdf(paste0('/home/campbell/yulee/project/result/plots/', strsplit(f9[ii], "_init_c.csv")[[1]], "_tonsil_heatmap_st.pdf"))
  draw(htt2, annotation_legend_list = lgd, row_title = "Cluster", column_title = "Marker", column_title_side = 'bottom')
  dev.off()
  
  pdf(paste0('/home/campbell/yulee/project/result/plots/', strsplit(f9[ii], "_init_c.csv")[[1]], "_tonsil_heatmap_combined.pdf"), width=18, height=6)
  draw(htt1 + htt2, annotation_legend_list = lgd, row_title = "Cluster", column_title = "Marker", column_title_side = 'bottom')
  dev.off()
}



