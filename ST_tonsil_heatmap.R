rm(list=ls())

library(tidyr)
library(dplyr)
library(ggplot2)
library(glue)
library(stringr)
library(ggpubr)
library(ggsci)

#library(rmote)
#start_rmote(port=4339)

library(circlize)
library(ComplexHeatmap)

## we fixed parameters based on benchmarking starling results
## we then run starling on every cells from the 16 tonsil images 5 times
## The paper reports the best starling model via negative log likelihood scores
##    - the script generates 5 cluster centroid plots
## we map clusters to cell types for further analysis

#setwd("~/Desktop")

minmax <- function(x, na.rm = TRUE) {
  return((x- min(x, na.rm = TRUE)) /(max(x, na.rm = TRUE)-min(x, na.rm = TRUE)))
}

ct_labels <- c('Epithelial', 'B', 'Macrophage', 'B', 'CD4 T', 'B', 'CD4 T', #'IFNg+ B'
               'Mix', 'Mix', 'CD4 T', 'CD4 T', 'Stromal', 'Epithelial', 'Monocytes',
               'B', 'Neutrophils', 'Mix', 'B', 'Endothelial', 'B', 'B',
               'B', 'B', 'Unknown', 'Unknown', 'Epithelial')

mat1 <- read.csv("~/Downloads/s250000_KM_n50_t3_r6_nm2_lv1_cs1_rr1_init_c.csv")[,-1]
mat2 <- read.csv("~/Downloads/s250000_KM_n50_t3_r6_nm2_lv1_cs1_rr1_st_c.csv")[,-1]
samp_lb <- read.csv("~/Downloads/s250000_KM_n50_t3_r6_nm2_lv1_cs1_rr1_st_l.csv")[,-1]

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

dmat <- apply(rbind(mat1, mat2[,-25]), 2, minmax)

i=0
ma1 <- t(as.matrix(dmat[(i*nc+1):((i+1)*nc),]))
colnames(ma1) <- ct_labels

f1 <- colorRamp2(seq(min(ma1), max(ma1), length = 3), c("blue", "#EEEEEE", "red"), space = "RGB")
haa1 = column_ha <- HeatmapAnnotation(" " = anno_barplot(as.integer(count[,4]), axis = TRUE))
htt1 <- Heatmap(ma1, 
                column_title = " ",
                col = f1, #rev(rainbow(3)),
                top_annotation = haa1,
                column_split = colnames(ma1),
                cluster_rows = TRUE, cluster_columns = TRUE, show_heatmap_legend = FALSE, border = TRUE)

pdf('tonsil_heatmap_init.pdf')
draw(htt1, column_title = "Cell Type", row_title = "Marker", column_title_side = 'bottom')
dev.off()

i=1
ma2 <- t(as.matrix(dmat[(i*nc+1):((i+1)*nc),]))
colnames(ma2) <- ct_labels

f2 <- colorRamp2(seq(min(ma2), max(ma2), length = 3), c("blue", "#EEEEEE", "red"), space = "RGB")
column_ha <- HeatmapAnnotation(" " = anno_barplot(as.integer(count[,5]), axis = TRUE))
htt2 <- Heatmap(ma2, column_title = " ",
                col = f2, #rev(rainbow(3)),
                top_annotation = column_ha,
                column_split = colnames(ma2),
                cluster_rows = TRUE, cluster_columns = TRUE, show_heatmap_legend = FALSE, border = TRUE)

pdf('tonsil_heatmap_st.pdf')
draw(htt2, column_title = "Cell Type", row_title = "Marker", column_title_side = 'bottom')
dev.off()

pdf('tonsil_heatmap_combined.pdf', width=18, height=6)
draw(htt1 + htt2, column_title = "Cell Type", row_title = "Marker", column_title_side = 'bottom')
dev.off()
