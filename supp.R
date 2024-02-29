rm(list=ls())

library(tidyr)
library(dplyr)
library(ggplot2)
library(glue)
library(stringr)
library(ggpubr)
library(ggsci)

library(circlize)
library(ComplexHeatmap)

setwd("~/Desktop/plot/supplementary/")

## Figure 1
basel <- read.csv("fig1_basel_exp_mat.csv")
meta <- read.csv("fig1_meta_exp_mat.csv")
tonsil <- read.csv("fig1_tonsil_exp_mat.csv")

names(basel) <- c("","sample_id","cell_id","x","y","area","area_convex","has_neighbor","nuclei_count","larger_than30","Argon","RT2","RT3","RT4","RT5","RT6","RT7","RT8","RT9","totHH3","H3K27me3","CK5","Fibronectin","CK19","CK8n18","Twist","CD68","CK14","SMA","Vimentin","cMyc","HER2","CD3","pHH3","ERK","Slug","RabbitIgGHL","PR1","PR2","p53","CD44","EpCAM","CD45","GATA3","CD20","BCatenin","CAIX","ECadherin","Ki67","EGFR","S6","Sox9","vWF","CD31","mTOR","CK7","pCK","kE","cPARP","cCaspase3","DNA1","DNA2")
names(meta) <- c("","sample_id","cell_id","x","y","area","area_convex","has_neighbor","nuclei_count","larger_than30","Argon","totHH3","Xe126","I127","Xe131","Xe134","H3K27me3","Ce140","CK5","Fibronectin","CK19","CK8n18","Twist","CD68","CK14","SMA","Vimentin","cMyc","HER2","CD3","pHH3","ERK","Slug","ER","PR","p53","CD44","EpCAM","CD45","GATA3","CD20","BCatenin","CAIX","ECadherin","Ki67","EGFR","pS6","Sox9","vWFCD31","mTOR","CK7","pCK","cPARP","DNA1","DNA2","Hg202","Pb204","Pb206","Pb207","Pb208")
names(tonsil) <- c("","sample_id","cell_id","x","y","area","area_convex","has_neighbor","nuclei_count","larger_than30","SMA","E_Cadherin","Cytokeratin","HLA_DR","Vimentin","CD28","CD15","CD45RA","CD66b","CD20","CD68","CD4","CD8","CD11c","CD45RO","CD3","IFNg","TCF1","CD14","CD56","PD1","DNA1","DNA2","CD45","PNAd","CD31")

basel[,"cohort"] <- "Jackson 2020"
meta[,"cohort"] <- "Ali 2020"
tonsil[,"cohort"] <- "Tonsil"

df <- rbind(basel[,c("cohort", "CD3", "CD20")], meta[,c("cohort", "CD3", "CD20")], tonsil[,c("cohort", "CD3", "CD20")])

ggplot(df, aes(x=asinh(CD3/5), y=asinh(CD20/5))) + geom_point() + facet_grid(cols = vars(cohort)) + 
  theme(plot.caption = element_text(hjust = 0.5)) +
  theme_pubr() + 
  ylim(c(0, 5)) + 
  labs(caption = "Supplementary Figure 1. Scatterplot of CD3 and CD20 across three IMC datasets.",
       x = "CD3", y = "CD20")
ggsave("fig1.png", width = 5, height = 5)

# Figure 4

rm(list=ls())

minmax <- function(x, na.rm = TRUE) {
  return((x- min(x, na.rm = TRUE)) /(max(x, na.rm = TRUE)-min(x, na.rm = TRUE)))
}

mat1 <- read.csv("fig4_supplement_pg_in_expression.csv")
mat2 <- read.csv("fig4_supplement_fs_in_expression.csv")
mat3 <- read.csv("fig4_supplement_km_in_expression.csv")
mat4 <- read.csv("fig4_supplement_km2_in_expression.csv")
mat5 <- read.csv("fig4_supplement_km2_st_expression.csv")

samp_lb1 <- read.csv("fig4_supplement_pg_in_label.csv")[,-1]
samp_lb2 <- read.csv("fig4_supplement_fs_in_label.csv")[,-1]
samp_lb3 <- read.csv("fig4_supplement_km_in_label.csv")[,-1]
samp_lb4 <- read.csv("fig4_supplement_km2_label.csv")[,-1]

nc <- dim(mat1)[1]
ct_labels <- c('T', 'B', 'Unknown', 'B', 'Unknown', 'T', 'B', 'T', 'B', 'B', 'B', 'Epithelial',
               'B', 'Unknown', 'B', 'B', 'B', 'Endothelial', 'Stromal', 'Macrophage', 'Immune',
               'Endothelial', 'Stromal', 'Neutrophils', 'Epithelial')

ct_labels1 <- 1:nc

count <- NULL
for (c in 0:(nc-1)) {
  count <- rbind(count, c(sum(samp_lb1[,"PG"] == c),
                          sum(samp_lb2[,"FS"] == c),
                          sum(samp_lb3[,"KM"] == c),
                          sum(samp_lb4[,"KM"] == c),
                          sum(samp_lb4[,"ST"] == c)))
}

dmat <- apply(rbind(mat1, mat2, mat3, mat4, mat5[,-26]), 2, minmax)[,-1]

i=0
ma1 <- t(as.matrix(dmat[(i*nc+1):((i+1)*nc),]))
colnames(ma1) <- ct_labels

f1 <- colorRamp2(seq(min(dmat), max(dmat), length = 3), c("blue", "#EEEEEE", "red"), space = "RGB")

haa1 <- HeatmapAnnotation(
  " " = anno_barplot(as.integer(count[,1]), axis = TRUE, ylim = c(0,65000)),
  '  ' = ct_labels,
  col = list('  ' = c(
    'Stromal' = '#dbdb8d', 'B' = '#279e68', 'T' = '#8c564b',
    'Immune' = "#999999", 'Unknown' = '#ffbb78', 'Epithelial' = '#17becf',
    'Macrophage' = '#aec7e8', 'Endothelial' = '#b5bd61', 
    'Monocytes' = '#c5b0d5','Neutrophils' = '#f7b6d2')),
  annotation_name_side = "left",
  annotation_name_rot = 0, 
  annotation_name_gp = gpar(fontsize = 8), 
  annotation_legend_param = list('  ' = list(title = "Cell Type", 
                                             labels_gp = gpar(fontsize = 10))))

htt1 <- Heatmap(ma1, 
                name = 'Scaled\nExpression',
                column_title = "PhenoGraph",
                col = f1, #rev(rainbow(3)),
                top_annotation = haa1,
                column_split = colnames(ma1),
                column_labels = ct_labels1,
                cluster_rows = TRUE, cluster_columns = TRUE, show_heatmap_legend = TRUE, border = TRUE)

pdf('fig4_supplement_pg_in_ct_heatmap.pdf')
draw(htt1, column_title = "Cluster", row_title = "Marker", column_title_side = 'bottom', merge_legend = TRUE)
dev.off()

i=1
ma2 <- t(as.matrix(dmat[(i*nc+1):((i+1)*nc),]))

ct_labels <- c('Epithelial', 'T', 'Endothelial', 'Endothelial', 'Endothelial', 'Endothelial', 'Unknown',
               'Epithelial', 'B', 'Unknown', 'Unknown', 'T', 'Unknown', 'Neutrophils', 'B', 'Epithelial',
               'Epithelial', 'Endothelial', 'T', 'B', 'Epithelial', 'Epithelial', 'Unknown', 'Macrophage',
               'Unknown')

colnames(ma2) <- ct_labels

haa2 <- HeatmapAnnotation(
  " " = anno_barplot(as.integer(count[,2]), axis = FALSE, ylim = c(0,65000)),
  '  ' = ct_labels,
  col = list('  ' = c(
    'Stromal' = '#dbdb8d', 'B' = '#279e68', 'T' = '#8c564b',
    'Immune' = "#999999", 'Unknown' = '#ffbb78', 'Epithelial' = '#17becf',
    'Macrophage' = '#aec7e8', 'Endothelial' = '#b5bd61', 
    'Monocytes' = '#c5b0d5','Neutrophils' = '#f7b6d2')),
  annotation_name_side = "left",
  annotation_name_rot = 0, 
  annotation_name_gp = gpar(fontsize = 8), 
  annotation_legend_param = list('  ' = list(title = "Cell Type", 
                                             labels_gp = gpar(fontsize = 10))))

htt2 <- Heatmap(ma2, 
                name = 'Scaled\nExpression',
                column_title = "FlowSOM",
                col = f1, #rev(rainbow(3)),
                top_annotation = haa2,
                column_split = colnames(ma2),
                column_labels = ct_labels1,
                cluster_rows = TRUE, cluster_columns = TRUE, show_heatmap_legend = TRUE, border = TRUE)

pdf('fig4_supplement_fs_in_ct_heatmap.pdf')
draw(htt2, column_title = "Cluster", row_title = "Marker", column_title_side = 'bottom', merge_legend = TRUE)
dev.off()

i=3
ma4 <- t(as.matrix(dmat[(i*nc+1):((i+1)*nc),]))
ct_labels <- c('Endothelial', 'B', 'T', 'Immune', 'B', 'B', 'T', 'B', 'Unknown', 'Stromal', 'Epithelial',
               'B', 'T', 'Macrophage', 'Endothelial', 'Stromal', 'B', 'Monocytes', 'B', 'Monocytes',
               'Neutrophils', 'B', 'T', 'B', 'Epithelial')

colnames(ma4) <- ct_labels

haa4 <- HeatmapAnnotation(
  " " = anno_barplot(as.integer(count[,4]), axis = FALSE, ylim = c(0,65000), axis_param = list(side = 'right')),
  '  ' = ct_labels,
  col = list('  ' = c(
    'Stromal' = '#dbdb8d', 'B' = '#279e68', 'T' = '#8c564b',
    'Immune' = "#999999", 'Unknown' = '#ffbb78', 'Epithelial' = '#17becf',
    'Macrophage' = '#aec7e8', 'Endothelial' = '#b5bd61', 
    'Monocytes' = '#c5b0d5','Neutrophils' = '#f7b6d2')),
  #annotation_name_side = "left",
  annotation_name_rot = 0, 
  annotation_name_gp = gpar(fontsize = 8), 
  annotation_legend_param = list('  ' = list(title = "Cell Type", 
                                             labels_gp = gpar(fontsize = 10))))

htt4 <- Heatmap(ma4, 
                name = 'Scaled\nExpression',
                column_title = "Kmeans",
                col = f1, #rev(rainbow(3)),
                top_annotation = haa4,
                column_split = colnames(ma4),
                column_labels = ct_labels1,
                cluster_rows = TRUE, cluster_columns = TRUE, show_heatmap_legend = TRUE, border = TRUE)

pdf('fig4_supplement_km_in_ct_heatmap.pdf')
draw(htt4, column_title = "Cluster", row_title = "Marker", column_title_side = 'bottom', merge_legend = TRUE)
dev.off()
