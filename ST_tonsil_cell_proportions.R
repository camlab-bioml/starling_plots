rm(list=ls())

library(tidyr)
library(dplyr)
library(ggplot2)
library(glue)
library(stringr)
library(ggpubr)
library(ggsci)

library(circlize)
library(hrbrthemes)
library(ComplexHeatmap)

## assignment/cell type proportions of each image (ROI) and patient

samp_lb <- as.data.frame(read.csv("/home/campbell/yulee/project/result/tonsil/_2bm/dc/s250000_KM_n50_t3_r6_nm2_lv1_cs1_rr1_st_l.csv"))
ct_labels <- c('Epithelial', 'B', 'Macrophage', 'B', 'CD4+ T', 'B', 'CD4+ T', #'IFNg+ B'
               'Mix', 'Mix', 'CD4+ T', 'CD4+ T', 'Stromal', 'Epithelial', 'Monocytes',
               'B', 'Neutrophils', 'Mix', 'B', 'Endothelial', 'B', 'B',
               'B', 'B', 'Unknown', 'Unknown', 'Epithelial')

tmp <- NULL
for (ROI in 4:19) {
  
  #i = 0
  select_image = filter(samp_lb, sample == ROI)
  
  if (ROI < 8) {
    patient <- 1
  } else if (ROI < 16) {
    patient <- 2
  } else {
    patient <- 3
  }
  
  for (c in 0:25) {
    
    #ii = 1
    tmp <- rbind(tmp, c(ROI, patient, ct_labels[c+1], 
                        c,
                        dim(filter(select_image, ST == c))[1],
                        dim(select_image)[1]-sum(select_image[,'ST'] == 'doublet')))
  }
}

colnames(tmp) <- c("ROI", 'P', "CT", 'C', 'Cell', 'Total')
res <- filter(as.data.frame(tmp), as.integer(ROI) > 3)
res[,"ROI"] <- as.integer(res[,"ROI"]) - 3
res[,"Cell"] <- as.numeric(res[,"Cell"])
res[,"Total"] <- as.numeric(res[,"Total"])
res[,"ROI_CP"] <- as.numeric(res[,"Cell"]) / as.numeric(res[,"Total"])

tmp1 <- spread(as.data.frame(res[,c("ROI", 'C', 'P', 'ROI_CP')]), C, ROI_CP)
roi <- tmp1[,"ROI"]
p <- tmp1[,'P'] #as.integer(tmp1[,'P'])
tmp1 <- apply(tmp1[,-c(1,2)], 2, as.numeric)
rownames(tmp1) <- roi

###
###

ha1 <- HeatmapAnnotation(Patient = as.character(p), 
                         col = list(Patient = c("1" = "#E69F00", "2" = "#0072B2", "3" = "#D55E00")),
                         annotation_name_side = "left", annotation_name_gp= gpar(fontsize = 8),
                         annotation_legend_param = list(labels_gp = gpar(fontsize = 8)))

row_ha <- rowAnnotation('Cell\nType' = ct_labels, 
                        col = list('Cell\nType'=c('B' = '#279e68', 'CD4+ T' = '#8c564b', "Endothelial" = '#b5bd61', 
                                                  "Epithelial" = '#17becf', "Macrophage" = '#aec7e8', "Mix" = '#ffbb78',
                                                  "Monocytes" = '#c5b0d5', "Neutrophils" = '#f7b6d2', "Stromal" = '#dbdb8d',
                                                  'Unknown' = '#8c6d31')),
                        annotation_name_rot = 0, annotation_name_gp = gpar(fontsize = 8), 
                        annotation_legend_param = list('Cell\nType' = list(title = "Cell Type", 
                                                                           labels_gp = gpar(fontsize = 8))))

col_fun <- colorRamp2(c(0, 1), c("white", "blue"))
htt1 <- Heatmap(t(tmp1), col = col_fun, top_annotation = ha1, left_annotation = row_ha,
                cluster_rows = TRUE, cluster_columns = TRUE, show_heatmap_legend = FALSE,
                row_names_gp = gpar(fontsize = 8),
                column_names_gp = gpar(fontsize = 8))

lgd <-  Legend(col_fun = col_fun, title = "Proportion", labels_gp = gpar(fontsize = 8))
pdf("/home/campbell/yulee/project/result/plots/tonsil_roi_prop_byLabel.pdf", width=7, height=7)
draw(htt1, annotation_legend_list = lgd, row_title = "Cluster", column_title = "ROI", 
     row_title_gp = gpar(fontsize = 8), column_title_gp = gpar(fontsize = 8),
     column_title_side = 'bottom')
dev.off()

###
###

tmp2 <- res %>% group_by(ROI, P, CT) %>% summarize_at(c('Cell', 'Total'), c(sum, mean)) %>% 
  mutate(ROI_CP=Cell_fn1/Total_fn2) %>% as.data.frame() %>% select(ROI, CT, P, ROI_CP) %>%
  spread(CT, ROI_CP)

roi <- tmp2[,"ROI"]
p <- tmp2[,'P'] #as.integer(tmp1[,'P'])
tmp2 <- apply(tmp2[,-c(1,2)], 2, as.numeric)
rownames(tmp2) <- roi

ha1 <- HeatmapAnnotation(Patient = as.character(p), 
                         col = list(Patient = c("1" = "#E69F00", "2" = "#0072B2", "3" = "#D55E00")),
                         annotation_name_side = "left", annotation_name_gp= gpar(fontsize = 8),
                         annotation_legend_param = list(labels_gp = gpar(fontsize = 8)))

row_ha <- rowAnnotation('Cell\nType' = colnames(tmp2), 
                        col = list('Cell\nType'=c('B' = '#279e68', 'CD4+ T' = '#8c564b', "Endothelial" = '#b5bd61', 
                                                  "Epithelial" = '#17becf', "Macrophage" = '#aec7e8', "Mix" = '#ffbb78',
                                                  "Monocytes" = '#c5b0d5', "Neutrophils" = '#f7b6d2', "Stromal" = '#dbdb8d',
                                                  'Unknown' = '#8c6d31')),
                        annotation_name_rot = 0,  annotation_name_gp = gpar(fontsize = 8), 
                        annotation_legend_param = list('Cell\nType' = list(title = "Cell Type", 
                                                                           labels_gp = gpar(fontsize = 8))))

htt1 <- Heatmap(t(tmp2), col = col_fun, top_annotation = ha1, left_annotation = row_ha,
                cluster_rows = TRUE, cluster_columns = TRUE, show_heatmap_legend = FALSE,
                row_names_gp = gpar(fontsize = 8),
                column_names_gp = gpar(fontsize = 8))

pdf("/home/campbell/yulee/project/result/plots/tonsil_roi_prop_byCelltype.pdf", width=7, height=3)
draw(htt1, annotation_legend_list = lgd, row_title = "Cluster", column_title = "ROI", 
     row_title_gp = gpar(fontsize = 8), column_title_gp = gpar(fontsize = 8),
     column_title_side = 'bottom')
dev.off()




