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

## enrichment analysis by ROIs via squidpy package
## two approaches were conducted:
##  - 1. compute p(z_n = c | d_n = 0), then define assignment and then enrichment analysis
##  - 2. compute doublets, and then conduct enrichment analysis
## this script takes the lower triangular and diagonal of the interaction matrix
## the results are stored in *_fs.csv

#counts <- read.csv('/home/campbell/yulee/project/result/plots/m1_fs.csv')
#counts <- read.csv('/home/campbell/yulee/project/result/plots/m2_fs.csv')

p <- c(1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3)

nhood_csd12 <- NULL
nhood_csn12 <- NULL
for ( ii in 1:2 ) {
  
  counts <- read.csv(paste0('/home/campbell/yulee/project/result/plots/m', ii, '_fs.csv'))
  
  nhood_cs <- c()
  for ( i in 4:19 ) {
    tmp = filter(counts, Sample == i)[,-1]
    nhood_cs <- rbind(nhood_cs, tmp[lower.tri(tmp, diag = TRUE)])
  }
  
  colnames(tmp)[2] <- "CD4+ T"
  rownames(tmp) <- colnames(tmp)
  colnames(nhood_cs) <- t(sapply(1:length(colnames(tmp)), function(ii) paste0(rownames(tmp)[ii], "-", colnames(tmp))))[lower.tri(tmp, diag = TRUE)]
  
  nhood_cs <- as.data.frame(nhood_cs) %>% select(-contains('doublet') & -contains('Mix') & -contains('Unknown') & -contains('IFNg+ B'))
  
  nhood_cs[nhood_cs == 'NaN'] <- 0
  rownames(nhood_cs) <- paste0("ROI", 1:16)
  
  nhood_csd <- select(nhood_cs, 'B-B', "CD4+ T-CD4+ T", "Endothelial-Endothelial", "Epithelial-Epithelial",
                      "Macrophage-Macrophage", "Monocytes-Monocytes", "Neutrophils-Neutrophils", "Stromal-Stromal")
  
  nhood_csn <- select(nhood_cs, -c('B-B', "CD4+ T-CD4+ T", "Endothelial-Endothelial", "Epithelial-Epithelial",
                                   "Macrophage-Macrophage", "Monocytes-Monocytes", "Neutrophils-Neutrophils", "Stromal-Stromal"))
  
  colnames(nhood_csd) <- do.call(rbind, strsplit(colnames(nhood_csd), "-"))[,1]
  
  col_fun = colorRamp2(c(0, max(nhood_csd, na.rm = TRUE)), c("white", "blue"))
  lgd = Legend(col_fun = col_fun, title = "Count", labels_gp = gpar(fontsize = 8))
  
  ha1 <- HeatmapAnnotation(Patient = as.character(p), 
                           col = list(Patient = c("1" = "#E69F00", "2" = "#0072B2", "3" = "#D55E00")),
                           annotation_name_side = "left", annotation_name_gp= gpar(fontsize = 8),
                           annotation_legend_param = list(labels_gp = gpar(fontsize = 8)))
  
  row_ha <- rowAnnotation('Cell\nType' = as.character(colnames(nhood_csd)), 
                          col = list('Cell\nType'=c('B' = '#279e68', 'CD4+ T' = '#8c564b', "Endothelial" = '#b5bd61', 
                                                    "Epithelial" = '#17becf', "Macrophage" = '#aec7e8', "Mix" = '#ffbb78',
                                                    "Monocytes" = '#c5b0d5', "Neutrophils" = '#f7b6d2', "Stromal" = '#dbdb8d',
                                                    'Unknown' = '#8c6d31')),
                          annotation_name_rot = 0, annotation_name_gp = gpar(fontsize = 8), 
                          annotation_legend_param = list('Cell\nType' = list(title = "Cell Type", 
                                                                             labels_gp = gpar(fontsize = 8))))
  
  htt1 <- Heatmap(t(nhood_csd), col = col_fun, top_annotation = ha1, 
                  left_annotation = row_ha, #show_row_name = FALSE, 
                  cluster_rows = TRUE, cluster_columns = TRUE, show_heatmap_legend = FALSE,
                  row_names_gp = gpar(fontsize = 8),
                  column_names_gp = gpar(fontsize = 8)) #, column_names_rot = 0
  
  pdf(paste0("/home/campbell/yulee/project/result/plots/tonsil_nhood_m", ii, "_countd.pdf"), width=7, height=3)
  draw(htt1, annotation_legend_list = lgd, 
       row_title = "Interaction (Count)", column_title = "ROI",
       row_title_gp = gpar(fontsize = 8), column_title_gp = gpar(fontsize = 8),
       column_title_side = 'bottom')
  dev.off()
  
  ###
  ###
  
  col_fun = colorRamp2(c(0, max(nhood_csn, na.rm = TRUE)), c("white", "blue"))
  lgd = Legend(col_fun = col_fun, title = "Count", labels_gp = gpar(fontsize = 8))
  
  ha1 <- HeatmapAnnotation(Patient = as.character(p), 
                           col = list(Patient = c("1" = "#E69F00", "2" = "#0072B2", "3" = "#D55E00")),
                           annotation_name_side = "left", annotation_name_gp= gpar(fontsize = 8),
                           annotation_legend_param = list(labels_gp = gpar(fontsize = 8)))
  
  row_ha <- rowAnnotation('Type 1' = as.character(do.call(rbind, strsplit(colnames(nhood_csn), "-"))[,1]),
                          'Type 2' = as.character(do.call(rbind, strsplit(colnames(nhood_csn), "-"))[,2]),
                          col = list('Type 1'= c('B' = '#279e68', 'CD4+ T' = '#8c564b', "Endothelial" = '#b5bd61', 
                                                 "Epithelial" = '#17becf', "Macrophage" = '#aec7e8', "Mix" = '#ffbb78',
                                                 "Monocytes" = '#c5b0d5', "Neutrophils" = '#f7b6d2', "Stromal" = '#dbdb8d',
                                                 'Unknown' = '#8c6d31'),
                                     'Type 2'= c('B' = '#279e68', 'CD4+ T' = '#8c564b', "Endothelial" = '#b5bd61', 
                                                 "Epithelial" = '#17becf', "Macrophage" = '#aec7e8', "Mix" = '#ffbb78',
                                                 "Monocytes" = '#c5b0d5', "Neutrophils" = '#f7b6d2', "Stromal" = '#dbdb8d',
                                                 'Unknown' = '#8c6d31')),
                          annotation_name_rot = 90, annotation_name_gp = gpar(fontsize = 8), 
                          annotation_legend_param = list('Type 1' = list(title = "Cell Type 1", 
                                                                         labels_gp = gpar(fontsize = 8)),
                                                         'Type 2' = list(title = "Cell Type 2", 
                                                                         labels_gp = gpar(fontsize = 8))))
  
  htt1 <- Heatmap(t(nhood_csn), col = col_fun, top_annotation=ha1, 
                  left_annotation = row_ha, #show_row_name = FALSE, 
                  #right_annotation = row_ha2,
                  cluster_rows = TRUE, cluster_columns = TRUE, show_heatmap_legend = FALSE,
                  row_names_gp = gpar(fontsize = 8),
                  column_names_gp = gpar(fontsize = 8)) #, column_names_rot = 0)
  
  pdf(paste0("/home/campbell/yulee/project/result/plots/tonsil_nhood_m", ii, "_countn.pdf"), width=10, height=9)
  draw(htt1, annotation_legend_list = lgd, row_title = "Interaction (Count)", column_title = "ROI", 
       row_title_gp = gpar(fontsize = 8), column_title_gp = gpar(fontsize = 8),
       column_title_side = 'bottom')
  dev.off()
  
  nhood_csd12 <- cbind(nhood_csd12, unlist(nhood_csd))
  nhood_csn12 <- cbind(nhood_csn12, unlist(nhood_csn))
}

colnames(nhood_csd12) <- c("method1", "method2")
colnames(nhood_csn12) <- c("method1", "method2")

ggplot(as.data.frame(nhood_csd12), aes(x=method1, y=method2)) + geom_point() + 
  ggtitle("Interaction frequency: Same cell type pairs") +
  xlab(expression('p('~z[n]~'= c|'~d[n]~'= 0)')) + ylab("Doublet")
ggsave("/home/campbell/yulee/project/result/plots/tonsil_nhood_m12_diagonal_freq_comparisons.pdf")

ggplot(as.data.frame(nhood_csn12), aes(x=method1, y=method2)) + geom_point() + 
  ggtitle("Interaction frequency: Different cell type pairs") +
  xlab(expression('p('~z[n]~'= c|'~d[n]~'= 0)')) + ylab("Doublet")
ggsave("/home/campbell/yulee/project/result/plots/tonsil_nhood_m12_offdiagonal_freq_comparisons.pdf")
