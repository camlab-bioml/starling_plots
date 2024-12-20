rm(list=ls())

options(width = 180)

library(tidyr)
library(dplyr)
library(ggplot2)
library(glue)
library(stringr)
library(ggpubr)
library(ggsci)
library(ggsignif)

library(tidyverse)

library(circlize)
library(hrbrthemes)
library(ComplexHeatmap)

library(reticulate)
np <- import("numpy")

###
###


##fig 2b
neighborhood_analysis <- function(cohort) {
  
  path = paste0("/2processed_files/" + cohort)
  fns <- list.files(path)
  nn_fns <- fns[grep("^nn_", fns)]
  n_fns <- fns[grep("^n_", fns)]
  
  n_df <- NULL
  nn_df <- NULL
  for ( i in 1:length(nn_fns) ) {
    n_df <- rbind(n_df, read.csv(paste0(path, n_fns[i])))
    nn_df <- rbind(nn_df, read.csv(paste0(path, nn_fns[i])))
  }
  n_df[,1] <- 'N'
  nn_df[,1] <- 'NN'
  res <- rbind(n_df, nn_df)
}

bdc <- neighborhood_analysis('basel')
mdc <- neighborhood_analysis('meta')
tdc <- neighborhood_analysis('tonsil')

res <- rbind(cbind(cohort='basel', MESMER='No', bcp), cbind(cohort='basel', MESMER='Yes', bdc), 
             cbind(cohort='meta', MESMER='Yes', mdc), cbind(cohort='tonsil', MESMER='Yes',tdc))

res[,"Arcsinh"] <- ifelse(res[,"Cofactor"] == 0, "No", "Yes")
res[,"Cofactor"] <- ifelse(res[,"Cofactor"] == 0, "N/A", res[,"Cofactor"])
res[res[,"cohort"] == 'basel',"cohort"] <- 'Jackson 2020'
res[res[,"cohort"] == 'meta',"cohort"] <- 'Ali 2020'
res[res[,"cohort"] == 'tonsil',"cohort"] <- 'Tonsil'

res[res[,"Cofactor"] == "N/A","type"] <- "No Arcsinh" 
res[res[,"Cofactor"] == "1","type"] <- "Arcsinh: CF=1"
res[res[,"Cofactor"] == "5","type"] <- "Arcsinh: CF=5"

res[res[,"Initial"] == "PG","Initial"] <- "PhenoGraph"
res[res[,"Initial"] == "KM","Initial"] <- "KMeans"
res[res[,"Initial"] == "FS","Initial"] <- "FlowSOM"
write.csv(filter(res, MESMER == 'Yes', Initial != 'GMM', Cofactor == 5), 
          "fig2b.csv")

rm(list=ls())
df <- read.csv("fig2b.csv")
df %>% ggplot(aes(x = X, y = Score, fill = X, color = X)) +
  facet_grid(rows=glue("Method: {Initial}") ~ cohort) +
  geom_hline(yintercept=.2, linetype="dashed") +
  geom_hline(yintercept=.4, linetype="dashed") +
  geom_hline(yintercept=.6, linetype="dashed") +
  geom_violin(trim = F) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  scale_color_manual(values = c("#E69F00", "#56B4E9")) +
  geom_boxplot(width = 0.1, fill = "white", color = "grey40") +
  geom_signif(comparisons = list(c("N", "NN")), map_signif_level=TRUE, color = 'black', y_position = 0.9) +
  labs(x = 'N: cell has neighbours \n NN: cell has no neighbours', y = 'Plausibility Score') +
  theme_pubr(legend = 'none') +
  theme(text=element_text(size=15)) +
  ylim(0, 1.05)
ggsave("fig2b.pdf", width = 4.5, height = 7)


####
####


##fig 2c
combind_data <- function(cohort, seg_type) {
  
  fns <- list.files(paste0(cohort, "/res/"))
  sfns <- fns[grepl("cs1", fns) & grepl("lv0.1", fns) & grepl("rr1", fns) & grepl("npy", fns)] ## statistics
  cfns <- fns[grepl("cs1", fns) & grepl("lv0.1", fns) & grepl("rr1", fns)] ## cell level predictions grepl("iaKM", fns) &  & grepl(".csv", fns)
  
  df <- NULL
  for ( i in 1:length(fns) ) {
    tmp1 <- np$load(paste0(cohort, "/res/", fns[i]), allow_pickle=TRUE)
    df <- rbind(df, tmp1)
  }
  
  df <- cbind(do.call(rbind, strsplit(df[,1], "_")), df[,-1])
  
  res <- as.data.frame(df)
  colnames(res) <- c("Initial", "Cluster", "Cofactor", "Dist", "Proportion", 
                     "Lambda", "Cell", "Relax", "Run", 
                     "inPS", "stPS", "tsPS", 'sD', 'dD', "stat", "pval", "stF1", "tsF1")
  
  res[,'Initial'] <- substr(res[,'Initial'], 3, 100)
  res[,'Cluster'] <- substr(res[,'Cluster'], 3, 100)
  res[,'Cofactor'] <- substr(res[,'Cofactor'], 3, 100)
  res[,'Dist'] <- substr(res[,'Dist'], 3, 100)
  res[,'Proportion'] <- substr(res[,'Proportion'], 3, 100)
  res[,'Lambda'] <- substr(res[,'Lambda'], 3, 100)
  res[,'Cell'] <- substr(res[,'Cell'], 3, 100)
  res[,'Relax'] <- substr(res[,'Relax'], 3, 100)
  res[,'Run'] <- substr(res[,'Run'], 2, 2)
  
  res[,'Cell'] <- ifelse(res[,'Cell'] == 1, 'Yes', 'No')
  res[,'Relax'] <- ifelse(res[,'Relax'] == 1, 'Yes', 'No')
  res[,'Dist'] <- ifelse(res[,'Dist'] == 0, 'Normal', 'Student-T')
  
  df2 <- NULL
  if (seg_type == 'dc') {
    
    for ( i in 1:length(cfns) ) {
      tmp2 <- read.csv(paste0(cohort, "/res/", cfns[i]))
      df2 <- rbind(df2, cbind(fn=cfns[i], select(tmp2, sample_id:larger_than30, doublet_prob:jett2)))
    }
    df2 <- cbind(do.call(rbind, strsplit(df2[,1], "_")), df2[,-1])
    
    df2[,'Initial'] <- substr(df2[,1], 3, 100)
    df2[,'Cluster'] <- substr(df2[,2], 3, 100)
    df2[,'Cofactor'] <- substr(df2[,3], 3, 100)
    df2[,'Dist'] <- substr(df2[,4], 3, 100)
    df2[,'Proportion'] <- substr(df2[,5], 3, 100)
    df2[,'Lambda'] <- substr(df2[,6], 3, 100)
    df2[,'Cell'] <- substr(df2[,7], 3, 100)
    df2[,'Relax'] <- substr(df2[,8], 3, 100)
    df2[,'Run'] <- substr(df2[,9], 2, 2)
    
    df2[,'Cell'] <- ifelse(df2[,'Cell'] == 1, 'Yes', 'No')
    df2[,'Relax'] <- ifelse(df2[,'Relax'] == 1, 'Yes', 'No')
    df2[,'Dist'] <- ifelse(df2[,'Dist'] == 0, 'Normal', 'Student-T')
  }
  return(list(res, df2))
}

bdc <- combind_data('basel', 'dc')
mdc <- combind_data('meta', 'dc')
tdc <- combind_data('tonsil', 'dc')

res <- rbind(
  cbind(cohort='basel', MESMER='Yes', bdc[[1]]),
  cbind(cohort='meta', MESMER='Yes', mdc[[1]]),
  cbind(cohort='tonsil', MESMER='Yes', tdc[[1]]))

res[,"Arcsinh"] <- ifelse(res[,"Cofactor"] == 0, "No", "Yes")
res[,"Cofactor"] <- ifelse(res[,"Cofactor"] == 0, "N/A", res[,"Cofactor"])

res[res[,"cohort"] == 'basel',"cohort"] <- 'Jackson 2020'
res[res[,"cohort"] == 'meta',"cohort"] <- 'Ali 2020'
res[res[,"cohort"] == 'tonsil',"cohort"] <- 'Tonsil'

res[res[,"Cofactor"] == "N/A","type"] <- "No Arcsinh" 
res[res[,"Cofactor"] == "1","type"] <- "Arcsinh: CF=1"
res[res[,"Cofactor"] == "5","type"] <- "Arcsinh: CF=5"

ps <- res %>%
  group_by(cohort, MESMER, Initial, Cluster, type,  Dist, Lambda, Cell, Relax, Run) %>% 
  summarise('At initialization' = mean(as.numeric(inPS)), 
            'After training' = mean(as.numeric(stPS)), 
            'Two step' = mean(as.numeric(tsPS)), .groups = 'drop') %>%
  pivot_longer(cols = 'At initialization':'Two step', names_to = 'Metric', values_to = 'Score') %>%
  group_by(cohort, MESMER, Initial, type, Dist, Lambda, Cell, Relax, Metric) %>% summarise(PS = mean(Score), .groups = 'drop') 

colnames(ps)[9] <- 'Method'
write.csv(filter(ps, MESMER == 'Yes', Dist == 'Student-T', Cell == 'Yes', Relax == 'Yes'), 
          "fig2c.csv")

ps1 <- as.data.frame(read.csv("fig2c.csv"))
ps1[,"Lambda"] <- as.factor(ps1[,"Lambda"])
df1 <- ps1 %>% filter(Method == "At initialization") 
df2 <- ps1 %>% filter(Method == "After training")
df3 <- ps1 %>% filter(Method == "Two step")

p1 <- ggplot(df2, aes(x = Initial, y = PS, shape=Method)) +
  labs(color = expression(lambda)) + 
  stat_summary(fun=mean, geom="line", aes(group = Lambda, col = Lambda)) +
  stat_summary(fun=mean, geom="point",  aes(group = Lambda, col = Lambda))

p2 <- p1 + 
  stat_summary(data=df1, fun=mean, geom="line", aes(group = 1), linetype = "dotted") + 
  stat_summary(data=df1, fun=mean, geom="point")

p3 <- p2 + 
  stat_summary(data=df3, fun=mean, geom="line", aes(group = 1), linetype = "dashed") + 
  stat_summary(data=df3, fun=mean, geom="point")

p3 + 
  facet_grid(rows = glue("{type}")~glue("{cohort}")) +
  geom_hline(yintercept=.1, linetype="dashed") +
  geom_hline(yintercept=.3, linetype="dashed") +
  geom_hline(yintercept=.5, linetype="dashed") +
  labs(x = "Initialization Algorithm", y = "Plausibility Score") +
  scale_colour_brewer(palette = "Set2") +
  cowplot::theme_cowplot() +
  theme(legend.position="bottom", 
        strip.background = element_rect(fill='white'), 
        strip.text = element_text(face='bold')) + 
  theme_pubr(base_size = 13)
ggsave("fig2c.pdf", height = 7 * (2/3), width = 7)


####
####


##fig 3a (use the benchmark results from 2c analysis)
res <- rbind(cbind(cohort='basel', bdc[[2]]), 
             cbind(cohort='meta', mdc[[2]]), 
             cbind(cohort='tonsil', tdc[[2]]))

res[res[,"larger_than30"] == 0, "larger_than30"] <- "1"
res[res[,"larger_than30"] == 2, "larger_than30"] <- ">1"

res[,"Arcsinh"] <- ifelse(res[,"Cofactor"] == 0, "No", "Yes")
res[,"Cofactor"] <- ifelse(res[,"Cofactor"] == 0, "N/A", res[,"Cofactor"])

res[res[,"cohort"] == 'basel',"cohort"] <- 'Jackson 2020'
res[res[,"cohort"] == 'meta',"cohort"] <- 'Ali 2020'
res[res[,"cohort"] == 'tonsil',"cohort"] <- 'Tonsil'

res[res[,"Cofactor"] == "N/A","type"] <- "No Arcsinh" 
res[res[,"Cofactor"] == "1","type"] <- "Arcsinh: CF=1"
res[res[,"Cofactor"] == "5","type"] <- "Arcsinh: CF=5"
write.csv(res, "fig3a.csv")

df <- read.csv("fig3a.csv")
df %>%
  ggplot(aes(x = factor(larger_than30, level = c("1", ">1")), 
             y = as.numeric(doublet_prob),
             color = larger_than30)) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_hline(yintercept=1, linetype="dashed") +
  ylim(NA, 1.3) + 
  scale_fill_manual(values = c("#D55E00", "#CC79A7", "#009E73")) +
  scale_color_manual(values = c("#D55E00", "#CC79A7", "#009E73")) + #"#999999", "#009E73"
  geom_boxplot(width=0.5) + 
  labs(x = 'Nuclear count', y = 'Segmentation error probability') + 
  geom_signif(comparisons = list(c("1", ">1")), color = 'black', y_position = 1.01) +
  facet_grid(rows = glue("{type}") ~ glue("{cohort}")) + #glue("Dataset: {cohort}")
  theme_pubr(legend = 'none', base_size = 11)
ggsave("fig3a.pdf", width = 3.5, height = 7 * (2/3))


####
####


## fig 3c
res <- rbind(
  cbind(cohort='basel', MESMER='Yes', bdc[[1]]),
  cbind(cohort='meta', MESMER='Yes', mdc[[1]]),
  cbind(cohort='tonsil', MESMER='Yes', tdc[[1]]))

res[res[,"larger_than30"] == 0, "larger_than30"] <- 1

res[,"Arcsinh"] <- ifelse(res[,"Cofactor"] == 0, "No", "Yes")
res[,"Cofactor"] <- ifelse(res[,"Cofactor"] == 0, "N/A", res[,"Cofactor"])

res[res[,"cohort"] == 'basel',"cohort"] <- 'Jackson 2020'
res[res[,"cohort"] == 'meta',"cohort"] <- 'Ali 2020'
res[res[,"cohort"] == 'tonsil',"cohort"] <- 'Tonsil'

res[res[,"Cofactor"] == "N/A","type"] <- "No Arcsinh" 
res[res[,"Cofactor"] == "1","type"] <- "Arcsinh: CF=1"
res[res[,"Cofactor"] == "5","type"] <- "Arcsinh: CF=5"

write.csv(filter(res, MESMER == 'Yes', Lambda == '1.0', Dist == 'Student-T', Cell == 'Yes', Relax == 'Yes'), 
          "fig3c.csv")

df <- read.csv("fig3c.csv")
df %>% 
  ggplot(aes(x = Initial, y = as.numeric(stat), fill = Initial, color = Initial)) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_violin(trim = F) +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9", "#009E73")) +
  scale_color_manual(values = c("#999999", "#E69F00", "#56B4E9", "#009E73")) +
  geom_boxplot(width = 0.1, fill = "white", color = "grey40") + 
  labs(x = 'Initialization algorithm', y = 'Estimate (t-test)') + 
  facet_grid(rows = glue("{type}") ~ glue("{cohort}")) +
  theme_pubr(legend = 'none', base_size = 11)
ggsave("fig3c.pdf", width = 3.5, height = 7 * (2/3))


## fig 3b
rm(list=ls())
mdf <- read_csv("metabrics_gamma.csv")
bdf <- read_csv("basel_gamma.csv")
tdf <- read_csv("tonsil_gamma.csv")

mdf[,"cohort"] <- 'Ali 2020'
bdf[,"cohort"] <- 'Jackson 2020'
tdf[,"cohort"] <- 'Tonsil'

df <- rbind(mdf, bdf, tdf)
write.csv(df, "fig3b.csv")

df <- read.csv("fig3b.csv")
df <- mutate(df,
             scenario = case_when(
               X == "no share" ~ "Cell type not in \n cell neighbours",
               X == "share" ~ "Cell type in \n cell neighbours"
             ))

cols <- c(
  "Cell type not in \n cell neighbours" = "grey50",
  "Cell type in \n cell neighbours" = scales::muted("red")
)

group_by(df, scenario) |> 
  ungroup() |> 
  ggplot(aes(x = scenario, y = prob, color = scenario)) +
  geom_boxplot() +
  scale_y_continuous(breaks=c(0,0.5,1), limits = c(NA, 1.1)) +
  geom_signif(comparisons = list(c("Cell type not in \n cell neighbours", "Cell type in \n cell neighbours")), 
              map_signif_level = TRUE, color = 'black', y_position = 1.01) +
  facet_wrap(.~cohort, ncol = 1) + 
  coord_flip() +
  theme_pubr(legend = 'none', base_size = 11) + 
  labs(y = "Probability of other cell type in \n segmentation error",
       x = "") +
  scale_color_manual(values = cols)
ggsave("fig3b.pdf", width = 3.5, height = 7 * (2/3))


###
###


## fig 5c/supplementary 16
minmax <- function(x, na.rm = TRUE) {
  return((x- min(x, na.rm = TRUE)) /(max(x, na.rm = TRUE)-min(x, na.rm = TRUE)))
}

mat1 <- read.csv("fig5c_init_centroids.csv")[,-1] ## initialization clusters (supplementary 16)
mat2 <- read.csv("fig5c_st_centroids.csv")[,-1]   ## starling clusters
samp_lb <- read.csv("fig5c_labels.csv")[,-1]      ## starling labels

nc <- dim(mat1)[1]

count <- NULL
for (c in 0:(nc-1)) {
  
  not_in_n_are_in  <- samp_lb[,"KM"] != c & samp_lb[,"ST"] == c ## werent in & are in
  were_in_n_are_in <- samp_lb[,"KM"] == c & samp_lb[,"ST"] == c ## were in & are in
  were_in_n_not_in <- samp_lb[,"KM"] == c & samp_lb[,"ST"] != c ## were in & not in
  
  count <- rbind(count, c(sum(not_in_n_are_in), sum(were_in_n_are_in), sum(were_in_n_not_in), 
                          sum(samp_lb[,"KM"] == c),
                          sum(samp_lb[,"ST"] == c)))
}

ct_labels <- c('Endothelial', 'B', 'T', 'Immune', 'B',
               'B', 'T', 'B', 'Mix', 'Stromal',
               'Epithelial', 'B', 'T', 'Macrophage', 'Endothelial',
               'Stromal', 'B', 'Monocytes', 'B', 'Monocytes',
               'Neutrophils', 'B', 'T', 'B', 'Epithelial')

dmat <- apply(rbind(mat1, mat2[,-25]), 2, minmax)

i=0
ma1 <- t(as.matrix(dmat[(i*nc+1):((i+1)*nc),]))
colnames(ma1) <- ct_labels ##  #1:nc #

f1 <- colorRamp2(seq(min(ma1), max(ma1), length = 3), c("blue", "#EEEEEE", "red"), space = "RGB")

haa2 <- HeatmapAnnotation(
  " " = anno_barplot(as.integer(count[,4]), axis = TRUE),
  'Cell Type' = ct_labels, ##rownames(tmp1),
  col = list('Cell Type' = c(
    'Stromal' = '#dbdb8d', 'B' = '#279e68', 'T' = '#8c564b',
    'Immune' = "#999999", 'Mix' = '#ffbb78', 'Epithelial' = '#17becf',
    'Macrophage' = '#aec7e8', 'Endothelial' = '#b5bd61', 
    'Monocytes' = '#c5b0d5','Neutrophils' = '#f7b6d2')),
  annotation_name_side = "left",
  annotation_name_rot = 0, 
  annotation_name_gp = gpar(fontsize = 8), 
  annotation_legend_param = list('Cell Type' = list(title = "Cell Type", 
                                                    labels_gp = gpar(fontsize = 10))))

ct_labels1 <- 1:nc

htt1 <- Heatmap(ma1, 
                name = 'Scaled\nExpression',
                column_title = " ",
                col = f1, #rev(rainbow(3)),
                top_annotation = haa2,
                column_split = colnames(ma1),
                column_labels = ct_labels1,
                cluster_rows = TRUE, cluster_columns = TRUE, show_heatmap_legend = TRUE, border = TRUE)

pdf('fig5c_init.pdf', width = 10) ## see supplementary 16
draw(htt1, column_title = "Cluster", row_title = "Marker", column_title_side = 'bottom', merge_legend = TRUE)
dev.off()

ct_labels <- c('Stromal', 'B', 'T', 'Immune', 'B',
               'B', 'T', 'B', 'Mix', 'Stromal',
               'Epithelial', 'B', 'T', 'Macrophage', 'Endothelial',
               'Stromal', 'B', 'Monocytes', 'B', 'Monocytes',
               'Neutrophils', 'B', 'T', 'B', 'Epithelial')

i=1
ma2 <- t(as.matrix(dmat[(i*nc+1):((i+1)*nc),]))
colnames(ma2) <- ct_labels # #

f2 <- colorRamp2(seq(min(ma2), max(ma2), length = 3), c("blue", "#EEEEEE", "red"), space = "RGB")

column_ha <- HeatmapAnnotation(
  " " = anno_barplot(as.integer(count[,5]), axis = TRUE),
  'Cell Type' = ct_labels, ##rownames(tmp1),
  col = list('Cell Type' = c(
    'Stromal' = '#dbdb8d', 'B' = '#279e68', 'T' = '#8c564b',
    'Immune' = "#999999", 'Mix' = '#ffbb78', 'Epithelial' = '#17becf',
    'Macrophage' = '#aec7e8', 'Endothelial' = '#b5bd61', 
    'Monocytes' = '#c5b0d5','Neutrophils' = '#f7b6d2')),
  annotation_name_side = "left",
  annotation_name_rot = 0, 
  annotation_name_gp = gpar(fontsize = 8), 
  annotation_legend_param = list('Cell Type' = list(title = "Cell Type", 
                                                    labels_gp = gpar(fontsize = 10))))
htt2 <- Heatmap(ma2, 
                name = 'Scaled\nExpression',,
                column_title = " ",
                col = f2, #rev(rainbow(3)),
                top_annotation = column_ha,
                column_split = colnames(ma2),
                column_labels = ct_labels1,
                cluster_rows = TRUE, cluster_columns = TRUE, show_heatmap_legend = TRUE, border = TRUE)

pdf('fig5c.pdf', width = 10)
draw(htt2, column_title = "Cluster", row_title = "Marker", column_title_side = 'bottom', merge_legend = TRUE)
dev.off()


###
###


##fig 5de
rm(list=ls())
col_fun <- colorRamp2(c(0, 1), c("white", "blue"))
lgd <-  Legend(col_fun = col_fun, title = "Proportion", labels_gp = gpar(fontsize = 8))

samp_lb <- read.csv("fig5de.csv")[,-1]
ct_labels <- c('CD4 T', 'CD4 T', 'Unknown', 'Unknown', 'Myeloid',
               'Myeloid', 'CD4 T', 'B', 'CD8 T', 'CD4 T',
               'B', 'Granulocyte', 'CD4 T', 'B', 'Endothelial',
               'B', 'Myeloid', 'CD4 T', 'Endothelial', 'Endothelial',
               'Stromal', 'CD4 T', 'B', 'Myeloid', 'Stromal',
               'Unknown', 'Epithelial', 'Epithelial')

tmp <- NULL
for (ROI in 1:16) {
  
  #i = 1
  select_image = filter(samp_lb, sample_id == ROI)
  
  if (ROI <= 4) {
    patient <- 1
  } else if (ROI <= 12) {
    patient <- 2
  } else {
    patient <- 3
  }
  
  for (c in 0:27) {
    
    tmp <- rbind(tmp, c(ROI, patient, ct_labels[c+1], 
                        c,
                        dim(filter(select_image, label_m1 == c))[1],
                        dim(select_image)[1],
                        dim(filter(select_image, label_m2 == c))[1],
                        dim(select_image)[1] - sum(select_image[,'label_m2'] == 'doublet'),
                        dim(filter(select_image, label_m1 == c))[1] - dim(filter(select_image, label_m2 == c))[1],
                        sum(select_image[,'label_m2'] == 'doublet')
    ))
    
  }
}

colnames(tmp) <- c("ROI", 'P', "CT", 'C', 'Cell_m1', 'Total_m1', 'Cell_m2', 'Total_m2', 'Cell_d', 'Total_d')
res <- as.data.frame(tmp)
res[,"ROI"] <- as.integer(res[,"ROI"])
res[,"C"] <- as.integer(res[,"C"]) + 1
res[,"Cell_m1"] <- as.numeric(res[,"Cell_m1"])
res[,"Cell_m2"] <- as.numeric(res[,"Cell_m2"])
res[,"Cell_d"] <- as.numeric(res[,"Cell_d"])
res[,"Total_m1"] <- as.numeric(res[,"Total_m1"])
res[,"Total_m2"] <- as.numeric(res[,"Total_m2"])
res[,"Total_d"] <- as.numeric(res[,"Total_d"])
res[,"ROI_CP_m1"] <- as.numeric(res[,"Cell_m1"]) / as.numeric(res[,"Total_m1"])
res[,"ROI_CP_m2"] <- as.numeric(res[,"Cell_m2"]) / as.numeric(res[,"Total_m2"])
res[,"ROI_CP_d"] <- as.numeric(res[,"Cell_d"]) / as.numeric(res[,"Total_d"])

tmp2 <- spread(as.data.frame(res[,c("ROI", 'C', 'P', 'ROI_CP_m2')]), C, ROI_CP_m2)
roi <- tmp2[,"ROI"]
p <- tmp2[,'P']
tmp2 <- apply(tmp2[,-c(1,2)], 2, as.numeric)
rownames(tmp2) <- roi

## create plot
ha1 <- rowAnnotation(Donor = as.character(p), 
                     col = list(Donor = c("1" = "#E69F00", "2" = "#0072B2", "3" = "#D55E00")),
                     annotation_name_gp = gpar(fontsize = 6),
                     annotation_name_rot = 0,
                     annotation_name_side = "bottom",
                     annotation_legend_param = list(labels_gp = gpar(fontsize = 6)))

row_ha <- HeatmapAnnotation('Cell Type' = ct_labels, ##rownames(tmp1),
                            col = list('Cell Type' = c(
                              'B' = '#279e68', 'CD4 T' = '#8c564b', 'CD8 T' = '#c5b0d5', 
                              'Endothelial' = '#b5bd61', 'Epithelial' = '#17becf',
                              'Granulocyte' = "#999999", 'Myeloid' = '#aec7e8',
                              'Stromal' = '#dbdb8d', 'Unknown' = '#ffbb78')), #,'Neutrophils' = '#f7b6d2',
                            annotation_name_side = "left",
                            annotation_name_rot = 0,
                            annotation_name_gp = gpar(fontsize = 6), 
                            annotation_legend_param = list('Cell Type' = list(title = "Cell Type", 
                                                                              labels_gp = gpar(fontsize = 6))))

htt1 <- Heatmap(tmp2, col = col_fun, top_annotation = row_ha, left_annotation = ha1,
                cluster_rows = TRUE, cluster_columns = TRUE, show_heatmap_legend = FALSE,
                row_names_gp = gpar(fontsize = 6),
                column_names_gp = gpar(fontsize = 6))

pdf("fig5d.pdf", height = 5, width = 8) 
draw(htt1, annotation_legend_list = lgd, row_title = "ROI", column_title = "Cluster", 
     row_title_gp = gpar(fontsize = 6), column_title_gp = gpar(fontsize = 6),
     column_title_side = 'bottom', merge_legend = TRUE)
dev.off()

tmp2 <- res %>% group_by(ROI, P, CT) %>% summarize_at(c('Cell_m2', 'Total_m2'), c(sum, mean)) %>% 
  mutate(ROI_CP_m2=Cell_m2_fn1/Total_m2_fn2) %>% as.data.frame() %>% select(ROI, CT, P, ROI_CP_m2) %>%
  spread(CT, ROI_CP_m2)

roi <- tmp2[,"ROI"]
p <- tmp2[,'P'] #as.integer(tmp1[,'P'])
tmp2 <- apply(tmp2[,-c(1,2)], 2, as.numeric)
rownames(tmp2) <- roi

ha1 <- rowAnnotation(Donor = as.character(p), 
                     col = list(Donor = c("1" = "#E69F00", "2" = "#0072B2", "3" = "#D55E00")),
                     annotation_name_gp = gpar(fontsize = 6),
                     annotation_name_rot = 0,
                     annotation_name_side = "bottom",
                     annotation_legend_param = list(labels_gp = gpar(fontsize = 6)))

row_ha <- HeatmapAnnotation('Cell Type' = unique(ct_labels), ##rownames(tmp1),
                            col = list('Cell Type' = c(
                              'B' = '#279e68', 'CD4 T' = '#8c564b', 'CD8 T' = '#c5b0d5', 
                              'Endothelial' = '#b5bd61', 'Epithelial' = '#17becf',
                              'Granulocyte' = "#999999", 'Myeloid' = '#aec7e8',
                              'Stromal' = '#dbdb8d', 'Unknown' = '#ffbb78')), #,'Neutrophils' = '#f7b6d2',
                            annotation_name_side = "left",
                            annotation_name_rot = 0, 
                            annotation_name_gp = gpar(fontsize = 6), 
                            annotation_legend_param = list('Cell Type' = list(title = "Cell Type", 
                                                                              labels_gp = gpar(fontsize = 6))))

htt1 <- Heatmap(tmp2, col = col_fun, top_annotation = row_ha, left_annotation = ha1,
                cluster_rows = TRUE, cluster_columns = TRUE, show_heatmap_legend = FALSE,
                row_names_gp = gpar(fontsize = 6),
                column_names_gp = gpar(fontsize = 6))

pdf("fig5e.pdf", height = 5, width = 4)
draw(htt1, annotation_legend_list = lgd, row_title = "ROI", column_title = "Cluster", 
     row_title_gp = gpar(fontsize = 6), column_title_gp = gpar(fontsize = 6),
     column_title_side = 'bottom', merge_legend = TRUE)
dev.off()


###
###


##fig 6ab
rm(list=ls())
p <- c(1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3)
nhood_zsd12 <- NULL
nhood_zsn12 <- NULL
for ( ii in 1:2 ) {
  
  #ii = 1
  zscores <- read.csv(paste0('m', ii, '_zs_roi1.csv'))
  zscores[zscores == 'NaN' | is.na(zscores)] <- 0
  colnames(zscores)[3:4] <- c("CD4 T", "CD8 T")
  
  nhood_zs <- c()
  for ( i in 1:16 ) {
    tmp = filter(zscores, ROI == i)[,-1]
    nhood_zs <- rbind(nhood_zs, tmp[lower.tri(tmp, diag = TRUE)])
  }
  
  rownames(tmp) <- colnames(tmp)
  colnames(nhood_zs) <- t(sapply(1:length(colnames(tmp)), function(ii) paste0(rownames(tmp)[ii], "-", colnames(tmp))))[lower.tri(tmp, diag = TRUE)]
  
  nhood_zs <- as.data.frame(nhood_zs) %>% select(-contains('doublet') & -contains('Mix') & -contains('Unknown') & -contains('IFNg+ B'))
  
  nhood_zs[nhood_zs == 'NaN' | is.na(nhood_zs)] <- 0
  nhood_zs[is.na(nhood_zs)] <- 0
  
  rownames(nhood_zs) <- 1:16 #paste0("ROI", 1:16)
  
  nhood_zsd <- select(nhood_zs, 'B-B', 'CD4 T-CD4 T', 'CD8 T-CD8 T', 'Endothelial-Endothelial', 'Epithelial-Epithelial', 
                      'Granulocyte-Granulocyte', 'Myeloid-Myeloid', 'Stromal-Stromal')
  
  nhood_zsn <- select(nhood_zs, -c('B-B', 'CD4 T-CD4 T', 'CD8 T-CD8 T', 'Endothelial-Endothelial', 'Epithelial-Epithelial', 
                                   'Granulocyte-Granulocyte', 'Myeloid-Myeloid', 'Stromal-Stromal'))
  
  colnames(nhood_zsd) <- do.call(rbind, strsplit(colnames(nhood_zsd), "-"))[,1]
  
  #col_fun <- colorRamp2(seq(0, max(nhood_zsd, na.rm = TRUE), length = 5), c("#fde725", "#5ec962", "#21918c", "#3b528b", "#440154"), space = 'RGB')
  col_fun <- colorRamp2(seq(0, max(nhood_zsd, na.rm = TRUE), length = 5), c("#440154", "#3b528b", "#21918c", "#5ec962", "#fde725"), space = 'RGB')
  
  lgd <- Legend(col_fun = col_fun, title = "Z-Score", labels_gp = gpar(fontsize = 8))
  
  ha1 <- rowAnnotation(Donor = as.character(p), 
                       col = list(Donor = c("1" = "#E69F00", "2" = "#0072B2", "3" = "#D55E00")),
                       annotation_name_rot = 0, annotation_name_gp = gpar(fontsize = 8),
                       annotation_legend_param = list(labels_gp = gpar(fontsize = 8)))
  
  row_ha <- HeatmapAnnotation('Cell Type' = as.character(colnames(nhood_zsd)), 
                              col = list('Cell Type' = c(
                                'B' = '#279e68', 'CD4 T' = '#8c564b', 'CD8 T' = '#c5b0d5', 
                                'Endothelial' = '#b5bd61', 'Epithelial' = '#17becf',
                                'Granulocyte' = "#999999", 'Myeloid' = '#aec7e8',
                                'Stromal' = '#dbdb8d', 'Unknown' = '#ffbb78')), #,'Neutrophils' = '#f7b6d2',
                              annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 8), 
                              annotation_legend_param = list('Cell Type' = list(title = "Cell Type", 
                                                                                labels_gp = gpar(fontsize = 8))))
  
  htt1 <- Heatmap(as.matrix(nhood_zsd), col = col_fun, left_annotation = ha1, 
                  top_annotation = row_ha, #show_row_name = FALSE, 
                  cluster_rows = TRUE, cluster_columns = TRUE, show_heatmap_legend = FALSE,
                  row_names_gp = gpar(fontsize = 8),
                  column_names_gp = gpar(fontsize = 8)) #, column_names_rot = 0
  
  pdf(paste0("fig6a_tonsil_nhood_m", ii, "_zsd.pdf"), width = 4, height = 5)
  draw(htt1, annotation_legend_list = lgd, 
       column_title = "Cell type pair", row_title = "ROI",
       row_title_gp = gpar(fontsize = 8), column_title_gp = gpar(fontsize = 8),
       column_title_side = 'bottom', merge_legend = TRUE)
  dev.off()
  
  ###
  ###
  
  col_fun <- colorRamp2(seq(min(nhood_zsn, na.rm = TRUE), max(nhood_zsn, na.rm = TRUE),, length = 5), c("#440154", "#3b528b", "#21918c", "#5ec962", "#fde725"), space = 'RGB')
  
  lgd <- Legend(col_fun = col_fun, title = "Z-Score", labels_gp = gpar(fontsize = 8))
  
  ha1 <- rowAnnotation(Donor = as.character(p), 
                       col = list(Donor = c("1" = "#E69F00", "2" = "#0072B2", "3" = "#D55E00")),
                       annotation_name_rot = 0, annotation_name_gp= gpar(fontsize = 8),
                       annotation_legend_param = list(labels_gp = gpar(fontsize = 8)))
  
  row_ha <- HeatmapAnnotation('Type 1' = as.character(do.call(rbind, strsplit(colnames(nhood_zsn), "-"))[,1]),
                              'Type 2' = as.character(do.call(rbind, strsplit(colnames(nhood_zsn), "-"))[,2]),
                              col = list('Type 1'= c(
                                'B' = '#279e68', 'CD4 T' = '#8c564b', 'CD8 T' = '#c5b0d5', 
                                'Endothelial' = '#b5bd61', 'Epithelial' = '#17becf',
                                'Granulocyte' = "#999999", 'Myeloid' = '#aec7e8',
                                'Stromal' = '#dbdb8d', 'Unknown' = '#ffbb78'),
                                'Type 2'= c(
                                  'B' = '#279e68', 'CD4 T' = '#8c564b', 'CD8 T' = '#c5b0d5', 
                                  'Endothelial' = '#b5bd61', 'Epithelial' = '#17becf',
                                  'Granulocyte' = "#999999", 'Myeloid' = '#aec7e8',
                                  'Stromal' = '#dbdb8d', 'Unknown' = '#ffbb78')),
                              annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 8), 
                              annotation_legend_param = list('Type 1' = list(title = "Cell Type 1", 
                                                                             labels_gp = gpar(fontsize = 8)),
                                                             'Type 2' = list(title = "Cell Type 2", 
                                                                             labels_gp = gpar(fontsize = 8))))
  
  htt1 <- Heatmap(as.matrix(nhood_zsn), col = col_fun, left_annotation = ha1, 
                  top_annotation = row_ha, #show_row_name = FALSE, 
                  #right_annotation = row_ha2,
                  cluster_rows = TRUE, cluster_columns = TRUE, show_heatmap_legend = FALSE,
                  row_names_gp = gpar(fontsize = 8),
                  column_names_gp = gpar(fontsize = 8)) #, column_names_rot = 0)
  
  pdf(paste0("fig6a_tonsil_nhood_m", ii, "_zsn.pdf"), width = 10, height = 6)
  draw(htt1, annotation_legend_list = lgd, column_title = "Cell type pair", row_title = "ROI", 
       row_title_gp = gpar(fontsize = 8), column_title_gp = gpar(fontsize = 8),
       column_title_side = 'bottom', merge_legend = TRUE)
  dev.off()
  
  nhood_zsd12 <- cbind(nhood_zsd12, unlist(nhood_zsd))
  nhood_zsn12 <- cbind(nhood_zsn12, unlist(nhood_zsn))
}

colnames(nhood_zsd12) <- c("method1", "method2")
colnames(nhood_zsn12) <- c("method1", "method2")

ggplot(as.data.frame(nhood_zsd12), aes(x=method1, y=method2)) + geom_point(alpha=0.5) + 
  stat_cor(method = "pearson") + 
  labs(x = "Enrichment assigning all cells", y = "Enrichment removing poorly\nsegmented cells") +
  theme_pubr(base_size = 8) #base_size = 6
ggsave("fig6b_tonsil_nhood_m12_diagonal_zs_comparisons.pdf", height = 2.06, width = 2.32)

ggplot(as.data.frame(nhood_zsn12), aes(x=method1, y=method2)) + geom_point(alpha=0.5) + 
  stat_cor(method = "pearson") + 
  labs(x = "Enrichment assigning all cells", y = "Enrichment removing poorly\nsegmented cells") +
  theme_pubr(base_size = 8) #
ggsave("fig6b_tonsil_nhood_m12_offdiagonal_zs_comparisons.pdf", height = 2.06, width = 2.32)