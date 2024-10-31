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

## fig 1
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

write.csv(df, "sup_fig1.csv")

ggplot(df, aes(x=asinh(CD3/5), y=asinh(CD20/5))) + geom_point() + facet_grid(cols = vars(cohort)) + 
  theme(plot.caption = element_text(hjust = 0.5)) +
  theme_pubr() + 
  ylim(c(0, 5)) + 
  labs(caption = "Supplementary Figure 1. Scatterplot of CD3 and CD20 across three IMC datasets.",
       x = "CD3", y = "CD20")
ggsave("sfig1.png", width = 5, height = 5)


###
###


## fig 2 & 3
rm(list=ls())
neighborhood_analysis <- function(cohort, seg_type, type) {
  
  path <- paste0(cohort, "/analysis/_1nh/", seg_type, "/")
  fns <- list.files(path)
  if ( type == 'random' ) {
    nn_fns <- fns[grepl("^nn_", fns) & grepl("_ctr.csv", fns)]
    n_fns <- fns[grepl("^n_", fns) & grepl("_ctr.csv", fns)]
    
  } else {
    nn_fns <- fns[grepl("^nn_", fns) & grepl("_ct.csv", fns)]
    n_fns <- fns[grepl("^n_", fns) & grepl("_ct.csv", fns)]
  }
  
  n_df <- NULL
  nn_df <- NULL
  for ( i in 1:length(n_fns) ) {
    n_df <- rbind(n_df, cbind(i, read.csv(paste0(path, n_fns[i]))))
    nn_df <- rbind(nn_df, cbind(i, read.csv(paste0(path, nn_fns[i]))))
  }
  n_df[,2] <- 'N'
  nn_df[,2] <- 'NN'
  res <- rbind(n_df, nn_df)
}

## fig 2
meta <- neighborhood_analysis('meta', 'dc', 's')
basel <- neighborhood_analysis('basel', 'dc', 's')
tonsil <- neighborhood_analysis('tonsil', 'dc', 's')

## fig 3
meta <- neighborhood_analysis('meta', 'dc', 'r')
basel <- neighborhood_analysis('basel', 'dc', 'r')
tonsil <- neighborhood_analysis('tonsil', 'dc', 'r')

res <- rbind(cbind(cohort='Jackson 2020', basel), cbind(cohort='Ali 2020', meta), cbind(cohort='Tonsil', tonsil))

out <- NULL
for (i in 1:nrow(res)) {
  
  tmp <- res[i,]
  cts <- strsplit(str_replace_all(tmp['ct'], "[[:punct:]]", ""), " ")[[1]]
  row1 <- tmp; row1['ct'] <- cts[1]
  
  if ( str_count(tmp['ct'], "[A-Z]") == 2 ) {
    row2 <- tmp; row2['ct'] <- cts[2]
    out <- rbind(out, rbind(row1, row2))
  } else {
    out <- rbind(out, row1)
  }
}

colnames(out)[8:10] <- c("Cluster", "Cofactor", "Initial")
write.csv(out, "fig2bCT.csv") ## write.csv(out, "fig2bCTr.csv")

rm(list=ls())
out <- read.csv("fig2bCT.csv") ## out <- read.csv("fig2bCTr.csv")
out[out[,"X"] == "N",'X'] <- "cell has neighbours"
out[out[,"X"] == "NN",'X'] <- "cell has no neighbours"

as.data.frame(out) %>% filter(Cofactor == 5) %>%
  ggplot(aes(x = ct, y = Score, fill = X)) +
  #geom_boxplot(width = 0.1, fill = "white", color = "grey40") + 
  geom_boxplot() + 
  facet_grid(rows=glue("Method: {Initial}") ~ cohort) +
  geom_hline(yintercept=.25, linetype="dashed") +
  geom_hline(yintercept=.5, linetype="dashed") +
  geom_hline(yintercept=.75, linetype="dashed") +
  #geom_violin(trim = F) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  scale_color_manual(values = c("#E69F00", "#56B4E9")) +
  #geom_signif(comparisons = list(c("N", "NN")), map_signif_level = TRUE, color = "black") + 
  labs(x = 'Cell Type', y = 'Plausibility Score') + 
  #labs(x = 'N: cell has neighbours \n NN: cell has no neighbours', y = 'Plausibility Score') + 
  theme_pubr(legend = 'none') +
  theme(text=element_text(size=25), legend.position="top", legend.title=element_blank())
ggsave("sfig2.pdf", width = 40, height = 15) ## ggsave("sfig3.pdf", width = 40, height = 15)


###
###


##fig 4
combind_data_revision <- function(cohort, seg_type) {
  
  fns <- list.files(paste0(cohort, "/analysis/_2bm/", seg_type, "/res"))
  cfns <- fns[grep(".csv", fns)]
  cfns <- cfns[grep("ia", cfns)]
  fns <- fns[grep("npy", fns)]
  
  df <- NULL
  for ( i in 1:length(fns) ) {
    tmp1 <- np$load(paste0(cohort, "/analysis/_2bm/", seg_type, "/res/", fns[i]), allow_pickle=TRUE)
    df <- rbind(df, tmp1)
  }
  
  df <- cbind(do.call(rbind, strsplit(df[,1], "_")), df[,-1])
  
  res <- as.data.frame(df)
  colnames(res) <- c("Initial", "Cluster", "Cofactor", "Dist", "Proportion", 
                     "Lambda", "Cell", "Relax", "Run", "threshold", "inPS", "stPS", "tsPS")
  
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
  
  return(res)
}

bdc <- combind_data_revision('basel', 'dc')
mdc <- combind_data_revision('meta', 'dc')
tdc <- combind_data_revision('tonsil', 'dc')

res <- rbind(
  cbind(cohort='basel', MESMER='Yes', bdc),
  cbind(cohort='meta', MESMER='Yes', mdc),
  cbind(cohort='tonsil', MESMER='Yes', tdc))

res[,"Arcsinh"] <- ifelse(res[,"Cofactor"] == 0, "No", "Yes")
res[,"Cofactor"] <- ifelse(res[,"Cofactor"] == 0, "N/A", res[,"Cofactor"])

res[res[,"cohort"] == 'basel',"cohort"] <- 'Jackson 2020'
res[res[,"cohort"] == 'meta',"cohort"] <- 'Ali 2020'
res[res[,"cohort"] == 'tonsil',"cohort"] <- 'Tonsil'

res[res[,"Cofactor"] == "N/A","type"] <- "No Arcsinh" 
res[res[,"Cofactor"] == "1","type"] <- "Arcsinh: CF=1"
res[res[,"Cofactor"] == "5","type"] <- "Arcsinh: CF=5"

res[res[,"Cell"] == "Yes","Cell"] <- "Size: Yes" 
res[res[,"Cell"] == "No","Cell"] <- "No Size"

res[res[,"Relax"] == "Yes","Relax"] <- "Relax: Yes" 
res[res[,"Relax"] == "No","Relax"] <- "No Relax"

ps <- res %>%
  group_by(cohort, MESMER, Initial, Cluster, type,  Dist, Lambda, Cell, Relax, Run, threshold) %>% 
  summarise('At initialization' = mean(as.numeric(inPS)), 
            'After training' = mean(as.numeric(stPS)), 
            'Two step' = mean(as.numeric(tsPS)), .groups = 'drop') %>%
  pivot_longer(cols = 'At initialization':'Two step', names_to = 'Metric', values_to = 'Score') %>%
  group_by(cohort, MESMER, Initial, type, Dist, Lambda, Cell, Relax, Metric, threshold) %>% summarise(PS = mean(Score), .groups = 'drop')

colnames(ps)[9] <- 'Method'

df1 <- ps %>% filter(Method == "At initialization") 
df2 <- ps %>% filter(Method == "After training")
df3 <- ps %>% filter(Method == "Two step")

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
  facet_grid(rows = glue("{type}")+glue("{Dist}")+glue("{Cell}")+glue("{Relax}")~glue("{cohort}")+glue("{threshold}")) +
  labs(x = "Initialization Algorithm", y = "Plausibility Score") +
  scale_colour_brewer(palette = "Set2") +
  cowplot::theme_cowplot() +
  theme(legend.position="bottom", 
        strip.background = element_rect(fill='white'), 
        strip.text = element_text(face='bold')) + 
  theme_pubr(base_size = 12)
ggsave("sfig4.pdf", width = 18, height = 10)


###
###


##fig 5
ps5 <- read.csv("fig2cCP.csv")
ps5[,"PS"] <- as.numeric(ps5[,"PS"])
ps5[,"Lambda"] <- as.character(ps5[,"Lambda"])
df1 <- ps5 %>% filter(Method == "At initialization") 
df2 <- ps5 %>% filter(Method == "After training")
df3 <- ps5 %>% filter(Method == "Two step")

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
  geom_hline(yintercept=.2, linetype="dashed") +
  geom_hline(yintercept=.4, linetype="dashed") +
  geom_hline(yintercept=.6, linetype="dashed") +
  labs(x = "Initialization Algorithm", y = "Plausibility Score") +
  scale_colour_brewer(palette = "Set2") +
  cowplot::theme_cowplot() +
  theme(legend.position="bottom", 
        strip.background = element_rect(fill='white'), 
        strip.text = element_text(face='bold')) + 
  theme_pubr(base_size = 13)
ggsave("sfig5.pdf", height = 7 * (2/3), width = 7)


###
###


##fig 6
thresholds <- as.data.frame(matrix(c("Ali 2020", "PG", "N/A", 0.1014217, 0.3330025,
                                     "Ali 2020", "FS", "N/A", 0.1826533, 0.283738,
                                     "Ali 2020", "KM", "N/A", 0.1687988, 0.2594652,
                                     "Ali 2020", "PG", "1", 0.1886188, 0.3974577,
                                     "Ali 2020", "FS", "1", 0.2622018, 0.459409,
                                     "Ali 2020", "KM", "1", 0.2047815, 0.325212,
                                     "Ali 2020", "PG", "5", 0.1621243, 0.40095,
                                     "Ali 2020", "FS", "5", 0.1897675, 0.4377308,
                                     "Ali 2020", "KM", "5", 0.1393599, 0.2751844,
                                     "Jackson 2020", "PG", "N/A", 0.1051123, 0.2670303,
                                     "Jackson 2020", "FS", "N/A", 0.1442789, 0.2213176,
                                     "Jackson 2020", "KM", "N/A", 0.1880208, 0.3096692,
                                     "Jackson 2020", "PG", "1", 0.2054446, 0.3968784,
                                     "Jackson 2020", "FS", "1", 0.2292674, 0.3554764,
                                     "Jackson 2020", "KM", "1", 0.1720473, 0.3442412,
                                     "Jackson 2020", "PG", "5", 0.1301132, 0.3541478,
                                     "Jackson 2020", "FS", "5", 0.1556333, 0.2749599,
                                     "Jackson 2020", "KM", "5", 0.1478738, 0.2599797,
                                     "Tonsil", "PG", "N/A", 0.1007285, 0.5655002,
                                     "Tonsil", "FS", "N/A", 0.1510936, 0.3930774,
                                     "Tonsil", "KM", "N/A", 0.09918653, 0.4334093,
                                     "Tonsil", "PG", "1", 0.2456372, 0.7363159,
                                     "Tonsil", "FS", "1", 0.3891866, 0.719849,
                                     "Tonsil", "KM", "1", 0.2511514, 0.7116285,
                                     "Tonsil", "PG", "5", 0.2065799, 0.7082042,
                                     "Tonsil", "FS", "5", 0.3066057, 0.6283785,
                                     "Tonsil", "KM", "5", 0.2358396, 0.5832926), ncol=5, byrow=TRUE))

colnames(thresholds) <- c("cohort", "Initial", "Cofactor", "lower", "upper")
thresholds[thresholds[,"Cofactor"] == "N/A","type"] <- "No Arcsinh" 
thresholds[thresholds[,"Cofactor"] == "1","type"] <- "Arcsinh: CF=1"
thresholds[thresholds[,"Cofactor"] == "5","type"] <- "Arcsinh: CF=5"

thresholds <- pivot_longer(thresholds, cols = 'lower':'upper', names_to = 'Method', values_to = 'PS')[,-3]
ps1 <- as.data.frame(read.csv("fig2c.csv"))
ps1[,"Lambda"] <- as.factor(ps1[,"Lambda"])
ps2 <- merge(ps1, thresholds, by = c("cohort", "Initial", 'type'), all.x = TRUE)
ps2[,"lower"] <- as.numeric(ps2[,"lower"])
ps2[,"upper"] <- as.numeric(ps2[,"upper"])
ps3 <- as.data.frame(ps2) %>% 
  group_by(cohort, Initial, type, Lambda, Cell, Relax, Method) %>% 
  summarise(x_scaled = ((PS - lower) / (upper - lower)), .groups = 'drop') %>%
  mutate(xs = (x_scaled - min(x_scaled)) / (max(x_scaled) - min(x_scaled)))

#((PS - lower) / (upper - lower)) * (upper - lower) + min
df1 <- ps3 %>% filter(Method == "At initialization") 
df2 <- ps3 %>% filter(Method == "After training")
df3 <- ps3 %>% filter(Method == "Two step")

p1 <- ggplot(df2, aes(x = Initial, y = xs, shape=Method)) +
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
  geom_hline(yintercept=0, linetype="dashed") +
  geom_hline(yintercept=.5, linetype="dashed") +
  geom_hline(yintercept=1, linetype="dashed") +
  labs(x = "Initialization Algorithm", y = "Plausibility Score") +
  scale_colour_brewer(palette = "Set2") +
  cowplot::theme_cowplot() +
  theme(legend.position="bottom", 
        strip.background = element_rect(fill='white'), 
        strip.text = element_text(face='bold')) + 
  theme_pubr(base_size = 13)
ggsave("sfig6.pdf", height = 7 * (2/3), width = 7)


###
###


##fig 7
ps4 <- ps1[,-1] %>% 
  pivot_wider(names_from = Method, values_from = PS) %>%
  mutate('STARLING' = 100 * (`After training` - `At initialization`) / `At initialization`, 
         'Two Step' = 100 * (`Two step` - `At initialization`) / `At initialization`) %>%
  select(cohort:Relax, 'STARLING':'Two Step') %>%
  pivot_longer(cols = 'STARLING':'Two Step', names_to = 'Method', values_to = 'PS')

df2 <- ps4 %>% filter(Method == "STARLING")
df3 <- ps4 %>% filter(Method == "Two Step")

p1 <- ggplot(df2, aes(x = Initial, y = PS, shape=Method)) +
  labs(color = expression(lambda)) + 
  stat_summary(fun=mean, geom="line", aes(group = Lambda, col = Lambda)) +
  stat_summary(fun=mean, geom="point",  aes(group = Lambda, col = Lambda))

p2 <- p1 + 
  stat_summary(data=df3, fun=mean, geom="line", aes(group = 1, color = 'Two step')) + 
  stat_summary(data=df3, fun=mean, geom="point", aes(group = 1, color = 'Two step'))

p2 + 
  facet_grid(rows = glue("{type}")~glue("{cohort}")) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_hline(yintercept=.5, linetype="dashed") +
  geom_hline(yintercept=1, linetype="dashed") +
  labs(x = "Initialization Algorithm", y = "% Increase") +
  scale_colour_brewer(palette = "Set2") +
  cowplot::theme_cowplot() +
  theme(legend.position="bottom", 
        strip.background = element_rect(fill='white'), 
        strip.text = element_text(face='bold')) + 
  theme_pubr(base_size = 13)
ggsave("sfig7.pdf", height = 7 * (2/3), width = 7)


###
###


##fig 8
rm(list=ls())
path <- "/home/campbell/yulee/project/st/codex/analysis/_2bm/cp/"
fns <- list.files(path)
wfns <- fns[grepl("r[0-9].csv", fns)]
res <- NULL
for (i in 1:length(wfns)) { res <- rbind(res, read.csv(paste0(path, wfns[i])))}
mdn <- do.call(rbind, strsplit(res[,"X0"], "_"))
mdn <- cbind(substr(mdn[,1:8], 3, 10), substr(mdn[,9], 2, 10))
res <- cbind(mdn, select(res, X1:X5))
colnames(res) <- c("Initial", "k", "Cofactor", "Dist", "singlet.prop", "Lambda", 
                   "cellsize.name", "overlap", "Run", "At initialization", "After training", "Two step", "rf.f1", "st.f1")


res[res[,"Cofactor"] == "0","type"] <- "No Arcsinh" 
res[res[,"Cofactor"] == "1","type"] <- "Arcsinh: CF=1"
res[res[,"Cofactor"] == "5","type"] <- "Arcsinh: CF=5"
#res <- res[(res[,'At initialization'] < res[,'After training']),]

ps <- as.data.frame(res) %>% #filter(Initial != 'FS') %>%
  select(c(Initial, type, Lambda, Run, "At initialization":"Two step")) %>%
  pivot_longer(cols = "At initialization":"Two step", names_to = 'Method', values_to = 'PS') %>%
  group_by(Initial, type, Lambda, Method) %>% summarise(PS = mean(PS), .groups='drop')

write.csv(ps, "sup_fig8.csv")
head(ps[order(ps$PS, decreasing=TRUE),], 15)

df1 <- ps %>% filter(Method == "At initialization") 
df2 <- ps %>% filter(Method == "After training")
df3 <- ps %>% filter(Method == "Two step")

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
  facet_grid(rows = glue("{type}")~.) +
  geom_hline(yintercept=.2, linetype="dashed") +
  geom_hline(yintercept=.4, linetype="dashed") +
  geom_hline(yintercept=.6, linetype="dashed") +
  labs(x = "Initialization Algorithm", y = "Plausibility Score") +
  scale_colour_brewer(palette = "Set2") +
  cowplot::theme_cowplot() +
  theme(legend.position="bottom", 
        strip.background = element_rect(fill='white'), 
        strip.text = element_text(face='bold')) + 
  theme_pubr(base_size = 13)
ggsave("sfig8.pdf")


###
###


##fig 9 and 10
it <- read.csv("sfig9_heatmap.csv")
st <- read.csv("sfig10_heatmap.csv")
asg <- read.csv("sfig9n10_labels.csv")

nc <- dim(it)[1]

count <- NULL
for (c in 0:(nc-1)) { count <- rbind(count, cbind(sum(asg[,"FS"] == c), sum(asg[,"FS_st"] == c))) }

dmat <- rbind(apply(it[,-1], 2, minmax), apply(st[,-1], 2, minmax))

i=0
ma1 <- t(as.matrix(dmat[(i*nc+1):((i+1)*nc),]))
colnames(ma1) <- 1:nc

f1 <- colorRamp2(seq(min(ma1), max(ma1), length = 3), c("blue", "#EEEEEE", "red"), space = "RGB")

haa1 <- HeatmapAnnotation(
  " " = anno_barplot(as.integer(count[,1]), axis = TRUE, ylim = c(0, 3500)),
  annotation_name_side = "left",
  annotation_name_rot = 0, 
  annotation_name_gp = gpar(fontsize = 8))

htt1 <- Heatmap(ma1, 
                name = 'Scaled\nExpression',
                column_title = " ",
                col = f1, #rev(rainbow(3)),
                top_annotation = haa1,
                #column_split = colnames(ma1),
                column_labels =  1:nc,
                cluster_rows = TRUE, cluster_columns = TRUE, show_heatmap_legend = TRUE, border = TRUE)

pdf('sfig9.pdf', height = 6)
draw(htt1, column_title = "Cluster", row_title = "Marker", column_title_side = 'bottom', merge_legend = TRUE)
dev.off()

i=1
ma2 <- t(as.matrix(dmat[(i*nc+1):((i+1)*nc),]))
colnames(ma2) <- 1:nc
f2 <- colorRamp2(seq(min(ma2), max(ma2), length = 3), c("blue", "#EEEEEE", "red"), space = "RGB")

haa2 <- HeatmapAnnotation(
  " " = anno_barplot(as.integer(count[,2]), axis = TRUE, ylim = c(0, 3500)),
  annotation_name_side = "left",
  annotation_name_rot = 0, 
  annotation_name_gp = gpar(fontsize = 8))

htt2 <- Heatmap(ma2, 
                name = 'Scaled\nExpression',
                column_title = " ",
                col = f2, #rev(rainbow(3)),
                top_annotation = haa2,
                column_labels =  1:nc,
                cluster_rows = TRUE, cluster_columns = TRUE, show_heatmap_legend = TRUE, border = TRUE)

pdf('sfig10.pdf', height = 6)
draw(htt2, column_title = "Cluster", row_title = "Marker", column_title_side = 'bottom', merge_legend = TRUE)
dev.off()


###
###


##fig 11
df <- read.csv("sfig11.csv")
res2 %>%
  ggplot(aes(x = factor(larger_than30, level = c("1", ">1")), 
             y = as.numeric(rf1_FS_dp),
             color = larger_than30)) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_hline(yintercept=1, linetype="dashed") +
  #geom_violin(trim = F) +
  ylim(NA, 1.3) + 
  scale_fill_manual(values = c("#D55E00", "#CC79A7", "#009E73")) +
  scale_color_manual(values = c("#D55E00", "#CC79A7", "#009E73")) + #"#999999", "#009E73"
  geom_boxplot(width=0.5) + 
  labs(x = 'Nuclear count', y = 'Segmentation error probability') + 
  geom_signif(comparisons = list(c("1", ">1")), color = 'black', y_position = 1.01) +
  facet_grid(rows = glue("{type}") ~ glue("{cohort}")) + #glue("Dataset: {cohort}") #
  theme_pubr(legend = 'none', base_size = 11)
ggsave("sfig11.pdf", width = 3.5, height = 7 * (2/3))


###
###


##fig 13
rm(list=ls())
mc <- read.csv("sfig13.csv")
mc$Cluster <- as.factor(mc$Cluster)
ggplot(mc, aes(x=Cluster, y=Stat)) + geom_bar(stat = "identity") + labs(y = "Estimate (t-test)")
ggsave("sfig13.pdf", height = 4, width = 5)