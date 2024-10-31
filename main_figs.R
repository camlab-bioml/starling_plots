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

# Figure 2b (_1nh_paper.pdf)
neighborhood_analysis <- function(cohort, seg_type) {
  
  cohort = 'meta'
  #seg_type = 'dc'
  if (cohort == 'basel') {
    path <- paste0("/home/campbell/yulee/project/st/", cohort, "/analysis/_1nh/", seg_type, "/v1/")
  } else {
    path <- paste0("/home/campbell/yulee/project/st/", cohort, "/analysis/_1nh/", seg_type, "/")
  }
  
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

#bcp <- neighborhood_analysis('basel', 'cp') ##v1
bdc <- neighborhood_analysis('basel', 'dc') ##v1
mdc <- neighborhood_analysis('meta', 'dc')
tdc <- neighborhood_analysis('tonsil', 'dc')

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
          "/home/campbell/yulee/project/st/plot/final/fig2b.csv")

ch = 'Tonsil' 
cf = "5"

df <- filter(res, MESMER == 'Yes', Initial != 'GMM', Cofactor == cf)

mean(filter(df, cohort == ch & X == 'N' & Initial == 'PhenoGraph')[,"Score"])
mean(filter(df, cohort == ch & X == 'NN' & Initial == 'PhenoGraph')[,"Score"])

mean(filter(df, cohort == ch & X == 'N' & Initial == 'FlowSOM')[,"Score"])
mean(filter(df, cohort == ch & X == 'NN' & Initial == 'FlowSOM')[,"Score"])

mean(filter(df, cohort == ch & X == 'N' & Initial == 'KMeans')[,"Score"])
mean(filter(df, cohort == ch & X == 'NN' & Initial == 'KMeans')[,"Score"])

## create figure
rm(list=ls())
setwd("~/Desktop/plot/")
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

##Figure 2b (revise)
rm(list=ls())
neighborhood_analysis <- function(cohort, seg_type, type) {
  
  #cohort = 'tonsil'
  #seg_type = 'dc'
  path <- paste0("/home/campbell/yulee/project/st/", cohort, "/analysis/_1nh/", seg_type, "/")
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

meta <- neighborhood_analysis('meta', 'dc', 'r')
basel <- neighborhood_analysis('basel', 'dc', 'r')
tonsil <- neighborhood_analysis('tonsil', 'dc', 'r')
res <- rbind(cbind(cohort='Jackson 2020', basel), cbind(cohort='Ali 2020', meta), cbind(cohort='Tonsil', tonsil))

#ct <- c("B", "Endothelial", "Epithelial", "Immune", "Myeloid", "Stromal", "T")
#unique(do.call(c, strsplit(str_replace_all(basel[,'ct'], "[[:punct:]]", ""), " ")))

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
write.csv(out, "fig2bCT.csv")
#write.csv(out, "fig2bCTr.csv")

rm(list=ls())
out <- read.csv("~/Downloads/fig2bCTr.csv")
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
ggsave("fig2bCT.pdf", width = 40, height = 15)


# Figure 2c (_2bm_ps_paper.pdf)
combind_data <- function(cohort, seg_type) {
  
  #cohort = 'basel'
  #seg_type = 'cp'
  fns <- list.files(paste0("/home/campbell/yulee/project/st/", cohort, "/analysis/_2bm/", seg_type, "/res"))
  #fns <- list.files(paste0("/ddn_exa/campbell/jlee/project/1st/st/", cohort, "/analysis/_2bm/", seg_type, "/res"))
  sfns <- fns[grepl("cs1", fns) & grepl("lv0.1", fns) & grepl("rr1", fns) & grepl("npy", fns)] ## statistics
  cfns <- fns[grepl("cs1", fns) & grepl("lv0.1", fns) & grepl("rr1", fns)] ## cell level predictions grepl("iaKM", fns) &  & grepl(".csv", fns)
  
  df <- NULL
  for ( i in 1:length(fns) ) {
    tmp1 <- np$load(paste0("/home/campbell/yulee/project/st/", cohort, "/analysis/_2bm/", seg_type, "/res/", fns[i]), allow_pickle=TRUE)
    #tmp1 <- np$load(paste0("/ddn_exa/campbell/jlee/project/1st/st/", cohort, "/analysis/_2bm/", seg_type, "/res/", fns[i]), allow_pickle=TRUE)
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
      #tmp2 <- read.csv(paste0("/home/campbell/yulee/project/st/", cohort, "/analysis/_2bm/", seg_type, "/res/", cfns[i]))
      tmp2 <- read.csv(paste0("/ddn_exa/campbell/jlee/project/1st/st/", cohort, "/analysis/_2bm/", seg_type, "/res/", cfns[i]))
      df2 <- rbind(df2, cbind(fn=cfns[i], select(tmp2, sample_id:larger_than30, doublet_prob:jett2)))
    }
    df2 <- cbind(do.call(rbind, strsplit(df2[,1], "_")), df2[,-1])
    #res2 <- as.data.frame(df2)
    #colnames(res2) <- c("Initial", "Cluster", "Cofactor", "Dist", "Proportion", "Lambda", "Cell", "Relax", "Run", 
    #                    "index", "sample_id", "cell_id", 'kieran', 'jett', "doublet_prob", "jett2")
    
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

res <- rbind(
  cbind(cohort='basel', MESMER='No', bdc[[1]]),
  cbind(cohort='meta', MESMER='No', mdc[[1]]))

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
          "/home/campbell/yulee/project/st/plot/final/fig2c.csv")

write.csv(ps, "/home/campbell/yulee/project/st/plot/final/fig2cCP.csv")

## create figure
#rm(list=ls())
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

#rm(list=ls())

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
#thresholds[,"MESMER"] <- "Yes"
#thresholds[,"Lambda"] <- "N/A"
#thresholds[,"Dist"] <- "Student-T"
#thresholds[,"Cell"] <- "Yes"
#thresholds[,"Relax"] <- "Yes"

thresholds <- pivot_longer(thresholds, cols = 'lower':'upper', names_to = 'Method', values_to = 'PS')[,-3]
ps1 <- as.data.frame(read.csv("fig2c.csv"))

filter(ps1, Method == 'After training')

ps1[,"Lambda"] <- as.factor(ps1[,"Lambda"])

ps1 <- merge(ps1, thresholds, by = c("cohort", "Initial", "type", "Method", "PS"), all = TRUE)
ps1[,"PS"] <- as.numeric(ps1[,"PS"])

#filter(df2, cohort == 'Ali 2020', Initial == 'KM', type == 'Arcsinh: CF=5', Lambda == 1)

#mean(filter(df2, cohort == 'Tonsil', Initial %in% c('KM', "PG"), type == 'Arcsinh: CF=5')[,"PS"])
#mean(filter(df2, cohort == 'Ali 2020', Initial %in% c('KM', "PG"), type == 'Arcsinh: CF=5')[,"PS"])
#mean(filter(df2, cohort == 'Jackson 2020', Initial %in% c('KM', "PG"), type == 'Arcsinh: CF=5')[,"PS"])

#mean(filter(df2, cohort == 'Tonsil', Initial %in% 'KM', type == 'Arcsinh: CF=5')[,"PS"])
#mean(filter(df2, cohort == 'Ali 2020', Initial %in% 'PG', type == 'Arcsinh: CF=5')[,"PS"])
#mean(filter(df2, cohort == 'Jackson 2020', Initial %in% "KM", type == 'Arcsinh: CF=5')[,"PS"])

#stat_summary(data=thresholds, fun=mean, geom="line", aes(group = Method, col = Method)) +
#  stat_summary(data=thresholds, fun=mean, geom="point",  aes(group = Method, col = Method))

df1 <- ps1 %>% filter(Method == "At initialization") 
df2 <- ps1 %>% filter(Method == "After training")
df3 <- ps1 %>% filter(Method == "Two step")
df4 <- ps1 %>% filter(Method == "lower") %>% select(cohort:PS)
df5 <- ps1 %>% filter(Method == "upper") %>% select(cohort:PS)

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

p4 <- p3 + 
  stat_summary(data=df4, fun=mean,  colour = "green", size = 2, geom="point")
  
p5 <- p4 + 
  stat_summary(data=df5, fun=mean,  colour = "red", size = 2, geom="point") 
  
p5 + 
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
ggsave("fig2c.pdf", height = 7 * (2/3), width = 7)


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
ggsave("fig2cAbsolute.pdf", height = 7 * (2/3), width = 7)

############
############

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
ggsave("fig2cPercentInc.pdf", height = 7 * (2/3), width = 7)

###
###

ps5 <- read.csv("~/Downloads/fig2cCP.csv")
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
ggsave("fig2cCP.pdf", height = 7 * (2/3), width = 7)



##########################
##########################

# Figure 3a (_2bm_mesmer_vs_st_paper_v2.pdf)
res2 <- rbind(cbind(cohort='basel', bdc[[2]]), 
              cbind(cohort='meta', mdc[[2]]), 
              cbind(cohort='tonsil', tdc[[2]]))

res2[res2[,"larger_than30"] == 0, "larger_than30"] <- "1"
res2[res2[,"larger_than30"] == 2, "larger_than30"] <- ">1"

res2[,"Arcsinh"] <- ifelse(res2[,"Cofactor"] == 0, "No", "Yes")
res2[,"Cofactor"] <- ifelse(res2[,"Cofactor"] == 0, "N/A", res2[,"Cofactor"])

res2[res2[,"cohort"] == 'basel',"cohort"] <- 'Jackson 2020'
res2[res2[,"cohort"] == 'meta',"cohort"] <- 'Ali 2020'
res2[res2[,"cohort"] == 'tonsil',"cohort"] <- 'Tonsil'

res2[res2[,"Cofactor"] == "N/A","type"] <- "No Arcsinh" 
res2[res2[,"Cofactor"] == "1","type"] <- "Arcsinh: CF=1"
res2[res2[,"Cofactor"] == "5","type"] <- "Arcsinh: CF=5"

write.csv(res2, "/home/campbell/yulee/project/st/plot/final/fig3a.csv")

write.csv(filter(res2, Lambda == '1.0', Dist == 'Student-T', Cell == 'Yes', Relax == 'Yes'), 
          "/home/campbell/yulee/project/st/plot/final/fig3at.csv")

dim(df[df[,"larger_than30"] == '1',])
dim(df[df[,"larger_than30"] != '1',])

dim(df[df[,"larger_than30"] == '1' & df[,"cohort"] == 'Jackson 2020' ,])
dim(df[df[,"larger_than30"] != '1' & df[,"cohort"] == 'Jackson 2020' ,])

## create figure
df <- read.csv("fig3a.csv")
df %>%
  ggplot(aes(x = factor(larger_than30, level = c("1", ">1")), 
             y = as.numeric(doublet_prob),
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
  facet_grid(rows = glue("{type}") ~ glue("{cohort}")) + #glue("Dataset: {cohort}")
  theme_pubr(legend = 'none', base_size = 11)
ggsave("fig3a.pdf", width = 3.5, height = 7 * (2/3))


# Figure 3b (gamma_prob.R)
#rm(list = ls())

mdf <- read_csv("fig3b_meta_gamma.csv")
bdf <- read_csv("fig3b_basel_gamma.csv")
tdf <- read_csv("fig3b_tonsil_gamma.csv")

mdf[,"cohort"] <- 'Ali 2020'
bdf[,"cohort"] <- 'Jackson 2020'
tdf[,"cohort"] <- 'Tonsil'

df <- rbind(mdf, bdf, tdf)
write.csv(df, "fig3b.csv")

## create figure
#rm(list=ls())
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
  #filter(prob > quantile(prob, 0.98)) |> 
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
  #theme(legend.position = "none") + #, , panel.spacing.x=unit(0.5, "lines")
  scale_color_manual(values = cols)
ggsave("fig3b.pdf", width = 3.5, height = 7 * (2/3))

# Figure 3c (_2bm_mc_stat_paper.pdf)
res <- rbind(
  cbind(cohort='basel', MESMER='Yes', bdc[[1]]),
  cbind(cohort='meta', MESMER='Yes', mdc[[1]]),
  cbind(cohort='tonsil', MESMER='Yes', tdc[[1]]))

res[res2[,"larger_than30"] == 0, "larger_than30"] <- 1

res[,"Arcsinh"] <- ifelse(res[,"Cofactor"] == 0, "No", "Yes")
res[,"Cofactor"] <- ifelse(res[,"Cofactor"] == 0, "N/A", res[,"Cofactor"])

res[res[,"cohort"] == 'basel',"cohort"] <- 'Jackson 2020'
res[res[,"cohort"] == 'meta',"cohort"] <- 'Ali 2020'
res[res[,"cohort"] == 'tonsil',"cohort"] <- 'Tonsil'

res[res[,"Cofactor"] == "N/A","type"] <- "No Arcsinh" 
res[res[,"Cofactor"] == "1","type"] <- "Arcsinh: CF=1"
res[res[,"Cofactor"] == "5","type"] <- "Arcsinh: CF=5"

write.csv(filter(res, MESMER == 'Yes', Lambda == '1.0', Dist == 'Student-T', Cell == 'Yes', Relax == 'Yes'), 
          "/home/campbell/yulee/project/st/plot/final/fig3c.csv")

## create figure
#rm(list=ls())

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
#tmp <- filter(res2, Lambda == '1.0', Dist == 'Student-T', Cell == 'Yes', Relax == 'Yes', Initial == 'PG', Cofactor == 5, cohort == 'Jackson 2020')

# Figure 4de
#rm(list=ls())
minmax <- function(x, na.rm = TRUE) {
  return((x- min(x, na.rm = TRUE)) /(max(x, na.rm = TRUE)-min(x, na.rm = TRUE)))
}

mat1 <- read.csv("fig4d_centroids.csv")
mat2 <- read.csv("fig4e_centroids.csv")
samp_lb <- read.csv("fig4de_labels.csv")[,-1]

ct_labels <- c('CD20', 'Unsure', 'Unsure', 'CD68', 'Unsure', 
               'panCK', 'Unsure', 'CD20', 'panCK', 'CD3', 
               'CD20', 'CD68', 'CD3')

nc <- dim(mat1)[1]
ct_labels1 <- 1:nc

count <- NULL
for (c in 0:(nc-1)) {
  
  not_in_n_are_in  <- samp_lb[,"init_label"] != c & samp_lb[,"st_label_w_doublet_label"] == c ## werent in & are in
  were_in_n_are_in <- samp_lb[,"init_label"] == c & samp_lb[,"st_label_w_doublet_label"] == c ## were in & are in
  were_in_n_not_in <- samp_lb[,"init_label"] == c & samp_lb[,"st_label_w_doublet_label"] != c ## were in & not in
  
  count <- rbind(count, c(sum(not_in_n_are_in), sum(were_in_n_are_in), sum(were_in_n_not_in), 
                          sum(samp_lb[,"init_label"] == c),
                          sum(samp_lb[,"st_label_w_doublet_label"] == c)))
}

dmat <- apply(rbind(mat1[,-1], mat2[,-1]), 2, minmax)

i=0
ma1 <- t(as.matrix(dmat[(i*nc+1):((i+1)*nc),]))
colnames(ma1) <- ct_labels

f1 <- colorRamp2(seq(min(ma1), max(ma1), length = 3), c("blue", "#EEEEEE", "red"), space = "RGB")

haa1 <- HeatmapAnnotation(
  " " = anno_barplot(as.integer(count[,4]), axis = TRUE),
  'Cell Type' = ct_labels, ##rownames(tmp1),
  col = list('Cell Type' = c(
    'CD20' = '#dbdb8d', 'CD68' = '#279e68', 'CD3' = '#8c564b',
    'panCK' = "#999999", 'Unsure' = '#ffbb78')),
  annotation_name_side = "left",
  annotation_name_rot = 0, 
  annotation_name_gp = gpar(fontsize = 8), 
  annotation_legend_param = list('Cell Type' = list(title = "Cell Type", 
                                                    labels_gp = gpar(fontsize = 10))))

htt1 <- Heatmap(ma1, 
                name = 'Scaled\nExpression',
                column_title = " ",
                col = f1, #rev(rainbow(3)),
                top_annotation = haa1,
                column_split = colnames(ma1),
                column_labels = ct_labels1,
                cluster_rows = TRUE, cluster_columns = TRUE, show_heatmap_legend = TRUE, border = TRUE)

pdf('fig4d.pdf', height = 4.5)
draw(htt1, column_title = "Cluster", row_title = "Marker", column_title_side = 'bottom', merge_legend = TRUE)
dev.off()

ct_labels <- c('CD20', 'CD20', 'Unsure', 'CD68', 'CD3', 
               'panCK', 'panCK', 'CD20', 'panCK', 'CD3', 
               'CD20', 'CD68', 'CD3')

i=1
ma2 <- t(as.matrix(dmat[(i*nc+1):((i+1)*nc),]))
colnames(ma2) <- ct_labels

f2 <- colorRamp2(seq(min(ma2), max(ma2), length = 3), c("blue", "#EEEEEE", "red"), space = "RGB")

haa2 <- HeatmapAnnotation(
  " " = anno_barplot(as.integer(count[,5]), axis = TRUE),
  'Cell Type' = ct_labels, ##rownames(tmp1),
  col = list('Cell Type' = c(
    'CD20' = '#dbdb8d', 'CD68' = '#279e68', 'CD3' = '#8c564b',
    'panCK' = "#999999", 'Unsure' = '#ffbb78')),
  annotation_name_side = "left",
  annotation_name_rot = 0, 
  annotation_name_gp = gpar(fontsize = 8), 
  annotation_legend_param = list('Cell Type' = list(title = "Cell Type", 
                                                    labels_gp = gpar(fontsize = 10))))

htt2 <- Heatmap(ma2, 
                name = 'Scaled\nExpression',
                column_title = " ",
                col = f1, #rev(rainbow(3)),
                top_annotation = haa2,
                column_split = colnames(ma2),
                column_labels = ct_labels1,
                cluster_rows = TRUE, cluster_columns = TRUE, show_heatmap_legend = TRUE, border = TRUE)

pdf('fig4e.pdf', height = 4.5)
draw(htt2, column_title = "Cluster", row_title = "Marker", column_title_side = 'bottom', merge_legend = TRUE)
dev.off()

# Figure 5 (tonsil_st_ct_heatmap.pdf)
#rm(list=ls())

minmax <- function(x, na.rm = TRUE) {
  return((x- min(x, na.rm = TRUE)) /(max(x, na.rm = TRUE)-min(x, na.rm = TRUE)))
}

mat1 <- read.csv("fig5c_init_centroids.csv")[,-1]
mat2 <- read.csv("fig5c_st_centroids.csv")[,-1]
samp_lb <- read.csv("fig5c_labels.csv")[,-1]

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

pdf('fig5c_init.pdf', width = 10)
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

pdf('fig5c_st.pdf', width = 10)
draw(htt2, column_title = "Cluster", row_title = "Marker", column_title_side = 'bottom', merge_legend = TRUE)
dev.off()

###
###

rm(list=ls())
mat1 <- read.csv("mibi_in_centroids.csv")[,-1]
mat2 <- read.csv("mibi_st_centroids.csv")[,-1]
samp_lb <- read.csv("mibi_labels.csv")[,-1]
nc <- dim(mat1)[1]

count <- NULL
for (c in 0:(nc-1)) {
  count <- rbind(count, c(sum(samp_lb[,"init_label"] == c), sum(samp_lb[,"st_label"] == c)))
}

dmat <- apply(rbind(mat1, mat2), 2, minmax)

i=0
ma1 <- t(as.matrix(dmat[(i*nc+1):((i+1)*nc),]))
colnames(ma1) <- 1:nc

f1 <- colorRamp2(seq(min(ma1), max(ma1), length = 3), c("blue", "#EEEEEE", "red"), space = "RGB")

haa2 <- HeatmapAnnotation(
  " " = anno_barplot(as.integer(count[,2]), axis = TRUE),
  annotation_name_side = "left",
  annotation_name_rot = 0, 
  annotation_name_gp = gpar(fontsize = 8))

ct_labels1 <- 1:nc

htt1 <- Heatmap(ma1, 
                name = 'Scaled\nExpression',
                column_title = " ",
                col = f1, #rev(rainbow(3)),
                top_annotation = haa2,
                column_split = colnames(ma1),
                column_labels = ct_labels1,
                cluster_rows = TRUE, cluster_columns = TRUE, show_heatmap_legend = TRUE, border = TRUE)

pdf('mibi_init.pdf', width = 10)
draw(htt1, column_title = "Cluster", row_title = "Marker", column_title_side = 'bottom', merge_legend = TRUE)
dev.off()

i=1
ma2 <- t(as.matrix(dmat[(i*nc+1):((i+1)*nc),]))
colnames(ma2) <- 1:nc

f2 <- colorRamp2(seq(min(ma2), max(ma2), length = 3), c("blue", "#EEEEEE", "red"), space = "RGB")

column_ha <- HeatmapAnnotation(
  " " = anno_barplot(as.integer(count[,2]), axis = TRUE),
  annotation_name_side = "left",
  annotation_name_rot = 0, 
  annotation_name_gp = gpar(fontsize = 8))

htt2 <- Heatmap(ma2, 
                name = 'Scaled\nExpression',,
                column_title = " ",
                col = f2, #rev(rainbow(3)),
                top_annotation = column_ha,
                column_split = colnames(ma2),
                column_labels = ct_labels1,
                cluster_rows = TRUE, cluster_columns = TRUE, show_heatmap_legend = TRUE, border = TRUE)

pdf('mibi_st.pdf', width = 10)
draw(htt2, column_title = "Cluster", row_title = "Marker", column_title_side = 'bottom', merge_legend = TRUE)
dev.off()


# Figure 5d (tonsil_roi_prop_byLabel_m2_v2.pdf)

#rm(list=ls())
col_fun <- colorRamp2(c(0, 1), c("white", "blue"))
lgd <-  Legend(col_fun = col_fun, title = "Proportion", labels_gp = gpar(fontsize = 8))

samp_lb <- read.csv("fig5de_labels.csv")[,-1]

ct_labels <- c('Stromal', 'B', 'T', 'Immune', 'B',
               'B', 'T', 'B', 'Mix', 'Stromal',
               'Epithelial', 'B', 'T', 'Macrophage', 'Endothelial',
               'Stromal', 'B', 'Monocytes', 'B', 'Monocytes',
               'Neutrophils', 'B', 'T', 'B', 'Epithelial')

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
  
  for (c in 0:24) {
    
    #c = 0
    #tmp <- rbind(tmp, c(ROI, patient, ct_labels[c+1], 
    #                    c,
    #                    dim(filter(select_image, ST == c))[1],
    #                    dim(select_image)[1]-sum(select_image[,'ST'] == 'doublet')))
    
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

row_ha <- HeatmapAnnotation('Cell Type' = ct_labels,
                            col = list('Cell Type' = c(
                              'Stromal' = '#dbdb8d', 'B' = '#279e68', 'T' = '#8c564b',
                              'Immune' = "#999999", 'Mix' = '#ffbb78', 'Epithelial' = '#17becf',
                              'Macrophage' = '#aec7e8', 'Endothelial' = '#b5bd61', 
                              'Monocytes' = '#c5b0d5','Neutrophils' = '#f7b6d2')),
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

# Figure 5e (tonsil_roi_prop_byCelltype_m2_v2.pdf)
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
                              'Stromal' = '#dbdb8d', 'B' = '#279e68', 'T' = '#8c564b',
                              'Immune' = "#999999", 'Mix' = '#ffbb78', 'Epithelial' = '#17becf',
                              'Macrophage' = '#aec7e8', 'Endothelial' = '#b5bd61', 
                              'Monocytes' = '#c5b0d5','Neutrophils' = '#f7b6d2')),
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

# Figure 6ab (tonsil_nhood_m1_zsn.pdf)
#rm(list=ls())

p <- c(1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3)

nhood_zsd12 <- NULL
nhood_zsn12 <- NULL
for ( ii in 1:2 ) {
  
  zscores <- read.csv(paste0('enrichment_analysis/m', ii, '_zs_roi1.csv'))
  zscores[zscores == 'NaN' | is.na(zscores)] <- 0
  
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
  
  nhood_zsd <- select(nhood_zs, 'B-B', "T-T", "Endothelial-Endothelial", "Epithelial-Epithelial",
                      "Immune-Immune", "Macrophage-Macrophage", "Monocytes-Monocytes", "Neutrophils-Neutrophils", "Stromal-Stromal")
  
  nhood_zsn <- select(nhood_zs, -c('B-B', "T-T", "Endothelial-Endothelial", "Epithelial-Epithelial",
                                   "Immune-Immune", "Macrophage-Macrophage", "Monocytes-Monocytes", "Neutrophils-Neutrophils", "Stromal-Stromal"))
  
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
                                'Stromal' = '#dbdb8d', 'B' = '#279e68', 'T' = '#8c564b',
                                'Immune' = "#999999", 'Mix' = '#ffbb78', 'Epithelial' = '#17becf',
                                'Macrophage' = '#aec7e8', 'Endothelial' = '#b5bd61', 
                                'Monocytes' = '#c5b0d5','Neutrophils' = '#f7b6d2')),
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
                                'Stromal' = '#dbdb8d', 'B' = '#279e68', 'T' = '#8c564b',
                                'Immune' = "#999999", 'Mix' = '#ffbb78', 'Epithelial' = '#17becf',
                                'Macrophage' = '#aec7e8', 'Endothelial' = '#b5bd61', 
                                'Monocytes' = '#c5b0d5','Neutrophils' = '#f7b6d2'),
                                'Type 2'= c(
                                  'Stromal' = '#dbdb8d', 'B' = '#279e68', 'T' = '#8c564b',
                                  'Immune' = "#999999", 'Mix' = '#ffbb78', 'Epithelial' = '#17becf',
                                  'Macrophage' = '#aec7e8', 'Endothelial' = '#b5bd61', 
                                  'Monocytes' = '#c5b0d5','Neutrophils' = '#f7b6d2')),
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
  labs(x = "Enrichment assigning all cells", y = "Enrichment removing poorly\nsegmented cells") +
  theme_pubr(base_size = 8) #base_size = 6
ggsave("fig6b_tonsil_nhood_m12_diagonal_zs_comparisons.pdf", height = 2.06, width = 2.32)

ggplot(as.data.frame(nhood_zsn12), aes(x=method1, y=method2)) + geom_point(alpha=0.5) + 
 labs(x = "Enrichment assigning all cells", y = "Enrichment removing poorly\nsegmented cells") +
  theme_pubr(base_size = 8) #
ggsave("fig6b_tonsil_nhood_m12_offdiagonal_zs_comparisons.pdf", height = 2.06, width = 2.32)