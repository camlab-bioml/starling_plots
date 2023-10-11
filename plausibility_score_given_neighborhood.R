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

## neighborhood analysis
neighborhood_analysis <- function(cohort, seg_type) {
  
  path <- paste0("/home/campbell/yulee/project/result/", cohort, "/_1nh/", seg_type, "/")
  fns <- list.files(path)
  
  n <- fns[grep("^n", fns)]
  m <- fns[grep("^m", fns)]
  nm <- c(n, m)
  wfs <- nm[grep("s10000", nm)]
  
  res <- NULL
  for (i in 1:length(wfs)) {
    fn <- strsplit(wfs[i], "_")[[1]]
    tmp <- read.csv(paste0(path, wfs[i]))
    tmp[,"X"] <- toupper(fn[1])
    tmp[,"n"] <- substr(fn[4], 2, 100)
    tmp[,"r"] <- strsplit(fn[5], "r|.csv")[[1]][2]
    tmp[,'cohort'] <- cohort
    tmp[,'seg'] <- seg_type
    res <- rbind(res, tmp)
  }
  res
}

## Basel (Jackson et al. 2020) and Metabrics (Ali et al. 2020) cohorts
## segmentation is given by these publications
bcp <- neighborhood_analysis('basel', 'cp')
mcp <- neighborhood_analysis('meta', 'cp')

## Deepcell segmentation was applied on Basel and Metabrics cohorts
bdc <- neighborhood_analysis('basel', 'dc')
mdc <- neighborhood_analysis('meta', 'dc')

## Cells from 16 unpublished Tonsil ROIs were retraived via Deepcell segmentation
tdc <- neighborhood_analysis('tonsil', 'dc')

## Combined results and make plots
res <- rbind(bcp, bdc, mcp, mdc, tdc)

res[,"Arcsinh"] <- ifelse(res[,"Cofactor"] == 0, "No", "Yes")
res[,"Cofactor"] <- ifelse(res[,"Cofactor"] == 0, "N/A", res[,"Cofactor"])
res[res[,"cohort"] == 'basel',"cohort"] <- 'Jackson 2020'
res[res[,"cohort"] == 'meta',"cohort"] <- 'Ali 2020'
res[res[,"cohort"] == 'tonsil',"cohort"] <- 'Tonsil'

## Overall by cohort
as.data.frame(filter(res, X != "M")) %>%
  ggplot(aes(x = X, y = Score, fill = X, color = X)) +
  facet_grid(rows=glue("Arcsinh: {Arcsinh}")+glue("CF: {Cofactor}")~glue("Dataset: {cohort}")) + 
  geom_violin(trim = F) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  scale_color_manual(values = c("#E69F00", "#56B4E9")) +
  geom_boxplot(width = 0.1, fill = "white", color = "grey40") +
  labs(x = 'N (Cell has neighbours) \n NN (Cell has no neighbours)', y = 'Plausibility Score') + 
  theme_pubr(legend = 'none') +
  theme(text=element_text(size=15)) +
  scale_y_continuous(breaks = c(0,0.5,1))
ggsave(paste0("/home/campbell/yulee/project/result/plots/n_vs_nn_ds.pdf"), width=10, height=10)

## Overall by initialization method
as.data.frame(filter(res, X != "M")) %>% 
  ggplot(aes(x = X, y = Score, fill = X, color = X)) + 
  facet_grid(rows=glue("NS: {type}")~glue("Initial: {Initial}")) + 
  geom_violin(trim = F) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  scale_color_manual(values = c("#E69F00", "#56B4E9")) +
  geom_boxplot(width = 0.1, fill = "white", color = "grey40") +
  labs(x = 'N (Cell has neighbours) \n NN (Cell has no neighbours)', y = 'Plausibility Score') + 
  theme_pubr(legend = 'none') +
  theme(text=element_text(size=15)) +
  scale_y_continuous(breaks = c(0,0.5,1))
ggsave(paste0("/home/campbell/yulee/project/result/plots/n_vs_nn_init.pdf"), width=10, height=10)

