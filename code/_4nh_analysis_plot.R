rm(list=ls())

library(tidyr)
library(dplyr)
library(ggplot2)
library(glue)
library(stringr)
library(ggpubr)
library(ggsci)

neighborhood_analysis <- function(cohort, seg_type) {
  
  #cohort = 'meta'
  #seg_type = 'cp'
  path <- paste0("/home/campbell/yulee/project/st/", cohort, "/analysis/_1nh/", seg_type, "/")
  fns <- list.files(path)
  
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

bcp <- neighborhood_analysis('basel', 'cp')
bdc <- neighborhood_analysis('basel', 'dc')
mdc <- neighborhood_analysis('meta', 'dc')
tdc <- neighborhood_analysis('tonsil', 'dc')

res <- rbind(cbind(cohort='basel', MESMER='No', bcp), cbind(cohort='basel', MESMER='Yes', bdc), 
             cbind(cohort='meta', MESMER='Yes', mdc), cbind(cohort='tonsil', MESMER='Yes',tdc))

res[,"Arcsinh"] <- ifelse(res[,"Cofactor"] == 0, "No", "Yes")
res[,"Cofactor"] <- ifelse(res[,"Cofactor"] == 0, "N/A", res[,"Cofactor"])
res[res[,"cohort"] == 'basel',"cohort"] <- 'Jackson 2020'
res[res[,"cohort"] == 'meta',"cohort"] <- 'Ali 2020'
res[res[,"cohort"] == 'tonsil',"cohort"] <- 'Tonsil'

res[res[,"Cofactor"] == "N/A","type"] <- "None" 
res[res[,"Cofactor"] == "1","type"] <- "Arcsinh (CF=1)"
res[res[,"Cofactor"] == "5","type"] <- "Arcsinh (CF=5)"

as.data.frame(res) %>%
  ggplot(aes(x = X, y = Score, fill = X, color = X)) +
  #facet_grid(rows=glue("Transformation: {type}")~glue("Dataset: {cohort}")+glue("MESMER: {MESMER}")+glue("Initial: {Initial}")) +
  facet_grid(rows=glue("Transformation: {type}")~glue("Dataset: {cohort}")+glue("Initial: {Initial}")) +
  #facet_grid(rows=glue("Dataset: {cohort}")+glue("Initial: {Initial}")~glue("Transformation: {type}")) + 
  geom_hline(yintercept=.2, linetype="dashed") +
  geom_hline(yintercept=.4, linetype="dashed") +
  geom_hline(yintercept=.6, linetype="dashed") +
  geom_violin(trim = F) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  scale_color_manual(values = c("#E69F00", "#56B4E9")) +
  geom_boxplot(width = 0.1, fill = "white", color = "grey40") +
  labs(x = 'N (Cell has neighbours) \n NN (Cell has no neighbours)', y = 'Plausibility Score') + 
  theme_pubr(legend = 'none') +
  theme(text=element_text(size=15)) +
  scale_y_continuous(breaks = c(0,0.5,1))
ggsave(paste0("/home/campbell/yulee/project/st/plot/_1nh_overall.pdf"), width=25, height=9)

as.data.frame(res) %>%
  ggplot(aes(x = X, y = Score, fill = X, color = X)) +
  #facet_grid(rows=glue("Transformation: {type}")~glue("Dataset: {cohort}")+glue("MESMER: {MESMER}")+glue("Initial: {Initial}")) +
  #facet_grid(rows=glue("Transformation: {type}")~glue("Dataset: {cohort}")+glue("Initial: {Initial}")) +
  facet_grid(rows=glue("Dataset: {cohort}")+glue("Initial: {Initial}")~glue("Transformation: {type}")) + 
  geom_hline(yintercept=.2, linetype="dashed") +
  geom_hline(yintercept=.3, linetype="dashed") +
  geom_hline(yintercept=.4, linetype="dashed") +
  geom_violin(trim = F) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  scale_color_manual(values = c("#E69F00", "#56B4E9")) +
  geom_boxplot(width = 0.1, fill = "white", color = "grey40") +
  labs(x = 'N (Cell has neighbours) \n NN (Cell has no neighbours)', y = 'Plausibility Score') + 
  theme_pubr(legend = 'none') +
  theme(text=element_text(size=15)) +
  scale_y_continuous(breaks = c(0,0.5,1))
ggsave(paste0("/home/campbell/yulee/project/st/plot/_1nh_overall_t.pdf"), width=9, height=25)