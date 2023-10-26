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

## Benchmark Starling

## plausibility score
rm(list=ls())
for (cohort in c("basel", "meta", 'tonsil')) {
  
  print(cohort)
  
  if ( cohort == 'tonsil' ) {
    correlation <- read.csv(paste0('/home/campbell/yulee/project/starling/analysis/exp6/tonsil_dc_ps.csv'))[,-1]
  } else {
    correlation <- read.csv(paste0('/home/campbell/yulee/project/starling/analysis/exp5/', cohort, '_cp_ps.csv'))[,-1]
  }
  
  correlation[,'Cluster'] <- as.factor(correlation[,'nc'])
  correlation[,'Initial'] <- as.factor(correlation[,'ia'])
  correlation[,'cf'][correlation[,'cf'] == '0'] <- 'N/A'
  correlation[,'Cofactor'] <- as.factor(correlation[,'cf'])
  correlation[,'Dist'] <- as.factor(correlation[,'nm'])
  correlation[,'Lambda'] <- as.factor(correlation[,'lv'])
  correlation[,'Cell'] <- as.factor(correlation[,'cs'])
  correlation[,'Relax'] <- as.factor(correlation[,'rr'])
  correlation[,'Run'] <- as.factor(correlation[,'r'])
  correlation[,'Cell'] <- ifelse(correlation[,'Cell'] == 1, 'Yes', 'No')
  correlation[,'Relax'] <- ifelse(correlation[,'Relax'] == 1, 'Yes', 'No')
  correlation[,'Dist'] <- ifelse(correlation[,'Dist'] == 0, 'Normal', 'Student-T')
  correlation[,'AC'] <- ifelse(correlation[,'Cofactor'] == 'N/A', 'No', 'Yes')
  
  df11 <- correlation %>% group_by(Cluster, Initial, Dist, AC, Cofactor, Lambda, Cell, Relax, Run) %>% 
    summarise('At initialization'=mean(metric), 'After training'=mean(after), 'Two step'=mean(rf), .groups = 'drop') %>%
    pivot_longer(cols = 10:12, names_to = 'Metric', values_to = 'score') %>%
    group_by(Initial, Dist, AC, Cofactor, Lambda, Cell, Relax, Metric) %>% summarise(Correlation=mean(score), .groups = 'drop')
  
  df1 <- df11 %>% filter(Metric == "At initialization") 
  df2 <- df11 %>% filter(Metric == "After training")
  df3 <- df11 %>% filter(Metric == "Two step")
  
  p1 <- ggplot(df2, aes(x = Initial, y = Correlation, shape=Metric)) +
    #geom_point() +
    stat_summary(fun=mean, geom="line", aes(group = Lambda, col = Lambda)) +
    stat_summary(fun=mean, geom="point",  aes(group = Lambda, col = Lambda))
  
  p2 <- p1 + 
    stat_summary(data=df1, fun=mean, geom="line", aes(group = 1), linetype = "dotted") + 
    stat_summary(data=df1, fun=mean, geom="point")
  
  p3 <- p2 +
    stat_summary(data=df3, fun=mean, geom="line", aes(group = 1), linetype = "dashed") + 
    stat_summary(data=df3, fun=mean, geom="point")
  
  if (cohort == 'tonsil') {
    p3 + facet_grid(glue("AC: {AC}")+glue("CF: {Cofactor}")~glue("'Model cell overlap': {Relax}")+glue("'Distribution': {Dist}"), labeller = label_parsed) + 
      labs(x = "Initialization Algorithm", y = "Plausibility Score") +
      #ylim(8, 16) + 
      scale_colour_brewer(palette = "Set2") +
      cowplot::theme_cowplot() +
      theme(legend.position="bottom", 
            #legend.text=element_text(size=8),
            #legend.title = element_text(size=10), #element_blank(),
            #axis.text.x = element_text(size = 10),
            strip.background = element_rect(fill='white'), 
            strip.text = element_text(face='bold')) + 
      theme_pubr()
    
  } else {
    p3 + facet_grid(glue("AC: {AC}")+glue("CF: {Cofactor}")~glue("'Model cell size': {Cell}")+glue("'Model cell overlap': {Relax}")+glue("'Distribution': {Dist}"), labeller = label_parsed) + 
      labs(x = "Initialization Algorithm", y = "Plausibility Score") +
      #ylim(8, 16) + 
      scale_colour_brewer(palette = "Set2") +
      cowplot::theme_cowplot() +
      theme(legend.position="bottom", 
            #legend.text=element_text(size=8),
            #legend.title = element_text(size=10), #element_blank(),
            #axis.text.x = element_text(size = 10),
            strip.background = element_rect(fill='white'), 
            strip.text = element_text(face='bold')) + 
      theme_pubr()
  }
  ggsave(paste0("/home/campbell/yulee/project/result/plots/", cohort, "_benchmark_plausibility.pdf"), width = 12)
}

###
###

## f1 score
rm(list=ls())
for (cohort in c("basel", "meta", "tonsil")) {
  
  print(cohort)
  
  if ( cohort == 'tonsil' ) {
    doublet_perf <- read.csv(paste0('/home/campbell/yulee/project/starling/analysis/exp6/tonsil_dc_dp.csv'))[,-1]
  } else {
    doublet_perf <- read.csv(paste0('/home/campbell/yulee/project/starling/analysis/exp5/', cohort, '_cp_dp.csv'))[,-1]
  }
  
  doublet_perf[,'Cluster'] <- as.factor(doublet_perf[,'nc'])
  doublet_perf[,'Initial'] <- as.factor(doublet_perf[,'ia'])
  doublet_perf[,'cf'][doublet_perf[,'cf'] == '0'] <- 'N/A'
  doublet_perf[,'Cofactor'] <- as.factor(doublet_perf[,'cf'])
  doublet_perf[,'Dist'] <- as.factor(doublet_perf[,'nm'])
  doublet_perf[,'Lambda'] <- as.factor(doublet_perf[,'lv'])
  doublet_perf[,'Cell'] <- as.factor(doublet_perf[,'cs'])
  doublet_perf[,'Relax'] <- as.factor(doublet_perf[,'rr'])
  doublet_perf[,'Run'] <- as.factor(doublet_perf[,'r'])
  doublet_perf[,'Metric'] <- doublet_perf[,'metric']
  doublet_perf[,'Metric'] <- as.character(doublet_perf[,'Metric'])
  doublet_perf[,'Metric'][doublet_perf[,'Metric'] == 'rf_F1'] = "Random Forest"
  doublet_perf[,'Metric'][doublet_perf[,'Metric'] == 'st_F1'] = "Starling"
  doublet_perf[,'Cell'] <- ifelse(doublet_perf[,'Cell'] == 1, 'Yes', 'No')
  doublet_perf[,'Relax'] <- ifelse(doublet_perf[,'Relax'] == 1, 'Yes', 'No')
  doublet_perf[,'Dist'] <- ifelse(doublet_perf[,'Dist'] == 0, 'Normal', 'Student-T')
  doublet_perf[,'AC'] <- ifelse(doublet_perf[,'Cofactor'] == 'N/A', 'No', 'Yes')
  
  df21 <- doublet_perf %>% group_by(Cluster, Initial, AC, Dist, Cofactor, Lambda, Cell, Relax, Metric) %>% summarise(F1=mean(F1), .groups = 'drop')
  
  df1 <- df21 %>% filter(Metric == "Starling")
  df2 <- df21 %>% filter(Metric == "Random Forest") 
  
  p1 <- ggplot(df1, aes(x = Initial, y = F1, shape=Metric)) +
    stat_summary(fun=mean, geom="line", aes(group = Lambda, col = Lambda)) +
    stat_summary(fun=mean, geom="point",  aes(group = Lambda, col = Lambda))
  
  p2 <- p1 + 
    stat_summary(data=df2, fun=mean, geom="line", linetype = 2, aes(group = 1)) +
    stat_summary(data=df2, fun=mean, geom="point")
  
  if ( cohort == "tonsil" ) {
    p2 + facet_grid(glue("AC: {AC}")+glue("CF: {Cofactor}")~glue("'Model cell overlap': {Relax}")+glue("'Distribution': {Dist}"), labeller = label_parsed) + 
      labs(x = "Initialization Algorithm", y = "F1") +
      #ylim(8, 16) + 
      scale_colour_brewer(palette = "Set2") +
      cowplot::theme_cowplot() +
      theme(legend.position="bottom", 
            strip.background = element_rect(fill='white'), 
            strip.text = element_text(face='bold')) + 
      theme_pubr()
  } else {
    p2 + facet_grid(glue("AC: {AC}")+glue("CF: {Cofactor}")~glue("'Model cell size': {Cell}")+glue("'Model cell overlap': {Relax}")+glue("'Distribution': {Dist}"), labeller = label_parsed) + 
      labs(x = "Initialization Algorithm", y = "F1") +
      #ylim(8, 16) + 
      scale_colour_brewer(palette = "Set2") +
      cowplot::theme_cowplot() +
      theme(legend.position="bottom", 
            strip.background = element_rect(fill='white'), 
            strip.text = element_text(face='bold')) + 
      theme_pubr()
  }
  ggsave(paste0("/home/campbell/yulee/project/result/plots/", cohort, "_benchmark_f1.pdf"), width = 12)
}

###
### 

## morphological score

rm(list=ls())
for (cohort in c("basel", "meta", "tonsil")) {
  
  #cohort = 'basel'
  print(cohort)
  
  if ( cohort == 'tonsil' ) {
    ms <- read.csv(paste0('/home/campbell/yulee/project/starling/analysis/exp6/tonsil_dc_mc.csv'))[,-1]
  } else {
    ms <- read.csv(paste0('/home/campbell/yulee/project/starling/analysis/exp5/', cohort, '_cp_mc.csv'))[,-1]
  }
  
  ms[,'Cluster'] <- as.factor(ms[,'nc'])
  ms[,'Initial'] <- as.factor(ms[,'ia'])
  ms[,'cf'][ms[,'cf'] == '0'] <- 'N/A'
  ms[,'Cofactor'] <- as.factor(ms[,'cf'])
  ms[,'Dist'] <- as.factor(ms[,'nm'])
  ms[,'Lambda'] <- as.factor(ms[,'lv'])
  ms[,'Cell'] <- as.factor(ms[,'cs'])
  ms[,'Relax'] <- as.factor(ms[,'rr'])
  ms[,'Run'] <- as.factor(ms[,'r'])
  #ms[,'Metric'] <- sp[,'metric']
  #ms[,'Metric'] <- as.character(ms[,'Metric'])
  #ms[,'Metric'][ms[,'Metric'] == 'rf_F1'] = "Random Forest"
  #ms[,'Metric'][ms[,'Metric'] == 'st_F1'] = "Starling"
  ms[,'Cell'] <- ifelse(ms[,'Cell'] == 1, 'Yes', 'No')
  ms[,'Relax'] <- ifelse(ms[,'Relax'] == 1, 'Yes', 'No')
  ms[,'Dist'] <- ifelse(ms[,'Dist'] == 0, 'Normal', 'Student-T')
  ms[,'Initial'] <- substr(ms[,'Initial'], 6, 100)
  ms[,'AC'] <- ifelse(ms[,'Cofactor'] == 'N/A', 'No', 'Yes')
  
  if ( cohort == 'tonsil' ) {
    ms %>% ggplot( aes(x = Initial, y = Size, fill = Initial, color = Initial)) +
      geom_hline(yintercept=0, linetype="dashed") +
      geom_violin(trim = F) +
      scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9", "#009E73")) +
      scale_color_manual(values = c("#999999", "#E69F00", "#56B4E9", "#009E73")) +
      geom_boxplot(width = 0.1, fill = "white", color = "grey40") + 
      labs(x = 'Initialization Algorithm', y = 'Estimate (t-test)') + 
      facet_grid(glue("AC: {AC}")+glue("CF: {Cofactor}")+glue("lambda: {Lambda}")~glue("'Model cell overlap': {Relax}")+glue("'Distribution': {Dist}"), labeller = label_parsed) +
      theme_pubr(legend = 'none')
    ggsave(paste0("/home/campbell/yulee/project/result/plots/", cohort, "_benchmark_morphological_est.pdf"), width = 12, height = 8)
    
    ms %>% ggplot( aes(x = Initial, y = -log(Pval), fill = Initial, color = Initial)) +
      geom_violin(trim = F) +
      scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9", "#009E73")) +
      scale_color_manual(values = c("#999999", "#E69F00", "#56B4E9", "#009E73")) +
      geom_boxplot(width = 0.1, fill = "white", color = "grey40") + 
      labs(x = 'Initialization Algorithm', y = 'P-Value (-log scale)') + 
      #facet_grid(Cofactor~Dist+Cell)
      facet_grid(glue("AC: {AC}")+glue("CF: {Cofactor}")+glue("lambda: {Lambda}")~glue("'Model cell overlap': {Relax}")+glue("'Distribution': {Dist}"), labeller = label_parsed) +
      theme_pubr(legend = 'none')
    ggsave(paste0("/home/campbell/yulee/project/result/plots/", cohort, "_benchmark_morphological_pva.pdf"), width = 12, height = 8)
    
  } else {
    ms %>% ggplot(aes(x = Initial, y = Size, fill = Initial, color = Initial)) +
      geom_hline(yintercept=0, linetype="dashed") + 
      geom_violin(trim = F) +
      scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9", "#009E73")) +
      scale_color_manual(values = c("#999999", "#E69F00", "#56B4E9", "#009E73")) +
      geom_boxplot(width = 0.1, fill = "white", color = "grey40") + 
      labs(x = 'Initialization Algorithm', y = 'Estimate (t-test)') + 
      facet_grid(glue("AC: {AC}")+glue("CF: {Cofactor}")+glue("lambda: {Lambda}")~glue("'Model cell size': {Cell}")+glue("'Model cell overlap': {Relax}")+glue("'Distribution': {Dist}"), labeller = label_parsed) +
      theme_pubr(legend = 'none')
    ggsave(paste0("/home/campbell/yulee/project/result/plots/", cohort, "_benchmark_morphological_est.pdf"), width = 12, height = 8)
    
    ms %>% ggplot(aes(x = Initial, y = -log(Pval), fill = Initial, color = Initial)) +
      geom_violin(trim = F) +
      scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9", "#009E73")) +
      scale_color_manual(values = c("#999999", "#E69F00", "#56B4E9", "#009E73")) +
      geom_boxplot(width = 0.1, fill = "white", color = "grey40") + 
      labs(x = 'Initialization Algorithm', y = 'P-Value (-log scale)') + 
      #facet_grid(Cofactor~Dist+Cell)
      facet_grid(glue("AC: {AC}")+glue("CF: {Cofactor}")+glue("lambda: {Lambda}")~glue("'Model cell size': {Cell}")+glue("'Model cell overlap': {Relax}")+glue("'Distribution': {Dist}"), labeller = label_parsed) +
      theme_pubr(legend = 'none')
    ggsave(paste0("/home/campbell/yulee/project/result/plots/", cohort, "_benchmark_morphological_pva.pdf"), width = 12, height = 8)
    
  }
}

