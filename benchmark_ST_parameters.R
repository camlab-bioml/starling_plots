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
ps1 = read.csv(paste0('/home/campbell/yulee/project/result/basel/_2bm/dc/s10k_t3.csv'))[,-1]
ps2 = read.csv(paste0('/home/campbell/yulee/project/result/meta/_2bm/cp/s10k_t3.csv'))[,-1]
ps3 = read.csv(paste0('/home/campbell/yulee/project/result/tonsil/_2bm/cp/s10k_t3.csv'))[,-1]
ps <- rbind(ps1, ps2, ps3)

#ps = read.csv('/home/campbell/yulee/project/result/tonsil/_2bm/dc/misc/nm2_t3.csv')[,-1]
colnames(ps) <- c("s", 'Initial', "n", "t", "Run", "Noise", "Lambda", "Cell", "Relax", "Before", "After", "rf", "rf_F1", "st_F1", 'Type')

df11 <- group_by(ps, Initial, Lambda, Cell, Relax, Run, Noise, Type) %>% 
  summarise('At initialization'=mean(Before), 'After training'=mean(After), 'Two step'=mean(rf), .groups = 'drop') %>%
  pivot_longer(cols = 8:10, names_to = 'Metric', values_to = 'score') %>%
  group_by(Initial, Lambda, Cell, Relax, Noise, Type, Metric) %>% summarise(Correlation=mean(score), .groups = 'drop')

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

p3 + facet_grid(glue("{Cell}")+glue("Type: {Type}")~glue("'Cell overlap': {Relax}")+glue("'Model Noise': {Noise}"), labeller = label_parsed) +
  #p3 + facet_grid(glue("AC: {AC}")+glue("CF: {Cofactor}")~glue("'Model cell size': {Cell}")+glue("'Model cell overlap': {Relax}")+glue("'Distribution': {Dist}"), labeller = label_parsed) + 
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
#ggsave(paste0("/home/campbell/yulee/project/starling/analysis/exp6/", cohort, "_dc_ps5.png"), width = 12)

p3 + facet_grid(glue("'Model cell size': {Cell}")~glue("'Model cell overlap': {Relax}"), labeller = label_parsed) +
  #p3 + facet_grid(glue("AC: {AC}")+glue("CF: {Cofactor}")~glue("'Model cell size': {Cell}")+glue("'Model cell overlap': {Relax}")+glue("'Distribution': {Dist}"), labeller = label_parsed) + 
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


#ps2[,'type'] <- 2
#ps3[,'type'] <- 3
#ps <- rbind(ps2, ps3)
#ps2 = read.csv(paste0('/home/campbell/yulee/project/result/', cohort, '/_2bm/', seg_type, '/s250k_t2.csv'))[,-1]
#ps3 = read.csv(paste0('/home/campbell/yulee/project/result/', cohort, '/_2bm/', seg_type, '/s250k_t3.csv'))[,-1]

#ps = read.csv(paste0('/home/campbell/yulee/project/result/', cohort, '/_2bm/', seg_type, '/s250k_t2.csv'))[,-1]
#ps = read.csv(paste0('/home/campbell/yulee/project/result/', cohort, '/_2bm/', seg_type, '/s250k_t3.csv'))[,-1]
#ps = read.csv(paste0('/home/campbell/yulee/project/result/', cohort, '/_2bm/', seg_type, '/s10k_t', tr_type, '.csv'))[,-1]
