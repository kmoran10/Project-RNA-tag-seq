


library(lme4)
library(lmerTest)
library(tidyverse)
library(ggpubr)
library(car)
library(rstatix)



bnp <- read.csv("raw_data/id.csv")


#P42 weight - save as w400 h600

ggplot(bnp, aes(group, wt.day.15, fill=group)) +
  stat_compare_means(method = "t.test", size = 6, label.x = 0.8) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(size = 4, position=position_jitter(0.15)) + 
  scale_fill_manual(values=c("white", "darkgray")) +
  labs(title="A. P42 Weight",x="Group", y = "Grams") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) 

t.test(wt.day.15 ~ group, bnp)
bnp %>% group_by(group) %>% summarise_at(vars(wt.day.15), list(name=sd))
bnp %>% cohens_d(wt.day.15 ~ group)
# [Control: 83.8 ± 7.66; Stressed: 91.8 ± 11.9; t(15.35) = 1.78, p < 0.094, d = 0.80]




#rate of gain - save as w400 h600

ggplot(bnp, aes(group, wt.rate.pct, fill=group)) +
  stat_compare_means(method = "t.test", size = 6, label.x = 0.8) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(size = 4, position=position_jitter(0.15)) + 
  scale_fill_manual(values=c("white", "darkgray")) +
  labs(title="B. Rate of Weight Gain",x="Group", y = "% of P28 Weight Gained") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) 

t.test(wt.rate.pct ~ group, bnp)
bnp %>% group_by(group) %>% summarise_at(vars(wt.rate.pct), list(name=sd))
bnp %>% cohens_d(wt.rate.pct ~ group)
# [Control: 48.89 ± 7.08; Stressed: 59.38 ± 7.01; t(18.00) = 3.33, p < 0.01, d = 1.49]




#P42 fat - save as w400 h600

ggplot(bnp, aes(group, total.fat, fill=group)) +
  stat_compare_means(method = "t.test", size = 6, label.x = 0.8) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(size = 4, position=position_jitter(0.15)) + 
  scale_fill_manual(values=c("white", "darkgray")) +
  labs(title="C. Fat Mass",x="Group", y = "Grams") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) 

t.test(total.fat ~ group, bnp)
bnp %>% group_by(group) %>% summarise_at(vars(total.fat), list(name=sd))
bnp %>% cohens_d(total.fat ~ group)
# [Control: 5.21 ± 1.13; Stressed: 6.83 ± 2.25; t(13.27) = 2.04, p = 0.062, d = .91]





#food intake - save as w400 h600 --- p=.5 - no total food diff

ggplot(bnp, aes(group, food.total, fill=group)) +
  stat_compare_means(method = "t.test", size = 6) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(size = 4, position=position_jitter(0.15)) + 
  scale_fill_manual(values=c("white", "darkgray")) +
  labs(title="Total Food Intake",x="Group", y = "Grams") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) 

t.test(food.total ~ group, bnp)
bnp %>% group_by(group) %>% summarise_at(vars(food.total), list(name=sd))
bnp %>% cohens_d(food.total ~ group)




#food efficiency - save as w400 h600

ggplot(bnp, aes(group, food.eff, fill=group)) +
  stat_compare_means(method = "t.test", size = 6, label.x = 0.8) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(size = 4, position=position_jitter(0.15)) + 
  scale_fill_manual(values=c("white", "darkgray")) +
  labs(title="D. Food Efficiency",x="Group", y = "Grams Gained/Grams Eaten") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) 

t.test(food.eff ~ group, bnp)
bnp %>% group_by(group) %>% summarise_at(vars(food.eff), list(name=sd))
bnp %>% cohens_d(food.eff ~ group)
# [Control: 0.162 ± 0.025; Stressed: 0.191 ± 0.028; t(17.75) = 2.42, p < 0.05, d = 1.09]





