### --- sauger_hiphop_plot.R --- ###
### created by SPJ on 050925 ###
# purpose: plot the hiphop test f1 results as well as 


getwd()
setwd("/Users/samjohnson/Desktop/hiphop")

#install.packages("hiphop")
library(hiphop)
library(dplyr)
library(tidyverse)
library(purrr)

# misleading title here. only 100, 50, and 25. goal is to first reconcile wtf is going on with my file names.
top1_100_50_25 <- read.csv("top1_100_75_50_25.csv")
top1_100_75_50 <- read.csv("top1_100_75_50.csv")

# see which categories are in each df
colnames(top1_100_50_25)
colnames(top1_100_75_50)

# need to add 25 to top1_100_75_50.
# isolate the 25% information
top1_25 <- top1_100_50_25[ , c(1:2, 15:21)]

# merge the dfs. same information now for both ptps situations. 
top1_100_75_50_25 <- merge(top1_100_75_50, top1_25, by = "offspring")

# the merge caused the brood column to be repeated. this removes the duplicated one. 
top1_100_75_50_25 <- top1_100_75_50_25 %>% 
    rename(brood = brood.x)
top1_100_75_50_25 <- top1_100_75_50_25[ , c(1:21, 23:29)]

# THIS IS THE CORRECT ONE. 
write.csv(top1_100_75_50_25, "top1_100_75_50_25.csv")
      
to_plot <- read.csv("top1_100_75_50_25.csv", row.names = 1)  
top_f1 <- read.csv("top_f1.csv", header = TRUE, row.names = 1)

to_plot <- data.frame(to_plot, top_f1$hothiphop.parents)

# now we need to filter this to get each of the     
ptps100 <- to_plot$score_100

ptps75_0 <- to_plot %>% 
    filter(count_par_75 == 0)
ptps75_par0 <- ptps75_0$score_75

ptps75_1 <- to_plot %>% 
    filter(count_par_75 == 1)
ptps75_par1 <- ptps75_1$score_75

ptps75_2 <- to_plot %>% 
    filter(count_par_75 == 2)
ptps75_par2 <- ptps75_2$score_75


ptps50_0 <- to_plot %>% 
    filter(count_par_50 == 0)
ptps50_par0 <- ptps50_0$score_50

ptps50_1 <- to_plot %>% 
  filter(count_par_50 == 1)
ptps50_par1 <- ptps50_1$score_50

ptps50_2 <- to_plot %>% 
  filter(count_par_50 == 2)
ptps50_par2 <- ptps50_2$score_50

ptps25_0 <- to_plot %>% 
    filter(count_par_25 == 0)
ptps25_par0 <- ptps25_0$score_25

ptps25_1 <- to_plot %>% 
    filter(count_par_25 == 1)
ptps25_par1 <- ptps25_1$score_25

ptps25_2 <- to_plot %>% 
  filter(count_par_25 == 2)
ptps25_par2 <- ptps25_2$score_25

allf1 <- top_f1$hothiphop.parents

full_to_plot <- bind_rows(
  data.frame(score = ptps100, group = "TestF1-PTPS100"),
  data.frame(score = ptps75_par0, group = "TestF1-PTPS75-Par0"),
  data.frame(score = ptps75_par1, group = "TestF1-PTPS75-Par1"),
  data.frame(score = ptps75_par2, group = "TestF1-PTPS75-Par2"),
  data.frame(score = ptps50_par0, group = "TestF1-PTPS50-Par0"),
  data.frame(score = ptps50_par1, group = "TestF1-PTPS50-Par1"),
  data.frame(score = ptps50_par2, group = "TestF1-PTPS50-Par2"),
  data.frame(score = ptps25_par0, group = "TestF1-PTPS25-Par0"),
  data.frame(score = ptps25_par1, group = "TestF1-PTPS25-Par1"),
  data.frame(score = ptps25_par2, group = "TestF1-PTPS25-Par2"),
  data.frame(score = allf1, group = "WildF1")
)

full_to_plot$group <- factor(full_to_plot$group, 
                             levels = c("TestF1-PTPS100", 
                                        "TestF1-PTPS75-Par0", "TestF1-PTPS75-Par1", "TestF1-PTPS75-Par2",
                                        "TestF1-PTPS50-Par0", "TestF1-PTPS50-Par1", "TestF1-PTPS50-Par2",
                                        "TestF1-PTPS25-Par0", "TestF1-PTPS25-Par1", "TestF1-PTPS25-Par2",
                                        "WildF1"))

full_to_plot <- full_to_plot %>%
  mutate(
    legend_group = case_when(
      grepl("Par0", group) ~ "0 True Parents Sampled",
      grepl("Par1", group) ~ "1 True Parent Sampled",
      grepl("Par2", group) ~ "2 True Parents Sampled",
      group == "WildF1" ~ "Wild F1",
      group == "TestF1-PTPS100" ~ "Test F1: 100% Parents Sampled",
      TRUE ~ "Other"
    )
  )

legend_colors <- c(
  "Test F1: 100% Parents Sampled" = "forestgreen",
  "0 True Parents Sampled" = "lightgrey",
  "1 True Parent Sampled" = "dodgerblue3",
  "2 True Parents Sampled" = "dodgerblue4",
  "Wild F1" = "salmon"
)

# library(ggplot2)

ggplot(full_to_plot, aes(x = group, y = score, fill = legend_group)) +
  geom_violin(color = "black", trim = TRUE) +
  scale_fill_manual(
    values = legend_colors,
    name = "Legend\n\nTrue Parents Sampled Per Indiv."
  ) +
  guides(fill = guide_legend(
    override.aes = list(size = 5),
    title.position = "top",
    title.hjust = 0.5
  )) +
  labs(
    title = "HOTHIPHOP Scores for Test and Wild Individuals", 
    x = "Parentage Scenario", 
    y = "HOTHIPHOP Score"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 11),
    legend.title = element_text(face = "bold", size = 11)
)







       