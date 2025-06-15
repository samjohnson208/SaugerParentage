## PCA WITH SNPRelate ##

# 1. Load packages
library(vcfR)
library(SNPRelate)
library(tidyverse)
library(gdsfmt)
library(dplyr)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SNPRelate")
library(SNPRelate)

# 2. Upload vcf
setwd("/Users/samjohnson/Documents/Sauger_042225/GeneticData/vcfs/pflav_mem/")
vcf <- read.vcfR("hard_variants_pflav_bial_noindels_q20_mindep8_maxdep75_maf1_miss95.recode.vcf", verbose = FALSE)

# 3. From vcf to GDS file, needed for the pca
snpgdsVCF2GDS("/Users/samjohnson/Documents/Sauger_042225/GeneticData/vcfs/pflav_mem/hard_variants_pflav_bial_noindels_q20_mindep8_maxdep75_maf1_miss95.recode.vcf", 
              "/Users/samjohnson/Documents/Sauger_042225/GeneticData/vcfs/pflav_mem/hard_variants_pflav_bial_noindels_q20_mindep8_maxdep75_maf1_miss95.recode.gds") #first argument is the location of the vcf, the second is where you want to save the created GDS file

# then upload the created GDS 
genofile <- snpgdsOpen("/Users/samjohnson/Documents/Sauger_042225/GeneticData/vcfs/pflav_mem/hard_variants_pflav_bial_noindels_q20_mindep8_maxdep75_maf1_miss95.recode.gds")

# 4. PCA
# you only need to indicate the uploaded GDS. The default is that it will create 32 PC's
pca <- snpgdsPCA(genofile, autosome.only=F) ## with all SNPs, none excluded


# 5. Extra code for creating a df with the loadings for each sample and the proportion of variance explained by the PC's
df_pca <- as.data.frame(pca$eigenvect)
colnames(df_pca) <- paste("PC", 1:32, sep= "")

df_pca <- df_pca %>% 
  mutate("Individual.name" = colnames(vcf@gt)[-1]) %>% # here you are extracting the samples' names from the vcf, to keep the order (which was the same used for the PCA)
  relocate(Individual.name, .before = PC1) %>% 
  mutate("Var.Prop" = pca$varprop) %>% 
  relocate(Var.Prop, .after = PC32)

df_pca <- df_pca %>% 
    rename(Sample_ID = Individual.name)

sar <- read.csv(file = "SAR_Data_092424.csv", header = TRUE)

sar_red <- sar %>% 
    select(Sample_ID, Group)

df_pca_sar <- df_pca %>% 
    left_join(sar_red, by = "Sample_ID") 

df_pca_sar <- df_pca_sar %>% 
  select(Sample_ID, Group, everything())

table(is.na(df_pca_sar$Group))
df_pca_sar$Group[1:22] <- "WCR"
str(df_pca_sar)

table(df_pca_sar$Group)
group_vals <- sort(unique(df_pca_sar$Group))
color_map <- c("purple", "blue", "forestgreen", "darkred", "gold", "goldenrod" ,"grey", "red")
names(color_map) <- group_vals
df_pca_sar$plot_col <- color_map[as.character(df_pca_sar$Group)]

round(df_pca_sar$Var.Prop[1]*100, 1) # 1.8
round(df_pca_sar$Var.Prop[2]*100, 1) # 1.7
round(df_pca_sar$Var.Prop[3]*100, 1) # 1.7
round(df_pca_sar$Var.Prop[4]*100, 1) # 1.5
round(df_pca_sar$Var.Prop[5]*100, 1) # 1.4

df_pca_sar <- df_pca_sar %>% 
  select(Sample_ID, Group, plot_col, everything())

df_pca_sar <- df_pca_sar[-c(556, 950, 896, 1183, 1184, 1182), ]

plot(PC2 ~ PC1, data = df_pca_sar, 
     col = plot_col, pch = 19,
     xlab = "PC1 (1.8% variance)",
     ylab = "PC2 (1.7% variance)")

plot(PC3 ~ PC2, data = df_pca_sar, 
     col = plot_col, pch = 19,
     xlab = "PC2 (1.7% variance)",
     ylab = "PC3 (1.7% variance)")

plot(PC4 ~ PC3, data = df_pca_sar, 
     col = plot_col, pch = 19,
     xlab = "PC3 (1.7% variance)",
     ylab = "PC4 (1.5% variance)")

plot(PC5 ~ PC4, data = df_pca_sar, 
     col = plot_col, pch = 19,
     xlab = "PC4 (1.5% variance)",
     ylab = "PC5 (1.4% variance)")

legend("bottomright", legend = c(sort(unique(df_pca_sar$Group))), 
                      col = c("purple", "blue", "forestgreen", "darkred", "gold", "goldenrod" ,"grey", "red"),
                      pch = 19, cex = 0.9)

write.csv(df_pca_sar, "/Users/samjohnson/Desktop/df_pca_sar.csv",
          row.names = F, quote = F)


# # 6. Plot: an example. I added more metadata stuff to that df 
# 
# library(ggplot2)
# library(ggalt)
# library(ggmagnify)
# 
# ggplot(df_pca.hard.yct, aes(x=PC1,y=PC2)) + 
#   geom_point(aes(col = factor(lib2)), size = 4, alpha=0.8) +
#   labs(x = paste("PC1 (", round(df_pca.hard.yct$Var.Prop[1]*100, 1), "% variance explained)",sep=""),
#        y = paste("PC2 (", round(df_pca.hard.yct$Var.Prop[2]*100, 1), "% variance explained)",sep="")) +
#   labs(col = "Library") +
#   #xlim(-0.021, 0.005) + ## zooming on the concentrated region
#   #ylim(-0.025, 0.03)+ ## zooming on the concentrated region
#   theme_minimal()