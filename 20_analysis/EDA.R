# This file contains the code for Exploratory Data Analysis of the dataset
library(ggplot2)
library(gridExtra)
library(grid)
library(corrplot)

# Load the dataset
cannabis <- readRDS(paste(here::here(), "02_cleaned_data/cannabis-2.rds", sep = "/"))

#-----------------------------------------------------------------------------#
# Number of predictors and number of samples ####
#-----------------------------------------------------------------------------#
dim(cannabis)
# 27 predictors, 57 samples


#-----------------------------------------------------------------------------#
# Distribution of the outcome variable ####
#-----------------------------------------------------------------------------#
# Density and Histogram of raw t_mmr1
plot_hist <- ggplot(cannabis, aes(x = t_mmr1)) +
  geom_histogram(aes(y=after_stat(count / sum(count))), bins=25, fill="lightgrey", color="black") +  # scale histogram y
  geom_density(col = "maroon") +
  labs(title = "Distribution of Metabolite Molar Ratio",
       x = "Metabolite Molar Ratio",
       y = "Density") +
  theme_minimal()

# Log-transformed t_mmr1
plot_hist_log <- ggplot(cannabis, aes(x = log(t_mmr1+10e-6))) +
  geom_histogram(aes(y=after_stat(count / sum(count))), bins=25, fill="lightgrey", color="black") +  # scale histogram y
  geom_density(col = "maroon") +
  labs(title = "Distribution of Log-Transformed Metabolite Molar Ratio",
       subtitle = "Note: 10e-6 was added to avoid taking the logarithm of 0",
       x = "Log Metabolite Molar Ratio",
       y = "Density") +
  theme_minimal()

# Arrange Plots
grid.arrange(plot_hist, plot_hist_log, ncol=2)

#-----------------------------------------------------------------------------#
# Variable Correlation ####
#-----------------------------------------------------------------------------#
# Create a correlation grid
corr_matrix <- cor(cannabis[,-c(1)])
corrplot(corr_matrix, method = "circle", 
         type="lower", tl.cex=0.75, tl.col="black")

# Find top 10 correlated variable pairs
corr_matrix[lower.tri(corr_matrix)] <- NA # Remove the lower triangle of corr_matrix
corr_pairs <- as.data.frame(as.table(corr_matrix))
names(corr_pairs) <- c("Var1", "Var2", "Corr")
corr_pairs <- corr_pairs[corr_pairs$Var1 != corr_pairs$Var2,] # Remove corr==1
corr_pairs <- corr_pairs[!is.na(corr_pairs$Corr),] # Remove corr=NA
corr_pairs <- corr_pairs[order(-abs(corr_pairs$Corr)),]

