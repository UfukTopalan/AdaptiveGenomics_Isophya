#!/usr/bin/Rscript

#Rscript plot_diversity.R --infile all_pops_diversity.tsv --outfile theta.png

# Load necessary libraries
if (!require(optparse)) install.packages("optparse", repos='http://cran.us.r-project.org')
if (!require(ggplot2)) install.packages("ggplot2", repos='http://cran.us.r-project.org')
if (!require(psych)) install.packages("psych", repos='http://cran.us.r-project.org')
if (!require(pastecs)) install.packages("pastecs", repos='http://cran.us.r-project.org')
if (!require(viridis)) install.packages("viridis", repos='http://cran.us.r-project.org')
if (!require(dplyr)) install.packages("dplyr", repos='http://cran.us.r-project.org')

library(optparse)
library(ggplot2)
library(psych)
library(pastecs)
library(viridis)
library(dplyr)

# Command-line argument parsing
option_list <- list(
  make_option(c("-i", "--infile"), type="character", default="all_pops_diversity.tsv", 
              help="Input file containing diversity data [default %default]", metavar="INFILE"),
  make_option(c("-o", "--outfile"), type="character", default="thetaP_plot.png", 
              help="Output file for saving the plot [default %default]", metavar="OUTFILE")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Read the input file
Data <- read.table(opt$infile, header = TRUE, sep = '\t')

# Normalize columns 3 and 4 by dividing by Nsites
Data[,c(3,4)] <- Data[,c(3,4)] / Data$Nsites

# Filter data based on conditions
Data <- subset(Data, Nsites >= 200)
Data <- subset(Data, tW < 0.1 & tP < 0.1)
Data <- Data %>%
  filter((.[, 3] <= 0.07) & (.[, 4] <= 0.07))

# Calculate means for each population
Pops_means <- aggregate(Data[,3:5], list(Data$Pop), FUN=mean)

# Set population order for plotting
Data$Pop <- factor(Data$Pop, levels = c("IST06", "PLVYL", "IST12", "IST13", "IST14", "CANCK", "PIKNK", "ELEVT", "VRCNK", "KALEK", "CKLYU"))

# Save plot as a PNG file
png(filename = opt$outfile, width = 800, height = 600)

# Generate the boxplot with ggplot2
ggplot(Data, aes(Pop, tP)) + 
  geom_boxplot(fill="grey", col="black", notch = TRUE) + 
  coord_flip() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 2.5, color = "#FC4E07") + 
  ggtitle(label = "ThetaP values for I. rizeensis populations") +
  scale_y_continuous(breaks=seq(0, 0.09, 0.01)) + 
  theme_bw()

# Close the graphics device (save the plot)
dev.off()

# Optionally, clean up (use with caution)
# unlink(opt$infile)
