#!/usr/bin/Rscript

#./line_plot.R --group1 group1.tsv --group2 group2.tsv --group3 group3.tsv --outfile maf_plot.png


# Load necessary libraries
if (!require(optparse)) install.packages("optparse", repos='http://cran.us.r-project.org')
if (!require(ggplot2)) install.packages("ggplot2", repos='http://cran.us.r-project.org')
if (!require(ggrepel)) install.packages("ggrepel", repos='http://cran.us.r-project.org')

library(optparse)
library(ggplot2)
library(ggrepel)

# Command-line argument parsing
option_list <- list(
  make_option(c("-g1", "--group1"), type="character", default="group1.tsv", 
              help="Input file for group 1 [default %default]", metavar="GROUP1"),
  make_option(c("-g2", "--group2"), type="character", default="group2.tsv", 
              help="Input file for group 2 [default %default]", metavar="GROUP2"),
  make_option(c("-g3", "--group3"), type="character", default="group3.tsv", 
              help="Input file for group 3 [default %default]", metavar="GROUP3"),
  make_option(c("-o", "--outfile"), type="character", default="maf_plot.png", 
              help="Output file for saving the plot [default %default]", metavar="OUTFILE")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Read the input files for each group
group1 <- read.table(opt$group1, header = TRUE, sep = '\t')
group2 <- read.table(opt$group2, header = TRUE, sep = '\t')
group3 <- read.table(opt$group3, header = TRUE, sep = '\t')

# Add group labels to each dataset
group1$group <- rep("group1", nrow(group1))
group2$group <- rep("group2", nrow(group2))
group3$group <- rep("group3", nrow(group3))

# Combine all data into one dataframe
all_data <- rbind(group1, group2, group3)

# Save plot as a PNG file
png(filename = opt$outfile, width = 800, height = 600)

# Generate the plot with ggplot2 (first version with legend outside)
ggplot(all_data, aes(x = group, y = knownEM, color = chromo_position)) +
  geom_point(size = 2) +
  geom_line(aes(group = chromo_position), size = 0.8) +
  labs(x = "Group", y = "Minor Allele Frequencies", color = "SNP Regions", 
       title = "Minor Allele Frequencies Changes for Three Groups") +
  scale_x_discrete(expand = expansion(add = 0.05)) +
  theme_minimal() +
  theme(legend.position = "right",
        plot.margin = unit(c(1, 4, 1, 1), "lines"))

# Close the graphics device (save the plot)
dev.off()

# Optionally, clean up (use with caution)
# unlink(opt$group1)
# unlink(opt$group2)
# unlink(opt$group3)
