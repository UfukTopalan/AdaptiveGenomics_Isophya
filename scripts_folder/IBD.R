#!/usr/bin/Rscript

#Rscript IBD.R --file all_pops_fst.tsv --dist Mantel_TestGeo.tsv --output ibd_plot.png

library(optparse)
library(vegan)

# Command-line argument parsing
option_list <- list(
  make_option(c("-f", "--file"), type="character", default="all_pops_fst.tsv", 
              help="Input file containing Fst values [default %default]", metavar="FILE"),
  make_option(c("-d", "--dist"), type="character", default="Mantel_TestGeo.tsv", 
              help="Geographic distance file [default %default]", metavar="DIST"),
  make_option(c("-o", "--output"), type="character", default="ibd_plot.png", 
              help="Output file for saving the plot [default %default]", metavar="OUTPUT")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Read the input files
all_pops_fst <- read.table(opt$file, sep = "", header = TRUE)
geographic_dist <- read.table(opt$dist, sep = "", header = TRUE)

# Linearize Fst values
all_pops_fst$Fst_forIBD <- all_pops_fst$Fst / (1 - all_pops_fst$Fst)

# Create distance matrix
populations <- unique(c(all_pops_fst$Pop1, all_pops_fst$Pop2))
dist_matrix <- matrix(0, nrow=length(populations), ncol=length(populations))
rownames(dist_matrix) <- populations
colnames(dist_matrix) <- populations

# Save distance matrix
write.table(dist_matrix, "dist_matrix.tsv", sep = "\t", quote = FALSE, col.names = NA)

# Create the distance matrix for linearized Fst
fstforIBD <- dist(all_pops_fst$Fst_forIBD)
geo_matrix <- as.dist(geographic_dist$Distance)

# Perform Mantel test
mantel_result <- mantel(fstforIBD, geo_matrix, method = "pearson", permutations = 999)

# Save the plot
png(filename = opt$output, width = 800, height = 600)

# Plot
plot(geo_matrix, fstforIBD, main = "Isolation By Distance in Isophya Rizeensis", 
     xlab = "Geographic Distance (m)", ylab = "Fst/(1- Fst)")
abline(lm(fstforIBD ~ geo_matrix), lwd = 3, col ="red")

# Annotate with Mantel statistics
#s_mantel_r <- mantel_result$statistic
#s_significance <- mantel_result$signif
p_mantel_r <- mantel_result$statistic
p_significance <- mantel_result$signif

mtext(side = 1, line = -4.5, at = 5, adj = -1.8, text = paste("Pearson Mantel statistic r:", p_mantel_r))
mtext(side = 1, line = -3.5, at = 5, adj = -2.2, text = paste("Pearson Significance:", p_significance))
#mtext(side = 1, line = -2.5, at = 5, adj = -1.7, text = paste("Spearman Mantel statistic r:", s_mantel_r))
#mtext(side = 1, line = -1.5, at = 5, adj = -2.1, text = paste("Spearman Significance:", s_significance))

# Close the graphics device
dev.off()

# Optionally, remove temporary files (use with caution)
# unlink("dist_matrix.tsv")
