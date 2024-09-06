#!/usr/bin/Rscript

# Load necessary libraries
if (!require(optparse)) install.packages("optparse", repos='http://cran.us.r-project.org')
if (!require(LEA)) install.packages("LEA", repos='http://cran.us.r-project.org')
if (!require(lfmm)) install.packages("lfmm", repos='http://cran.us.r-project.org')
if (!require(vcfR)) install.packages("vcfR", repos='http://cran.us.r-project.org')
if (!require(qvalue)) install.packages("qvalue", repos='http://cran.us.r-project.org')

library(optparse)
library(LEA)
library(lfmm)
library(vcfR)
library(qvalue)

# Command-line argument parsing
option_list <- list(
  make_option(c("-d", "--datafile"), type="character", default="isophya71.lfmm", 
              help="Input LFMM genotype file [default %default]", metavar="DATAFILE"),
  make_option(c("-e", "--envfile"), type="character", default="environmental_variables.txt", 
              help="Input environmental data matrix [default %default]", metavar="ENVFILE"),
  make_option(c("-o", "--outfile"), type="character", default="result_LFMM.txt", 
              help="Output file for LFMM results [default %default]", metavar="OUTFILE")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Start data analysis

# Step 1: Load the LFMM genotype data
data <- read.table(opt$datafile, sep=" ", header=FALSE, stringsAsFactors=FALSE)

# Step 2: Set the first row as column names and remove the first row
colnames(data) <- as.character(data[1, ])
data <- data[-1, ]

# Step 3: Set the first column as row names and remove the first column
rownames(data) <- as.character(data[, 1])
data <- data[, -1]

# Step 4: Make the data numeric for LFMM analysis
data[] <- lapply(data, as.numeric)

# Load the environmental data matrix (columns: environmental variables, rows: population values)
env_data <- read.table(opt$envfile, header=TRUE)

# Repeat the environmental variables for each individual in the populations
# Adjust 'repetitions' based on your specific population structure
repetitions <- c(7, 6, 7, 7, 7, 5, 5, 7, 7, 7, 6)
env_repeated <- apply(env_data, 2, function(col) rep(col, times = repetitions))
env_repeated <- as.matrix(env_repeated)

# Run the LFMM analysis
Isophya_lfmm <- lfmm_ridge(Y=data, X=env_repeated, K=2)

# Identify LFMM candidates using False Discovery Rate (FDR) of 5%
Isophya_pval <- lfmm_test(Y=data, X=env_repeated, lfmm=Isophya_lfmm, calibrate="gif")

# Extract genomic inflation factor (GIF)
gif <- Isophya_pval$gif
print(gif)

# Extract q-values and filter significant loci (q-value < 0.05)
Isophya_qv <- qvalue(Isophya_pval$calibrated.pvalue)$qvalues
selected_qvalues <- Isophya_qv[which(Isophya_qv < 0.05)]
selected_columns <- colnames(data)[which(Isophya_qv < 0.05)]

# Create a result data frame with significant loci and q-values
result_df <- data.frame(Column=selected_columns, QValue=selected_qvalues)

# Save the results to an output file
write.table(result_df, file=opt$outfile, row.names=FALSE, quote=FALSE)

print("LFMM analysis completed successfully. Results saved to file.");
