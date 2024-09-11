#!/usr/bin/env Rscript

#You need to change the environmental variables for your own use. 

# Load required libraries
if (!require("corrplot")) {
  install.packages("corrplot", repos = "http://cran.us.r-project.org")
}

library(corrplot)

# Load your environmental data into a data frame
dataset <- data.frame(
  Elevation = elev$Elevation, 
  EVI = evi$EVI, 
  MTAS = mtas$MTAS, 
  NDVI = ndvi$NDVI,
  Precipitation = prec5$prec5, 
  SolarRad = solar5month$Solarad5, 
  TSA = ts5$TSA
)

# Check the dataset structure
print(dataset)

# Calculate the correlation matrix using Pearson method
cor_matrix <- cor(dataset, method = "pearson")

# Print the correlation matrix to the console
print(cor_matrix)

# Define variable names for better labeling in the plot
var_names <- c("Elevation", "EVI.AS", "MT.AS", "NDVI.AS", "P.AS", "SD.AS", "TS.AS")
colnames(cor_matrix) <- rownames(cor_matrix) <- var_names

# Define a color palette for the correlation plot
my_color_palette <- colorRampPalette(c("darkred", "white", "darkblue"))(n = 100)

# Generate the correlation plot
corrplot(
  cor_matrix,
  method = "color",         # Use colored squares to represent correlations
  type = "upper",           # Show only the upper triangle of the matrix
  order = "hclust",         # Cluster variables based on correlation
  tl.col = "black",         # Set text label color
  col = my_color_palette,   # Apply the custom color palette
  addCoef.col = "black"     # Add correlation coefficients to the plot
)
