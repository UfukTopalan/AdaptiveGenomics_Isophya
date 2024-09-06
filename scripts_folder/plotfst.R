#!/usr/bin/Rscript

#Rscript plot_fst.R --file all_pops_fst.tsv --title "Custom Title" --output "custom_fst_plot.png"

library(optparse)
library(ggplot2)
library(viridis)

# Command-line argument parsing
option_list <- list(
  make_option(c("-f", "--file"), type="character", default="all_pops_fst.tsv", 
              help="Input file containing Fst values [default %default]", metavar="FILE"),
  make_option(c("-t", "--title"), type="character", default="Fst between I. rizeensis Groups", 
              help="Plot title [default %default]", metavar="TITLE")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Read the input file
Fst <- read.table(opt$file, sep = "", header = TRUE)

# Define orders for x and y axes
x_order <- c( "PLVYL", "IST12", "IST13", "IST14", "CANCK", "PIKNK", "ELEVT", "VRCNK", "KALEK", "CKLYU" )
y_order <- c( "IST06", "PLVYL", "IST12", "IST13", "IST14", "CANCK", "PIKNK", "ELEVT", "VRCNK", "KALEK" )

Fst$Pop1 <- factor(Fst$Pop1, levels = y_order)
Fst$Pop2 <- factor(Fst$Pop2, levels = x_order)

# Create the plot
ggplot(Fst, aes(Pop2, Pop1)) + 
  geom_tile(aes(fill = Fst), color="white") +
  theme_minimal() +
  scale_fill_viridis(discrete=FALSE, option = "plasma", space = "Lab", direction = -1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45,vjust=1,size=10,hjust=1),
        panel.border=element_blank(), 
        panel.grid.major=element_blank(), 
        axis.title.y = element_blank(), 
        axis.title.x = element_blank()) +
  scale_y_discrete(position = "left") +
  geom_text(aes(label = round(Fst, 2)), size=3.5, colour = "white") +
  ggtitle(opt$title)
# Close the graphics device (save the plot)
dev.off()

# Remove the input file from the system if required (use with caution)
unlink(opt$file)
