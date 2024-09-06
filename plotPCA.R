#!/usr/bin/Rscript
# Usage: Rscript -i infile.covar -c component1-component2 -a annotation.file -o outfile.eps

.libPaths("~/R/library")

library(optparse)
library(ggplot2)

option_list <- list(make_option(c('-i','--in_file'), action='store', type='character', default=NULL, help='Input file (output from ngsCovar)'),
                    make_option(c('-c','--comp'), action='store', type='character', default=1-2, help='Components to plot'),
                    make_option(c('-a','--annot_file'), action='store', type='character', default=NULL, help='Annotation file with individual classification (2 column TSV with ID and ANNOTATION)'),
                    make_option(c('-o','--out_file'), action='store', type='character', default=NULL, help='Output file')
)
opt <- parse_args(OptionParser(option_list = option_list))

# Annotation file is in plink cluster format

#################################################################################

# Read input file
covar <- read.table(opt$in_file, stringsAsFact=F);

# Read annot file
annot <- read.table(opt$annot_file, sep=" ", header=T); # note that plink cluster files are usually tab-separated instead

# Parse components to analyze
comp <- as.numeric(strsplit(opt$comp, "-", fixed=TRUE)[[1]])



# Eigenvalues
eig <- eigen(covar, symm=TRUE);
eig$val <- eig$val/sum(eig$val);
cat(signif(eig$val, digits=3)*100,"\n");


# Write eigenvalues
#write.table(eig, file = "eigen_scores.txt", quote = FALSE)


# Plot
PC <- as.data.frame(eig$vectors)
colnames(PC) <- gsub("V", "PC", colnames(PC))
PC$Color <- factor(annot$CLUSTER)
PC$Pop <- factor(annot$IID)
PC$Lab <- factor(annot$FID)
# Write PC components
write.table(PC, file = "PC_scores.txt", quote = FALSE)


title <- paste("PC",comp[1]," (",signif(eig$val[comp[1]], digits=3)*100,"%)"," / PC",comp[2]," (",signif(eig$val[comp[2]], digits=3)*100,"%)",sep="",collapse="")

x_axis = paste("PC",comp[1],sep="")
y_axis = paste("PC",comp[2],sep="")

isoPalette=c("black","green2","dodgerblue","#6A3D9A", "#FF7F00", "dodgerblue2", "gold1", "skyblue2", "#FB9A99", "#E31A1C", "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon", "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise", "green1", "yellow4", "yellow3", "darkorange4", "brown") ### 25 Colors
### Basic plot ###
#ggplot() + geom_point(data=PC, aes_string(x=x_axis, y=y_axis, color="Pop")) + ggtitle(title) + scale_colour_manual(values=isoPalette)

### Basic plot with labeled individuals ###
##ggplot() + geom_point(data=PC, aes_string(x=x_axis, y=y_axis, color="Pop")) + ggtitle(title) + scale_colour_manual(values=isoPalette) + geom_text(data=PC, aes_string(x=x_axis, y=y_axis, label="Lab", vjust= -0.784), cex=2.5, hjust=-0.3)

### Mukti group with shapes and without labeled individuals ###

population_labels <- c("IST06","PLVYL","IST12","IST13", "IST14", "CANCK", "PIKNK", "ELEVT","VRCNK","KALEK","CKLYU")

ggplot() + geom_point(data=PC, aes_string(x=x_axis, y=y_axis, color="Color", shape= "Pop")) + ggtitle(title) + scale_colour_manual(values = isoPalette) + geom_point(size = 6)  + scale_shape_manual(values = c(1,2,3,4,5,6,7,8,9,10,11), labels = population_labels)

### Multi group plot with labeled individuals ###
#ggplot() + geom_point(data=PC, aes_string(x=x_axis, y=y_axis, color="Pop", shape="Tra")) + ggtitle(title) + scale_colour_manual(values=isoPalette) + scale_shape_manual(values = c(3,4,7,8,9,12,15,16,17)) + geom_text(data=PC, aes_string(x=x_axis, y=y_axis, label="Lab", vjust= -0.784), check_overlap=TRUE, cex=2.5, hjust=-0.3)

ggsave(opt$out_file)
unlink("Rplots.pdf", force=TRUE)
