#/usr/bin/Rscript


print(getwd())

#args= commandArgs(TRUE)


#infile = args[1]


#upload libraries

library(vcfR)

library(adegenet)

library(ade4)

gll=vcfR2genlight(read.vcfR("isophya71_final.vcf"))

#look at genlight object

x=as.matrix(gll)

gi=as.genind(x)

data <- data.frame(
     population = rep(c(1, 2), times = c(27,44))
   )
r <- factor(data$population, level = c(1,2))
print(r)

dp=dapc(gi, var.contrib = T,
        pop=r, perc.pca=80,
        scale = F,n.da = 1)

col= c("black","green2")

c <- snpzip(x, dp, plot = FALSE, loading.plot = TRUE, method = "average")

# Save the loading plot as an image file

pdf("loading_plotAVERAGE.pdf")
print(c$FS)
dev.off()

output <- capture.output(c$FS)
file_path <- "loadings_average.txt"
writeLines(output, file_path)

#For density plot
B <- dp$ind.coord[, 1]
B <- as.matrix(B)
Data <- cbind(B, data)
print(Data)
write.csv(Data, file = "Discriminate.csv", row.names = FALSE)
#Density Plot
pdf("density.pdf")
scatter(dp,1,1 ,col=col, bg="white",
        scree.da=F, legend= T,solid = 0.4)
dev.off()

#Generate Scatter plot
pdf("scatter.pdf")  # Set the file format and name (PDF in this case)
scatter(dp, legend = T, solid = .4, col = col, cex = 2, pch=c(18,19))
dev.off()

#Generate additional Scatter Plot
pdf("scatter1.pdf")
scatter(dp, scree.da=FALSE, bg="white",segcol = 5, pch=20, cell=0, cstar=0,
        col=col, posi.leg = "topright",solid=1 ,cex=3,clab=0, leg=TRUE, txt.leg=paste(c("Group 1", "Group 2")))
dev.off()

#Compoplot
pdf("compoplot.pdf")
compoplot(dp,
          lab="",
          ncol(r), xlab="individuals", col=col)
dev.off()

#Loading Plot
#pdf("loadings.pdf")
#loadingplot(dp$var.contr,main = "DAPC Loading Plot",adj = 1, srt = 270,cex.lab = 0.5, axis = 1, threshold = 0.000015)
#dev.off()

#loadings <- which(dp$var.contr[,1] > 0.000015)
write.table(loadings, "./loadings.txt")
