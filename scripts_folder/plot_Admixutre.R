args <- commandArgs(TRUE)
infile <- args[1]
info <- args[2]

pdf(file=paste(infile,".pdf", sep =""))
pdf(file="isophya_admxiture_altd.pdf")
pop<-read.table(info,as.is=T)
admix<-t(as.matrix(read.table(infile)))
admix<-admix[,order(pop[,3])]
pop<-pop[order(pop[,3]),]
isoPalette=c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black", "gold1", "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon", "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise", "green1", "yellow4", "yellow3", "darkorange4", "brown") ### 25 Colors
h<-barplot(admix,col=isoPalette,  las=2, space=0,border=NA,xlab="Individuals",ylab="admixture")
text(tapply(1:nrow(pop),pop[,3],mean), -0.07, unique(pop[,1]),xpd=T, srt = -90)
abline(v = c(5,12,19,26,33,40,47,54,61,67,74,81), lty = 2, col = "black")

unlink("Rplots.pdf", force=TRUE)
dev.off();
