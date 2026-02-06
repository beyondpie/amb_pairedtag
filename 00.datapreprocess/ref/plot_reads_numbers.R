setwd("/projects/ps-renlab2/zhw063/14.MouseBrainPTexp8_20230418/")

library(Matrix)

# Previous plotting function

# plot_BC<-function(dna, rna, c1=1000, c2=1000, title=title){
#   layout(matrix(c(1,1,1,1),2,2,byrow=TRUE))
#   m<-merge(dna, rna, by="BC")
#   plot(m[,2:3],pch=19, cex=0.1, log="xy", xlab="DNA reads", ylab="RNA reads", col="grey", xlim=c(1,100000), ylim=c(1,100000), main=title)
#   lines(c(c1,c1),c(1,5000000), lty=2);lines(c(1,5000000),c(c2,c2), lty=2)
#   pf<-m[m[,2] >= c1 & m[,3] >= c2,]
#   points(pf[,2:3], pch=19, cex=0.1, col="firebrick")
#   legend("bottomright", legend=c( paste("DNA cutoff:", c1, sep=" "), paste("RNA cutoff:", c2, sep = " "), paste("# PF:", dim(pf)[1], sep=" ")), bty="n")
#   data<-na.omit(pf); val.cell<-data[,1]
#   write.table(val.cell, file=paste(title,"filt.xls", sep="."), col.names=F, row.names=F, sep="\t", quote=F)
# }


# My new plotting function! 20221012

plot_BC<-function(dna, rna, c1=1000, c2 = 10000, c3=1000, c4=10000, title=title){
  layout(matrix(c(1,1,1,1),2,2,byrow=TRUE))
  m<-merge(dna, rna, by="BC")
  plot(m[,2:3],pch=19, cex=0.1, log="xy", xlab="DNA reads", ylab="RNA reads", col="grey", xlim=c(1,100000), ylim=c(1,100000), main=title)
  lines(c(c1,c1),c(1,5000000), lty=2);lines(c(c2,c2), c(1,5000000), lty=2);lines(c(1,5000000),c(c3,c3), lty=2);lines(c(1,5000000),c(c4,c4), lty=2)
  pf<-m[m[,2] >= c1 & m[,2] <= c2 & m[,3] >= c3 & m[,3] <= c4,]
  points(pf[,2:3], pch=19, cex=0.1, col="firebrick")
  legend("bottomright", legend=c( paste("DNA cutoff:", c1,"~",c2, sep=" "), paste("RNA cutoff:", c3,"~",c4, sep = " "), paste("# PF:", dim(pf)[1], sep=" ")), bty="n")
  data<-na.omit(pf); val.cell<-data[,1]
  write.table(val.cell, file=paste(title,"filt.xls", sep="."), col.names=F, row.names=F, sep="\t", quote=F)
}






dna.path="./04.matrices/ZW711_mm10_sorted_rmdup_mtx2/" ## DNA raw matrix file
rna.path="./04.matrices/ZW713_mm10_sorted_rmdup_mtx2/" ## RNA raw matrix file

dna.mtx<-readMM(paste(dna.path,"/matrix.mtx", sep=""))
dna.bc<-read.table(paste(dna.path,"/barcodes.tsv", sep=""), head=F)
dna<-cbind(dna.bc[,1,drop = F], colSums(dna.mtx))
rna.mtx<-readMM(paste(rna.path, "/matrix.mtx", sep=""))
rna.bc<-read.table(paste(rna.path, "/barcodes.tsv", sep=""), head=F)
rna<-cbind(rna.bc[,1,drop = F], colSums(rna.mtx))
colnames(dna)<-c("BC", "n_reads_DNA")
colnames(rna)<-c("BC", "n_reads_RNA")



## This function will plot the # of reads (DNA vs RNA) plot and automaticcally generate a file with pass-filter barcodes. Please run the function multiple times with suitable cutoffs.

plot_BC(dna, rna,
        c1=1000, ### lower cutoff for #-reads/nuclei of DNA library
        c2=20000, ### higher cutoff for #-reads/nuclei of DNA library
        c3=300, ### lower cutoff for #-reads/nuclei of RNA library
        c4=20000, ### higher cutoff for #-reads/nuclei of RNA library
        title="Exp8lib73_ZW711_ZW713" ## prefix of output barcode file
)



dna.path="./04.matrices/ZW712_mm10_sorted_rmdup_mtx2/" ## DNA raw matrix file
rna.path="./04.matrices/ZW714_mm10_sorted_rmdup_mtx2/" ## RNA raw matrix file

dna.mtx<-readMM(paste(dna.path,"/matrix.mtx", sep=""))
dna.bc<-read.table(paste(dna.path,"/barcodes.tsv", sep=""), head=F)
dna<-cbind(dna.bc[,1,drop = F], colSums(dna.mtx))
rna.mtx<-readMM(paste(rna.path, "/matrix.mtx", sep=""))
rna.bc<-read.table(paste(rna.path, "/barcodes.tsv", sep=""), head=F)
rna<-cbind(rna.bc[,1,drop = F], colSums(rna.mtx))
colnames(dna)<-c("BC", "n_reads_DNA")
colnames(rna)<-c("BC", "n_reads_RNA")



## This function will plot the # of reads (DNA vs RNA) plot and automaticcally generate a file with pass-filter barcodes. Please run the function multiple times with suitable cutoffs.

plot_BC(dna, rna,
        c1=1000, ### lower cutoff for #-reads/nuclei of DNA library
        c2=20000, ### higher cutoff for #-reads/nuclei of DNA library
        c3=300, ### lower cutoff for #-reads/nuclei of RNA library
        c4=20000, ### higher cutoff for #-reads/nuclei of RNA library
        title="Exp8lib74_ZW712_ZW714" ## prefix of output barcode file
)



