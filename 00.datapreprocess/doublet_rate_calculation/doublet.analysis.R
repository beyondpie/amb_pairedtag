# 2024-02-12 Paired-Tag Doublet Analysis
# Zhaoning Wang

setwd("~/Desktop/Paired-Tag_QC_reports/amb_doublet/")
###################################
mouse.barcodes <- read.table(file = "./A_barcodes_mm10.tsv")
mouse.barcodes.split <-matrix(ncol=4, byrow = T, unlist(strsplit(as.character(mouse.barcodes$V1),split=":")))
table(mouse.barcodes.split[,4])

human.barcodes <- read.table(file = "./A_Exp1.cat_hg38.xls")
human.barcodes.split <-matrix(ncol=4, byrow = T, unlist(strsplit(as.character(human.barcodes$V1),split=":")))
table(human.barcodes.split[,4])

overlapping.barcodes <- as.data.frame(intersect(mouse.barcodes$V1, human.barcodes$V1))
colnames(overlapping.barcodes) <- "V1"
overlapping.barcodes.split <-matrix(ncol=4, byrow = T, unlist(strsplit(as.character(overlapping.barcodes$V1),split=":")))

#################################
mouse.barcodes <- read.table(file = "./B_barcodes_mm10.tsv")
mouse.barcodes.split <-matrix(ncol=4, byrow = T, unlist(strsplit(as.character(mouse.barcodes$V1),split=":")))
table(mouse.barcodes.split[,4])

human.barcodes <- read.table(file = "./B_Exp2.cat_hg38.xls")
human.barcodes.split <-matrix(ncol=4, byrow = T, unlist(strsplit(as.character(human.barcodes$V1),split=":")))
table(human.barcodes.split[,4])

overlapping.barcodes <- as.data.frame(intersect(mouse.barcodes$V1, human.barcodes$V1))
colnames(overlapping.barcodes) <- "V1"
overlapping.barcodes.split <-matrix(ncol=4, byrow = T, unlist(strsplit(as.character(overlapping.barcodes$V1),split=":")))
table(overlapping.barcodes.split[,4])


#################################
mouse.barcodes <- read.table(file = "./C_barcodes_mm10.tsv")
mouse.barcodes.split <-matrix(ncol=4, byrow = T, unlist(strsplit(as.character(mouse.barcodes$V1),split=":")))
table(mouse.barcodes.split[,4])

human.barcodes <- read.table(file = "./C_Exp3.cat_hg38.xls")
human.barcodes.split <-matrix(ncol=4, byrow = T, unlist(strsplit(as.character(human.barcodes$V1),split=":")))
table(human.barcodes.split[,4])

overlapping.barcodes <- as.data.frame(intersect(mouse.barcodes$V1, human.barcodes$V1))
colnames(overlapping.barcodes) <- "V1"
overlapping.barcodes.split <-matrix(ncol=4, byrow = T, unlist(strsplit(as.character(overlapping.barcodes$V1),split=":")))
table(overlapping.barcodes.split[,4])



#################################
mouse.barcodes <- read.table(file = "./D_barcodes_mm10.tsv")
mouse.barcodes.split <-matrix(ncol=4, byrow = T, unlist(strsplit(as.character(mouse.barcodes$V1),split=":")))
table(mouse.barcodes.split[,4])

human.barcodes <- read.table(file = "./D_Exp4.cat_hg38.xls")
human.barcodes.split <-matrix(ncol=4, byrow = T, unlist(strsplit(as.character(human.barcodes$V1),split=":")))
table(human.barcodes.split[,4])

overlapping.barcodes <- as.data.frame(intersect(mouse.barcodes$V1, human.barcodes$V1))
colnames(overlapping.barcodes) <- "V1"
overlapping.barcodes.split <-matrix(ncol=4, byrow = T, unlist(strsplit(as.character(overlapping.barcodes$V1),split=":")))
table(overlapping.barcodes.split[,4])



#################################
mouse.barcodes <- read.table(file = "./E_barcodes_mm10.tsv")
mouse.barcodes.split <-matrix(ncol=4, byrow = T, unlist(strsplit(as.character(mouse.barcodes$V1),split=":")))
table(mouse.barcodes.split[,4])

human.barcodes <- read.table(file = "./E_Exp5.cat_hg38.xls")
human.barcodes.split <-matrix(ncol=4, byrow = T, unlist(strsplit(as.character(human.barcodes$V1),split=":")))
table(human.barcodes.split[,4])

overlapping.barcodes <- as.data.frame(intersect(mouse.barcodes$V1, human.barcodes$V1))
colnames(overlapping.barcodes) <- "V1"
overlapping.barcodes.split <-matrix(ncol=4, byrow = T, unlist(strsplit(as.character(overlapping.barcodes$V1),split=":")))
table(overlapping.barcodes.split[,4])

#################################
mouse.barcodes <- read.table(file = "./F_barcodes_mm10.tsv")
mouse.barcodes.split <-matrix(ncol=4, byrow = T, unlist(strsplit(as.character(mouse.barcodes$V1),split=":")))
table(mouse.barcodes.split[,4])

human.barcodes <- read.table(file = "./F_Exp6.cat_hg38.xls")
human.barcodes.split <-matrix(ncol=4, byrow = T, unlist(strsplit(as.character(human.barcodes$V1),split=":")))
table(human.barcodes.split[,4])

overlapping.barcodes <- as.data.frame(intersect(mouse.barcodes$V1, human.barcodes$V1))
colnames(overlapping.barcodes) <- "V1"
overlapping.barcodes.split <-matrix(ncol=4, byrow = T, unlist(strsplit(as.character(overlapping.barcodes$V1),split=":")))
table(overlapping.barcodes.split[,4])



#################################
mouse.barcodes <- read.table(file = "./F_barcodes_mm10.tsv")
mouse.barcodes.split <-matrix(ncol=4, byrow = T, unlist(strsplit(as.character(mouse.barcodes$V1),split=":")))
table(mouse.barcodes.split[,4])

human.barcodes <- read.table(file = "./F_Exp6.cat_hg38.xls")
human.barcodes.split <-matrix(ncol=4, byrow = T, unlist(strsplit(as.character(human.barcodes$V1),split=":")))
table(human.barcodes.split[,4])

overlapping.barcodes <- as.data.frame(intersect(mouse.barcodes$V1, human.barcodes$V1))
colnames(overlapping.barcodes) <- "V1"
overlapping.barcodes.split <-matrix(ncol=4, byrow = T, unlist(strsplit(as.character(overlapping.barcodes$V1),split=":")))
table(overlapping.barcodes.split[,4])



#################################
mouse.barcodes <- read.table(file = "./G_barcodes_mm10.tsv")
mouse.barcodes.split <-matrix(ncol=4, byrow = T, unlist(strsplit(as.character(mouse.barcodes$V1),split=":")))
table(mouse.barcodes.split[,4])

human.barcodes <- read.table(file = "./G_Exp7.cat_hg38.xls")
human.barcodes.split <-matrix(ncol=4, byrow = T, unlist(strsplit(as.character(human.barcodes$V1),split=":")))
table(human.barcodes.split[,4])

overlapping.barcodes <- as.data.frame(intersect(mouse.barcodes$V1, human.barcodes$V1))
colnames(overlapping.barcodes) <- "V1"
overlapping.barcodes.split <-matrix(ncol=4, byrow = T, unlist(strsplit(as.character(overlapping.barcodes$V1),split=":")))
table(overlapping.barcodes.split[,4])


#################################
mouse.barcodes <- read.table(file = "./H_barcodes_mm10.tsv")
mouse.barcodes.split <-matrix(ncol=4, byrow = T, unlist(strsplit(as.character(mouse.barcodes$V1),split=":")))
table(mouse.barcodes.split[,4])

human.barcodes <- read.table(file = "./H_Exp8.cat_hg38.xls")
human.barcodes.split <-matrix(ncol=4, byrow = T, unlist(strsplit(as.character(human.barcodes$V1),split=":")))
table(human.barcodes.split[,4])

overlapping.barcodes <- as.data.frame(intersect(mouse.barcodes$V1, human.barcodes$V1))
colnames(overlapping.barcodes) <- "V1"
overlapping.barcodes.split <-matrix(ncol=4, byrow = T, unlist(strsplit(as.character(overlapping.barcodes$V1),split=":")))
table(overlapping.barcodes.split[,4])



#################################
mouse.barcodes <- read.table(file = "./I_barcodes_mm10.tsv")
mouse.barcodes.split <-matrix(ncol=4, byrow = T, unlist(strsplit(as.character(mouse.barcodes$V1),split=":")))
table(mouse.barcodes.split[,4])

human.barcodes <- read.table(file = "./I_Exp9.cat_hg38.xls")
human.barcodes.split <-matrix(ncol=4, byrow = T, unlist(strsplit(as.character(human.barcodes$V1),split=":")))
table(human.barcodes.split[,4])

overlapping.barcodes <- as.data.frame(intersect(mouse.barcodes$V1, human.barcodes$V1))
colnames(overlapping.barcodes) <- "V1"
overlapping.barcodes.split <-matrix(ncol=4, byrow = T, unlist(strsplit(as.character(overlapping.barcodes$V1),split=":")))
table(overlapping.barcodes.split[,4])




#################################
mouse.barcodes <- read.table(file = "./JO_barcodes_mm10.tsv")
mouse.barcodes.split <-matrix(ncol=4, byrow = T, unlist(strsplit(as.character(mouse.barcodes$V1),split=":")))
table(mouse.barcodes.split[,4])

human.barcodes <- read.table(file = "./JO_Exp10.cat_hg38.xls")
human.barcodes.split <-matrix(ncol=4, byrow = T, unlist(strsplit(as.character(human.barcodes$V1),split=":")))
table(human.barcodes.split[,4])

overlapping.barcodes <- as.data.frame(intersect(mouse.barcodes$V1, human.barcodes$V1))
colnames(overlapping.barcodes) <- "V1"
overlapping.barcodes.split <-matrix(ncol=4, byrow = T, unlist(strsplit(as.character(overlapping.barcodes$V1),split=":")))
table(overlapping.barcodes.split[,4])


#################################
mouse.barcodes <- read.table(file = "./K_barcodes_mm10.tsv")
mouse.barcodes.split <-matrix(ncol=4, byrow = T, unlist(strsplit(as.character(mouse.barcodes$V1),split=":")))
table(mouse.barcodes.split[,4])

human.barcodes <- read.table(file = "./K_Exp11.cat_hg38.xls")
human.barcodes.split <-matrix(ncol=4, byrow = T, unlist(strsplit(as.character(human.barcodes$V1),split=":")))
table(human.barcodes.split[,4])

overlapping.barcodes <- as.data.frame(intersect(mouse.barcodes$V1, human.barcodes$V1))
colnames(overlapping.barcodes) <- "V1"
overlapping.barcodes.split <-matrix(ncol=4, byrow = T, unlist(strsplit(as.character(overlapping.barcodes$V1),split=":")))
table(overlapping.barcodes.split[,4])


#################################
mouse.barcodes <- read.table(file = "./L_barcodes_mm10.tsv")
mouse.barcodes.split <-matrix(ncol=4, byrow = T, unlist(strsplit(as.character(mouse.barcodes$V1),split=":")))
table(mouse.barcodes.split[,4])

human.barcodes <- read.table(file = "./L_Exp12.cat_hg38.xls")
human.barcodes.split <-matrix(ncol=4, byrow = T, unlist(strsplit(as.character(human.barcodes$V1),split=":")))
table(human.barcodes.split[,4])

overlapping.barcodes <- as.data.frame(intersect(mouse.barcodes$V1, human.barcodes$V1))
colnames(overlapping.barcodes) <- "V1"
overlapping.barcodes.split <-matrix(ncol=4, byrow = T, unlist(strsplit(as.character(overlapping.barcodes$V1),split=":")))
table(overlapping.barcodes.split[,4])


#################################
mouse.barcodes <- read.table(file = "./M_barcodes_mm10.tsv")
mouse.barcodes.split <-matrix(ncol=4, byrow = T, unlist(strsplit(as.character(mouse.barcodes$V1),split=":")))
table(mouse.barcodes.split[,4])

human.barcodes <- read.table(file = "./M_Exp13.cat_hg38.xls")
human.barcodes.split <-matrix(ncol=4, byrow = T, unlist(strsplit(as.character(human.barcodes$V1),split=":")))
table(human.barcodes.split[,4])

overlapping.barcodes <- as.data.frame(intersect(mouse.barcodes$V1, human.barcodes$V1))
colnames(overlapping.barcodes) <- "V1"
overlapping.barcodes.split <-matrix(ncol=4, byrow = T, unlist(strsplit(as.character(overlapping.barcodes$V1),split=":")))
table(overlapping.barcodes.split[,4])

