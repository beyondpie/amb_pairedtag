library(GenomicRanges)
library(rtracklayer)
library(tmpRpkg)

# * load GFF / GTF file
# https://www.gencodegenes.org/mouse/release_M25.html
gff <- file.path("/projects/ps-renlab2/szu/genome",
  "gencode.vM25.annotation.gff3")

geneAnnot <- rtracklayer::import(gff)

# * promoter
promoters <- aroundGeneTSS(
  geneAnnot, up = 1500, down = 500)

chrs.tmp <- GenomicRanges::seqnames(promoters)
chrs.tmp <- as.factor(chrs.tmp)
chrs <- levels(chrs.tmp)[chrs.tmp]

strand.tmp <- GenomicRanges::strand(promoters)
strand.tmp <- as.factor(strand.tmp)
strands <- levels(strand.tmp)[strand.tmp]

df <- data.frame(
  chroms = chrs,
  starts = promoters@ranges@start,
  ends = end(promoters),
  gene_names = mcols(promoters)$gene_name,
  strands = strands,
  gene_types = mcols(promoters)$gene_type,
  gene_ids = mcols(promoters)$gene_id
)

outf <- file.path("/projects/ps-renlab2/szu/genome",
  "mouseGenePromoter.GRCm38-M25.TSS-up1.5k-down500.tsv")
data.table::fwrite(x = df, file = outf,
  row.names = F, col.names = F, sep = "\t")

# * promoter used in cicero and GRN
promoters <- aroundGeneTSS(geneAnnot, up = 1000, down = 1000)

chrs.tmp <- GenomicRanges::seqnames(promoters)
chrs.tmp <- as.factor(chrs.tmp)
chrs <- levels(chrs.tmp)[chrs.tmp]

strand.tmp <- GenomicRanges::strand(promoters)
strand.tmp <- as.factor(strand.tmp)
strands <- levels(strand.tmp)[strand.tmp]

df <- data.frame(
  chroms = chrs,
  starts = promoters@ranges@start,
  ends = end(promoters),
  gene_names = mcols(promoters)$gene_name,
  strands = strands,
  gene_types = mcols(promoters)$gene_type,
  gene_ids = mcols(promoters)$gene_id
)

outf <- file.path("/projects/ps-renlab2/szu/genome",
  "mouseGenePromoter.GRCm38-M25.promoter.1kbAwayTSS.tsv")
data.table::fwrite(x = df, file = outf,
  row.names = F, col.names = F, sep = "\t")
