library(tidyverse)

# * meta
samtools <- "/home/szu/miniforge3/bin/samtools"
projd <- here::here()
cellMeta <- file.path(
  projd, "meta",
  "pairedtag.cell.meta.all.240626.csv"
) |>
  data.table::fread(file = _, sep = ",", header = TRUE, data.table = FALSE)
rownames(cellMeta) <- cellMeta$barcode
cellMeta <- cellMeta[cellMeta$annotQuality == "Good", ]

ptDNAbamd <- file.path(projd, "data", "ptDNAbam")
ptDNAbamMetaDir <- file.path(ptDNAbamd, "meta")
ptDNAbamsDir <- file.path(ptDNAbamd, "bam")
ptDNAbamTagDir <- file.path(ptDNAbamd, "tag")

all_sc <- unique(cellMeta$annot.sc)
neu_sc <- all_sc[vapply(all_sc, \(i) {
  str_split_1(i, " ")[3]
}, "NN") != "NN"]

all_sp <- unique(cellMeta$annot.sp)
nn_sp <- all_sp[grepl("NN", vapply(all_sp, \(i) {
  str_split_1(i, " ")[3]
}, "NN_1"))]

modality <- c("H3K4me1", "H3K9me3", "H3K27ac", "H3K27me3")
part <- c("MaleA", "MaleB", "FemaleA", "FemaleB")


# * check merge
# 1. use enough barcodes
getAllenAnnot2Str <- function(allenAnnot) {
  gsub("/", "_", allenAnnot) |>
    gsub(" +", "_", x = _) |>
    gsub("-", "_", x = _)
}

getBarcodes4NeuGroup <- function(annotsc, modality, part) {
  index <- (cellMeta$annot.sc == annotsc) &
    (cellMeta$exp.modularity == modality) &
    (cellMeta$sample.rep == part)
  return(cellMeta$barcode[index])
}


getBarcodes4NNGroup <- function(annotsp, modality, part) {
  index <- (cellMeta$annot.sp == annotsp) &
    (cellMeta$exp.modularity == modality) &
    (cellMeta$sample.rep == part)
  return(cellMeta$barcode[index])
}

getBarcodesFromMetaFile <- function(annotsc, modality, part) {
  fnm <- file.path(ptDNAbamMetaDir, paste0(paste(
    getAllenAnnot2Str(annotsc), modality, part,
    sep = "-"), ".barcode.bams.txt")
  )
  data.table::fread(
    file = fnm,
    header = FALSE, data.table = FALSE, sep = ","
  )[, 1] |> vapply(X = _, FUN = \(l) {
      str_split_1(basename(l), "\\.")[1]
    }, FUN.VALUE = "I01:BA")
}

checkNeuBarcodeMeta <- function(annotsc, modality, part) {
  b1 <- getBarcodes4NeuGroup(annotsc, modality, part)
  b2 <- getBarcodesFromMetaFile(annotsc, modality, part)
  if ((length(b1) == length(b2)) && (all(b1 %in% b2))) {
    TRUE
  } else {
    FALSE
  }
}

checkNNBarcodeMeta <- function(annotsp, modality, part) {
  b1 <- getBarcodes4NNGroup(annotsp, modality, part)
  b2 <- getBarcodesFromMetaFile(annotsp, modality, part)
  if ((length(b1) == length(b2)) && (all(b1 %in% b2))) {
    TRUE
  } else {
    FALSE
  }
}

# pass
for (sp in nn_sp) {
  for (m in modality) {
    for (p in part) {
      torf <- checkNNBarcodeMeta(sp, m, p)
      message(paste(sp, m, p, torf))
    }
  }
}

# some of them are missing (recored in test/resource)
for (sc in neu_sc) {
  for (m in modality) {
    for (p in part) {
      tryCatch(
        {
          torf <- checkNeuBarcodeMeta(sc, m, p)
        },
        error = function(e) {
          message(paste(sc, m, p, sep = "-"))
        }
      )
    }
  }
}


# 2. all the merging jobs are done
# check tag file and exists of bam file
checkTag <- function(sc, m, p) {
  f1 <- file.path(
    ptDNAbamTagDir, 
    paste0(paste(getAllenAnnot2Str(sc), m, p, sep = "-"), ".merge.done")
  )
  f2 <- file.path(
    ptDNAbamsDir, getAllenAnnot2Str(sc),
    paste0(paste(m, p, sep = "-"), ".srt.bam")
  )
  if (file.exists(f1) && file.exists(f2)) {
    TRUE
  } else {
    FALSE
  }
}

# pass
for (sp in nn_sp) {
  for (m in modality) {
    for (p in part) {
      torf <- checkTag(sp, m, p)
      message(paste(sp, m, p, torf))
    }
  }
}
# just the ones not detected as above.
for (sc in neu_sc) {
  for (m in modality) {
    for (p in part) {
      if (!checkTag(sc, m, p)) {
        message(paste(sc, m, p))
      }
    }
  }
}

## check neu sc
all_neu_scgroups <- cellMeta[cellMeta$annot.sc %in% neu_sc, ] |>
  with(data = _,
    paste(getAllenAnnot2Str(annot.sc), exp.modularity, sample.rep, sep = "-")
  ) |> unique()

# just the same 40 groups are missing.
for (sc in neu_sc) {
  for (m in modality) {
    for (p in part) {
      if (!paste(
        getAllenAnnot2Str(sc), m, p,
        sep = "-"
      ) %in%
        all_neu_scgroups) {
        message(paste(sc, m, p))
      }
    }
  }
}

# * get stats of # of cells and # of reads per sample
# ** neu sc
# Ref: https://www.biostars.org/p/138116/
# test:
# bamfnm <- file.path(ptDNAbamsDir,
#   "001_CLA_EPd_CTX_Car3_Glut", "H3K4me1-Female.srt.bam")
getReadsPerBam <- function(bamfnm) {
  cmd <- paste(samtools, "view -F 0x904 -c", bamfnm)
  as.integer(system(command = cmd, intern = TRUE))
}

t <- cellMeta[cellMeta$annot.sc %in% neu_sc, ] |>
  group_by(annot.sc, exp.modularity, sample.rep) |>
  summarise(ncell = length(barcode))

nreads <- vapply(seq_len(nrow(t)), \(i) {
  bamfnm <- file.path(
    ptDNAbamsDir, getAllenAnnot2Str(t[i, "annot.sc"]),
      paste0(paste(t[i, "exp.modularity"], t[i, "sample.rep"], sep = "-"),
        ".srt.bam"))
  message(bamfnm)
  getReadsPerBam(bamfnm)
  }, 1)

t$nReads <- nreads
t <- t[order(t$nReads), ]

data.table::fwrite(
  x = t,
  file = file.path(projd, "04.peakcalling",
    "src", "main","resource", "Neu.subclass.SexRep.stat.csv"),
  sep = ",",
  col.names = TRUE,
  row.names = FALSE
)

# ** nn sp
tnn <- cellMeta[cellMeta$annot.sp %in% nn_sp, ] |>
  group_by(annot.sp, exp.modularity, sample.rep) |>
  summarise(ncell = length(barcode))

nreadsNN <- vapply(seq_len(nrow(tnn)), \(i) {
  bamfnm <- file.path(
    ptDNAbamsDir, getAllenAnnot2Str(tnn[i, "annot.sp"]),
    paste0(paste(tnn[i, "exp.modularity"],
      tnn[i, "sample.rep"], sep = "-"),
        ".srt.bam"))
  message(bamfnm)
  getReadsPerBam(bamfnm)
}, 1)

tnn$nReads <- nreadsNN
tnn <- tnn[order(tnn$nReads), ]
data.table::fwrite(
  x = tnn,
  file = file.path(projd, "04.peakcalling",
    "src", "main","resource", "NN.supertype.SexRep.stat.csv"),
  sep = ",",
  col.names = TRUE,
  row.names = FALSE
)

