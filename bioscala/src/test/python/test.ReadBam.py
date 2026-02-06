import os
import sys
import pysam

projd = "/projects/ps-renlab2/szu/projects/amb_pairedtag"
rscd = os.path.join(projd, "src/test/resource")
bam1 = os.path.join(rscd, "H3K4me1.srt.bam")

bamF = pysam.Samfile(bam1, 'rb')

r = bamF.fetch(until_eof=True)

read = next(r)
read.reference_name
read.reference_start
read.is_reverse


import pandas as pd
import anndata as ad
cellmeta = pd.read_csv("/projects/ps-renlab2/szu/projects/amb_pairedtag/meta/pairedtag.cell.meta.all.240626.csv",
                       sep = ",", header = 0)


ptRNA = ad.read_h5ad("/projects/ps-renlab2/szu/projects/amb_pairedtag/data/pairedtag_ann/pt.rawRNAseq.sa2ann.genesymbol.h5ad", backed = 'r')
