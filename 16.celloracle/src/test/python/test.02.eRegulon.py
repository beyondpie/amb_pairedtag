import os

import pandas as pd

from tmpPypkg.celloracle import getCRE2tf, getCRE2gene, getTF2gene
from tmpPypkg.celloracle import geteRegulon
from tmpPypkg.celloracle import load_ppdcSum

import tmpPypkg.globalvar as gv
from tmpPypkg.globalvar import load_ATA_CPM_pbysc
from tmpPypkg.globalvar import load_AllenRNA_logCPM_scbyg
from tmpPypkg.globalvar import load_PairedTag_RNA_logCPM_scbyg

import celloracle as co
from celloracle.network_analysis.links_object import Links
from celloracle.motif_analysis.tfinfo_core import load_TFinfo, TFinfo

projd = gv.pt_projd
workd = os.path.join(projd, "16.celloracle")
outd = os.path.join(workd, "out")

grnd = gv.GRNd
pdcd = gv.pdcd
tfscand = gv.tfscand

atacCPMpbysc = load_ATA_CPM_pbysc()
ptRNAlogCPMscbyg = load_PairedTag_RNA_logCPM_scbyg()
allenRNAlogCPMscbyg = load_AllenRNA_logCPM_scbyg()
ppdcSum: pd.DataFrame = load_ppdcSum()

atacsc = "Astro-TE_NN"
allensc = "Astro-TE_NN"
ptsc = "319_Astro_TE_NN"

# load TF2target
linkfnm = os.path.join(grnd, f"GRN.{atacsc}.celloracle.links")

tf2gene = getTF2gene(
    linkfnm=linkfnm, ptsc=ptsc,
    allensc=allensc, atacsc=atacsc,
    allenRNAlogCPMscbyg=allenRNAlogCPMscbyg,
    ptRNAlogCPMscbyg=ptRNAlogCPMscbyg
)

pdcfnm = os.path.join(pdcd, f"sa2subclass.{atacsc}.pdc.bedpe")
CRE2gene = getCRE2gene(pdcfnm, ppdcSum=ppdcSum)

tfinfofnm = os.path.join(
    tfscand, f"sa2subclass.{atacsc}.pdc.celloracle.tfinfo")

CRE2tf = getCRE2tf(tfinfofnm, atacsc, atacCPMpbysc)

eRegulon = geteRegulon(CRE2gene, CRE2tf, tf2gene)
eRegulon.to_csv(os.path.join(outd, f"test.{atacsc}.eRegulon.csv"),
                header=True, sep=",", index=False)

# * double double check the results from gimme
tfinfofnm = os.path.join(
    tfscand, f"sa2subclass.{atacsc}.pdc.celloracle.tfinfo")
tfinfo: TFinfo = load_TFinfo(tfinfofnm)
