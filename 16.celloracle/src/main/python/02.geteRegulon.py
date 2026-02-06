"""
Recover eRegulons from Cell Oracle.
For each subclass, we recover the CRE-TF-gene relationships.
Particularly, we will keep the expression levels, association
scores, and the corresponding statistical significances if possible.
"""

import os
import re
import sys

from typing import Tuple
from multiprocessing import Pool

import numpy as np
import pandas as pd

import tmpPypkg.globalvar as gv
from tmpPypkg.globalvar import load_ptscmeta
from tmpPypkg.globalvar import load_ATA_CPM_pbysc
from tmpPypkg.globalvar import load_AllenRNA_logCPM_scbyg
from tmpPypkg.globalvar import load_PairedTag_RNA_logCPM_scbyg
from tmpPypkg.celloracle import getCRE2tf, getCRE2gene, getTF2gene
from tmpPypkg.celloracle import geteRegulon
from tmpPypkg.celloracle import load_ppdcSum


# * meta
projd = gv.pt_projd
workd = os.path.join(projd, "16.celloracle")
eRegulond = os.path.join(workd, "out", "eRegulon")

ptscmeta = load_ptscmeta()
atacscs = ptscmeta.CEMBATACname
ptscs = ptscmeta.PairedTagName

grnd = gv.GRNd
pdcd = gv.pdcd
tfscand = gv.tfscand

atacCPMpbysc = load_ATA_CPM_pbysc()
#ptRNAlogCPMscbyg = load_PairedTag_RNA_logCPM_scbyg()
#allenRNAlogCPMscbyg = load_AllenRNA_logCPM_scbyg()
#ppdcSum: pd.DataFrame = load_ppdcSum()

# * functions
def geteRegulonParallel(x: Tuple[str, str]) -> None:
    atacsc, ptsc = x
    print(f"Get eRegulon for {ptsc}.")
    eRegulonfnm = os.path.join(eRegulond,
                               f"cellOracle.{ptsc}.eRegulon.csv")
    linkfnm = os.path.join(grnd, f"GRN.{atacsc}.celloracle.links")
    tf2gene = getTF2gene(
        linkfnm=linkfnm, ptsc=ptsc,
        allensc=atacsc, atacsc=atacsc,
        allenRNAlogCPMscbyg=allenRNAlogCPMscbyg,
        ptRNAlogCPMscbyg=ptRNAlogCPMscbyg)
    pdcfnm = os.path.join(pdcd, f"sa2subclass.{atacsc}.pdc.bedpe")
    CRE2gene = getCRE2gene(pdcfnm, ppdcSum=ppdcSum)

    tfinfofnm = os.path.join(
        tfscand, f"sa2subclass.{atacsc}.pdc.celloracle.tfinfo")

    CRE2tf = getCRE2tf(tfinfofnm, atacsc, atacCPMpbysc)

    eRegulon = geteRegulon(CRE2gene, CRE2tf, tf2gene)
    eRegulon.to_csv(eRegulonfnm,
                    header=True, sep=",", index=False)

def runCRE2TF(x: Tuple[str, str]) -> None:
    atacsc, ptsc = x
    print(f"Get CRE2tf for {ptsc}.")
    tfinfofnm = os.path.join(
        tfscand, f"sa2subclass.{atacsc}.pdc.celloracle.tfinfo")
    CRE2tf: pd.DataFrame = getCRE2tf(tfinfofnm, atacsc, atacCPMpbysc)
    outfnm = os.path.join(eRegulond, f"gimmeMotif.{ptsc}.CRE2tf.csv")
    CRE2tf.to_csv(outfnm, header=True, sep=",", index=False)
    return None


# * main
# with Pool(20) as p:
#     p.map(geteRegulonParallel, zip(atacscs, ptscs))

with Pool(30) as p:
    p.map(runCRE2TF, zip(atacscs, ptscs))
