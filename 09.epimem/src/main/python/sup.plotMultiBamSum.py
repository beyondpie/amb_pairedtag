"""Plot multiBamSummary.
"""
import numpy as np
import os
from typing import List
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import tmpPypkg.globalvar as gv
from tmpPypkg.globalvar import reorder_scs
from tmpPypkg.globalvar import transform_allenlabel
from tmpPypkg.bam import call_plotcor
import pickle

# from importlib import reload
# reload(matplotlib)

projd = gv.pt_projd
workd = os.path.join(projd, "06.ChromHMM")
multiBamSumd = os.path.join(workd, "out/multiBamSum")
plotd = multiBamSumd
binSize:int = 10000
exclude_scs: List[str] = [
  "042 OB-out Frmd7 Gaba",
  "044 OB Dopa-Gaba",
  "070 LSX Prdm12 Slit2 Gaba",
  "105 TMd-DMH Foxd2 Gaba",
  "138 PH Pitx2 Glut",
  "140 PMd-LHA Foxb1 Glut",
  "143 MM-ant Foxb1 Glut",
  "191 PAG-MRN Rln3 Gaba"
]
exclude_scs = [transform_allenlabel(i) for i in exclude_scs]

# * how to read multiBamSum
# hm = "H3K27me3"
# binSize = 10000
# load read_multiBamSum from tmpPypkg.bam
# mat, labels = read_multiBamSum(multiBamSumd, hm, binSize)

# * plot cor
# call_plotcor from tmpPypkg.bam
for hm in ["H3K27me3", "H3K27ac", "H3K4me1", "H3K9me3"]:
  # for m in ["pearson", "spearman"]:
  for m in ["spearman"]:
    print(f"{hm} using {m}.")
    call_plotcor(
      npzfnm = os.path.join(multiBamSumd, hm, f"{hm}_binSize_{binSize}.npz"),
      plotd = plotd,
      hm = hm,
      cormthd = m
    )

# * Now plot by fixing the order
def cleanstr(t: str, hm: str)-> str:
  return t.replace("'", "").replace(".srt", "").replace(f".{hm}", "")

def tofloat(line: str)-> List[float]:
  return [float(i) for i in line.split("\t")[1:]]

def readCorMat(
  hm: str,
  multiBamSumd: str,
  cormthd: str = "pearson",
  exclude_scs: List[str] = []) -> pd.DataFrame:
  corfnm = os.path.join(
      multiBamSumd, hm,
      f"Heatmap_{hm}.{cormthd}.binSize{binSize}.corMatrix.tab")
  with open(corfnm, 'r') as f:
    lines = [line for line in f.readlines() if "#" not in line]
  heads = [cleanstr(t, hm) for t in lines[0].strip().split("\t") if len(t) > 0]
  colnms = [cleanstr(line.split("\t")[0], hm) for line in lines[1:]]
  r = [ tofloat(line) for line in lines[1:]]
  corm: pd.DataFrame = pd.DataFrame(
    data = r,
    index = colnms,
    columns = heads, dtype = np.float32)
  # reorder of dataframe
  scs_ord = reorder_scs(heads, sep = "_")
  if len(exclude_scs) > 0:
    scs_ord = [i for i in scs_ord if i not in set(exclude_scs)]
  corm_ord = corm.loc[scs_ord, scs_ord]
  return corm_ord

cormthd = "pearson"
corm_dict = {hm: readCorMat(hm, multiBamSumd, cormthd, exclude_scs)
             for hm in ["H3K27me3", "H3K27ac", "H3K4me1", "H3K9me3"]}

cormthd = "spearman"
corm_spearman_dict = {hm: readCorMat(hm, multiBamSumd, cormthd, exclude_scs)
             for hm in ["H3K27me3", "H3K27ac", "H3K4me1", "H3K9me3"]}

# * draw tri by seaborn
mask = np.triu(np.ones_like(corm, dtype=bool))
# f, axs = plt.subplots(figsize=(30 * 4, 20), ncols = 4, nrows=1)
# # diverging_colors = sns.color_palette("RdBu", 30)
cmap = sns.diverging_palette(230, 20, as_cmap=True)

def getHeatmap(
  title:str, corm: pd.DataFrame, outfnm,
  width:int = 20, height:int = 10,
  vmax: float = None, vmin: float = None,
  center: float = None):
  f, ax = plt.subplots(figsize=(width, height))
  sns.heatmap(data=corm,
              mask=mask,
              cmap=cmap,
              vmax=vmax,
              vmin=vmin,
              center=center,
              square=True,
              linewidths=.5,
              cbar_kws={"shrink": .5},
              cbar = True,
              xticklabels=True,
              yticklabels=True,
              ax=ax)
  ax.set_title(label = title)
  f.savefig(fname = outfnm, format = "pdf")
  return f

for hm in ["H3K27me3", "H3K27ac", "H3K4me1", "H3K9me3"]:
  getHeatmap(
    title = f"{hm} pearson of bams with bins of {binSize} size.",
    corm = corm_dict[hm],
    outfnm = os.path.join(plotd, f"Heatmap.{hm}.{cormthd}.{binSize}bin.pdf"),
    width = 30,
    height = 25
  )



for hm in ["H3K27me3", "H3K27ac", "H3K4me1", "H3K9me3"]:
  getHeatmap(
    title = f"{hm} spearman of bams with bins of {binSize} size.",
    corm = corm_spearman_dict[hm],
    outfnm = os.path.join(plotd, f"Heatmap.{hm}.{cormthd}.{binSize}bin.pdf"),
    width = 30,
    height = 25
  )






