from typing import Tuple, List
import pysam
import re
import numpy as np
import os
import subprocess

def flatten_list(xss: List[List]) -> List:
    return [x for xs in xss for x in xs]

def get_barcodes_from_bam_sync(
        bam_fnm:str,
        index_fnm:str,
        barcode_pattern: str = ':\w\w:\w\w:\w\w:'
) -> List[str]:
    with pysam.AlignmentFile(bam_fnm,
                               index_filename = index_fnm,
                               mode = "rb") as abam:
        barcodes = [re.search(
                    barcode_pattern, r.query_name).group().strip(":")
                for r in abam.fetch()]
    # ubarcodes = list(set(barcodes))
    return barcodes

def get_barcodes_from_intervals_async(
        intervals: List[Tuple[str, int, int]],
        bam_fnm: str,
        index_fnm: str,
        barcode_pattern: str) -> List[str]:
    with pysam.AlignmentFile(bam_fnm,
                               index_filename = index_fnm,
                               mode = "rb") as abam:
        barcode = [re.search(
                    barcode_pattern, r.query_name).group().strip(":")
                for i in intervals
                for r in abam.fetch(contig = i[0], start = i[1], stop = i[2])]
    return barcode


def read_multiBamSum(
    multiBamSumd: str, hm: str, binSize:
    int = 10000)-> Tuple[np.ndarray, np.ndarray]:
  with np.load(
    os.path.join(
        multiBamSumd, hm, f"{hm}_binSize_{binSize}.npz")) as f:
    mat = f['matrix']
    labels = f['labels']
  return (mat, labels)


def call_plotcor(
  npzfnm:str,
  plotd: str,
  hm: str = "H3K27me3",
  binSize: str = "10000",
  cormthd: str = "pearson",
  colorMap: str = "RdBu",
  plotNumber: bool = False,
  plotHeight: str = "20",
  plotWidth: str = "25"
)-> None:
  prefix = f"{hm}.{cormthd}.binSize{binSize}"
  raw_cmds: List[str] = [
    "plotCorrelation", "-in", npzfnm, "--corMethod", cormthd,
    "--skipZeros",
    "--plotTitle", prefix,
    "--whatToPlot", "heatmap",
    "--colorMap", colorMap,
    "--plotNumbers" if plotNumber else "",
    "-o",
    os.path.join(plotd, hm, f"Heatmap_{prefix}.pdf"),
    "--removeOutliers",
    "--outFileCorMatrix",
    os.path.join(plotd, hm, f"Heatmap_{prefix}.corMatrix.tab"),
    "--plotHeight", plotHeight, "--plotWidth", plotWidth
  ]
  cmds = [e for e in raw_cmds if len(e) > 0]
  subprocess.run(cmds)
