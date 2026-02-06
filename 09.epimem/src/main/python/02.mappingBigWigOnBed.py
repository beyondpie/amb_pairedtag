import os
import sys
from typing import Tuple, Dict, List
import re
import numpy as np
import pandas as pd
import pyBigWig as pbw

# * meta
projd = "/projects/ps-renlab2/szu/projects/amb_pairedtag"
outd = os.path.join(projd, "09.epimem", "out")
ATAC_bwd = os.path.join(projd,
                        "data", "snATAC", "subclass_bigwig_bamCoverage")
ATAC_peakd = os.path.join(projd,
                          "data", "snATAC", "sa2.subclassv3.final.peak.srt")
Histones_bwd = os.path.join(projd, "data", "ptDNAbam", "bigwig")
Histones_peakd = os.path.join(projd,
                              "data", "pairedtag_peak", "subclass_peak")
# ext_size = [500, 1500, 2500]
ext_size = [2500]
# sc = "326_OPC_NN"
sc: str = sys.argv[1]
raw_from_m: str = sys.argv[2]
from_m: List[str] = raw_from_m.split(",")
raw_on_m: str = sys.argv[3]
on_m: List[str] = raw_on_m.split(",")
print(f"Map Bigwig Signals on Bed for: {sc}.")
print(f"Histone signals from: {raw_from_m}.")
print(f"Histone signals on: {raw_on_m}.")


# * function
def get_mergepeak_tag(m: str) -> str:
    if m == "H3K27ac":
        return "BestSPM"
    return "bedtoolmerge"


def get_smsize(m: str) -> int:
    if m == "H3K9me3":
        return 1000
    return 300


def load_histone_bw_and_bed(m: str, sc: str) -> Tuple[pbw.pyBigWig, pd.DataFrame]:
    mtag = get_mergepeak_tag(m)
    sm = get_smsize(m)
    bw = pbw.open(os.path.join(Histones_bwd, f"{sc}.{m}.e100.bs100.sm{sm}.bw"))
    beds = pd.read_csv(
        os.path.join(Histones_peakd, f"{sc}-{m}.{mtag}.peak"), sep="\t", header=0
    )
    return (bw, beds)


def load_ATAC_bw_and_bed(sc: str) -> Tuple[pbw.pyBigWig, pd.DataFrame]:
    bw = pbw.open(
        os.path.join(ATAC_bwd, f"{sc}.ATAC.e100.bs100.sm300.bw"))
    msc = re.sub("^\d+_", "", sc)
    beds = pd.read_csv(os.path.join(ATAC_peakd, f"{msc}.bed"), header=None, sep="\t")
    beds.columns = ["chrom", "startFrom", "endTo", "name"]
    return (bw, beds)

def load_histone_bw(m: str, sc: str) -> pbw.pyBigWig:
    sm = get_smsize(m)
    bw = pbw.open(os.path.join(Histones_bwd, f"{sc}.{m}.e100.bs100.sm{sm}.bw"))
    return bw

def load_histone_bed(m: str, sc: str) -> pd.DataFrame:
    mtag = get_mergepeak_tag(m)
    beds = pd.read_csv(
        os.path.join(Histones_peakd, f"{sc}-{m}.{mtag}.peak"), sep="\t", header=0
    )
    return beds

def load_ATAC_bw(sc: str) -> pbw.pyBigWig:
    bw = pbw.open(
        os.path.join(ATAC_bwd, f"{sc}.ATAC.e100.bs100.sm300.bw"))
    return bw

def load_ATAC_bed(sc: str) -> pd.DataFrame:
    msc = re.sub("^\d+_", "", sc)
    beds = pd.read_csv(os.path.join(ATAC_peakd, f"{msc}.bed"), header=None, sep="\t")
    beds.columns = ["chrom", "startFrom", "endTo", "name"]
    return beds

def extend_narrow_peak(beds: pd.DataFrame,
                       extend_size: int = 500,
                       default_relaSummit: int = 249) -> pd.DataFrame:
    """
    Extend narrow peaks around their summits:

        summit = startFrom + relaSummit
        [summit - extend_size, summit + extend_size)
        If relaSummit is not in beds.columns, use default_relaSummit,
        which is especially used for ATAC-seq.
    """
    if "relaSummit" in beds.columns:
        summit = beds.startFrom + beds.relaSummit
    else:
        summit = beds.startFrom + default_relaSummit
    r = pd.DataFrame(
        data={
            "chrom": beds.chrom,
            "startFrom": summit - extend_size,
            "endTo": summit + extend_size
        }
    )
    return r

def extend_broad_peak(
        beds: pd.DataFrame, extend_size: int = 500) -> pd.DataFrame:
    r = pd.DataFrame(
        data={
            "chrom": beds.chrom,
            "startFrom": beds.startFrom - extend_size,
            "endTo": beds.endTo + extend_size
        }
    )
    return r

def cal_biwig_on_a_region(bw: pbw.pyBigWig, dfrow: pd.Series) -> float:
    r = bw.stats(
        dfrow["chrom"], dfrow["startFrom"], dfrow["endTo"], type="mean", exact=True
    )
    return r


def map_bigwig_on_regions(bw: pbw.pyBigWig, region: pd.DataFrame) -> np.ndarray:
    return region.apply(
        lambda x: cal_biwig_on_a_region(bw, x)[0], axis=1).to_numpy()


def get_outlier(
    x: np.ndarray, up_quantile: float = 0.995, fold_std: float = 2.0
) -> float:
    q = np.quantile(x, up_quantile)
    print(f"max: {x.max()}; quantile {up_quantile}: {q}.")
    v = x[x <= q]
    m = np.mean(v)
    sd = np.std(v)
    thres = m + fold_std * sd
    print(f"mean: {m}; std: {sd}; threshold: {thres}.")
    print(f"under threshold {thres}, {(x >= thres).sum()} outliers.")
    return thres


# * main
# load subclass peaks on Histone modifications and snATAC
# H3K27me3_bw, H3K27me3_bed = load_histone_bw_and_bed("H3K27me3", sc)
# H3K27ac_bw, H3K27ac_bed = load_histone_bw_and_bed("H3K27ac", sc)
# H3K9me3_bw, H3K9me3_bed = load_histone_bw_and_bed("H3K9me3", sc)
# H3K4me1_bw, H3K4me1_bed = load_histone_bw_and_bed("H3K4me1", sc)
# ATAC_bw, ATAC_bed = load_ATAC_bw_and_bed(sc)

# all_bws: Dict[str, pbw.pyBigWig] = {
#     "ATAC": ATAC_bw,
#     "H3K27me3": H3K27me3_bw,
#     "H3K27ac": H3K27ac_bw,
#     "H3K9me3": H3K9me3_bw,
#     "H3K4me1": H3K4me1_bw
# }

# all_beds: Dict[str, pd.DataFrame] = {
#     "ATAC": ATAC_bed,
#     "H3K27me3": H3K27me3_bed,
#     "H3K27ac": H3K27ac_bed,
#     "H3K9me3": H3K9me3_bed,
#     "H3K4me1": H3K4me1_bed
# }

# BUG: Ignore ATAC here
all_m = list(set(from_m + on_m))
all_bws = {k: load_histone_bw(k, sc) for k in all_m}
all_beds = {k: load_histone_bed(k, sc) for k in all_m}


# calculate the signals
def run_(sc: str, from_m: str, to_m: str, ext_size: int, skip: bool = False) -> None:
    print(f"{sc}: {from_m} onto {to_m} with extension {ext_size}...")
    outfnm = os.path.join(
        outd, f"{from_m}_on_{to_m}",
        f"{sc}_{from_m}-on-{to_m}_e{ext_size}.tsv")
    if os.path.exists(outfnm) and skip:
        print(f"{outfnm} exists, and skipt it.")
        return None
    if not os.path.exists(os.path.dirname(outfnm)):
        os.makedirs(os.path.dirname(outfnm), exist_ok=True)
    cur_bw = all_bws[from_m]
    cur_bed = all_beds[to_m]
    if to_m in ["ATAC", "H3K27ac"]:
        ext_bed = extend_narrow_peak(cur_bed, ext_size)
    else:
        ext_bed = extend_broad_peak(cur_bed, ext_size)
    r = map_bigwig_on_regions(cur_bw, ext_bed)
    out_bed = cur_bed.copy()
    out_bed.insert(out_bed.shape[1], f"{from_m}_e{ext_size}", r)
    out_bed.to_csv(outfnm, sep="\t", header=True, index=False)


for i in from_m:
    for j in on_m:
        for e in ext_size:
            run_(sc, i, j, e)

# for j in ["H3K4me1", "ATAC", "H3K27ac"]:
#     for i in ["H3K27me3", "H3K9me3"]:
#         for e in ext_size:
#             run_(sc, j, i, e)
