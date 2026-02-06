import os
import sys
from typing import List, Dict, Tuple, Any
from collections import OrderedDict
import re
from multiprocessing import Pool
import functools
import logging


import pandas as pd
import pysam

from sinto import filterbarcodes
from sinto import utils

import pyprojroot
proj_dir = str(pyprojroot.here())
py_dir = f"{proj_dir}/package/python"
sys.path.insert(0, py_dir)
from mybam import get_barcodes_from_intervals_async
from myutils import reimport, flatten_list, unique_list
work_dir = os.path.join(proj_dir, "02.track")
outdir = os.path.join(work_dir, "out")
out_bam_dir = os.path.join(outdir, "cell_bams")
indexdir = os.path.join(outdir, "index_dir")

# * snakemake
cell_meta_fnm: str = os.path.join(
    work_dir, "src/main/resource/cellMeta.20240321.csv")
sublib_fnm: str = os.path.join(
    proj_dir, "meta/amb_pairedtag_sublibrary_v2.csv")
meta_sublib: pd.DataFrame = pd.read_csv(sublib_fnm,
                                        sep = ",", header = 0,
                                        names = None)
sublib_DNAs = meta_sublib['DNA']
sublib_DNA_bams = meta_sublib['DNA_bam']
sublib_ids = meta_sublib['sublibraryid']
d2bam: Dict[str, str] = {
    d: b for d, b in zip(sublib_DNAs, sublib_DNA_bams)}
d2sublib: Dict[str, str] = {
    d: s for d, s in zip(sublib_DNAs, sublib_ids)}
def get_index_fnm(sublib: str) -> str:
    return os.path.join(indexdir, f"{sublib}.index")
def get_bam_fnm(sublib:str) -> str:
    return d2bam[sublib]
def get_sublib_id(sublib: str) -> str:
    return d2sublib[sublib]

nthread: int = 2
overwrite: int = 0

# * failed samples
failed_sublibs: List[str] = ["ZW315",
"ZW317",
"ZW287",
"ZW319",
"ZW321",
"ZW735",
"ZW278",
"ZW280",
"ZW316",
"ZW318",
"ZW288",
"ZW320",
"ZW322",
"ZW277",
"ZW279"]

sublib = "ZW317"
sublib_id = get_sublib_id(sublib)
# * functions
def get_barcodes_from_bam_async(
        bam_fnm:str,
        index_fnm: str,
        barcode_pattern: str = ':\w\w:\w\w:\w\w:',
        nproc: int = 2) -> List[str]:
    with pysam.AlignmentFile(
        bam_fnm, index_filename = index_fnm, mode = "rb") as abam:
        intervals: Dict[int, List[Tuple[str, int, int]]] = utils.chunk_bam(
            abam, nproc, unmapped = True)
    p = Pool(nproc)
    barcodes = p.map_async(
        functools.partial(get_barcodes_from_intervals_async,
                          bam_fnm = bam_fnm,
                          index_fnm = index_fnm,
                          barcode_pattern = barcode_pattern),
        intervals.values()
    ).get(9999999)
    barcodes = flatten_list(barcodes)
    return(barcodes)

def write_oneread(oneread,
                  sublib_id,
                  barcodes: set,
                  outfhandles: Dict[str, Any],
                  logger: None|logging.Logger = None,
                  barcode_pattern: str = ':\w\w:\w\w:\w\w:',
                  verbose: bool = False) -> None:
    """
    1. get barcode, then get outfnm
    2. update_query_name by adding sublib_id
    3. write to the file
    """
    b = re.search(barcode_pattern,
                  oneread.query_name).group().strip(":")
    b2 = f"{sublib_id}:{b}"
    if b2 not in barcodes:
        if verbose:
            l = f"{b2} is not in the barcode list, ignored."
            if logger is None:
                print(l)
            else:
                logger.info(l)
        return None
    q = re.sub(
        barcode_pattern, f":{b2}:", oneread.query_name)

    if verbose:
        l = f"query name from {oneread.query_name} to {q}."
        if logger is None:
            print(l)
        else:
            logger.info(l)
    oneread.query_name = q
    try:
        outfhandles[b2].write(oneread)
    except KeyError:
        l = f"{b2} is not found in outfhandles."
        if logger is None:
            print(l)
        else:
            logger.info(l)
    finally:
        return None

# * main
barcode_pattern = ':\w\w:\w\w:\w\w:'
cell_meta = pd.read_csv(cell_meta_fnm,
                        sep = ",", header = 0,
                        names = None)
cell_meta = cell_meta[
    cell_meta['barcode'].str.contains(sublib_id, na = False)]
barcode2outdir: Dict[str, str] = OrderedDict(
    [(row['barcode'] , os.path.join(out_bam_dir, "{r}_{h}_{s}".format(
        r = row['brainregion'],
        h = row['modality'],
        s = row['rep']))) for _, row in cell_meta.iterrows()]
)

uoutdirs = list(set(barcode2outdir.values()))
# for d in uoutdirs:
#     os.makedirs(d, exist_ok = True)

barcode2outfnm: Dict[str, str] = OrderedDict(
    [k, os.path.join(v, f"{k}.srt.rmdup.bam")]
    for k, v in barcode2outdir.items()
)

barcodes: List[str] = get_barcodes_from_bam_async(
    bam_fnm = bam_fnm,
    index_fnm = index_fnm,
    barcode_pattern = barcode_pattern,
    nproc = nthread
)
ubarcodes = list(set(barcodes))

uubarcodes = [f"{sublib_id}:{u}" for u in ubarcodes]
uubarcodes = [u for u in uubarcodes if u in barcode2outfnm]
assert len(uubarcodes) == len(set(uubarcodes))

    
# get new head
with pysam.AlignmentFile(
    bam_fnm, index_filename = index_fnm, mode = "rb") as abam:
    header = abam.header.to_dict()
    validKeys = [x for x in ["HD", "SQ", "RG"] if x in header.keys()]
    newhead = dict((k, header[k]) for k in validKeys)

# set up file handlers
barcode2outfhs: Dict[str, Any] = OrderedDict(
    [(k,
      pysam.AlignmentFile(barcode2outfnm[k], 'wb', header = newhead))
     for k in uubarcodes]
)
with pysam.AlignmentFile(bam_fnm, index_filename = index_fnm,
                         mode = 'rb') as abam:
    # for r in abam.fetch('chr1',0 , 195471971):
    for r in abam.fetch():
        write_oneread(r, sublib_id = sublib_id,
                      barcodes = set(uubarcodes),
                      outfhandles = barcode2outfhs,
                      logger = None,
                      barcode_pattern = ':\w\w:\w\w:\w\w:',
                      verbose = True)
for _, fh in barcode2outfhs.items():
    fh.close()


