import os
import sys
from typing import List, Dict, Tuple, Any
from collections import OrderedDict
import string
import random
import re
import logging
from multiprocessing import Pool
import functools


import pandas as pd
import pysam

from sinto import filterbarcodes
from sinto import utils

import pyprojroot
proj_dir = str(pyprojroot.here())
py_dir = f"{proj_dir}/package/python"
sys.path.insert(0, py_dir)
from mylog import StreamToLogger, set_file_logger
from mybam import get_barcodes_from_intervals_async
from myutils import reimport, flatten_list, unique_list


# * snakemake
logfnm: str = snakemake.log[0]
bam_fnm: str = snakemake.input['bam_fnm']
index_fnm: str = snakemake.input['index_fnm']
cell_meta_fnm: str = snakemake.input['cell_meta_fnm']
sublib_id: str = snakemake.params['sublib_id']
out_bam_dir: str = snakemake.params['out_bam_dir']
nthread: int = snakemake.threads
debug: int = snakemake.params['debug']
overwrite: int = snakemake.params['overwrite']

# * set debug status
if debug > 0:
    print("Under debug mode ...")
    out_dir = os.path.join(proj_dir, "02.track", "test")
    logfnm = os.path.join(out_dir, "test.02.track.log")
    out_bam_dir = os.path.join(out_dir, "bams")
    for d in [out_dir, out_bam_dir]:
        os.makedirs(d, exist_ok = True)
    
    # column 9
    sexs = ['MaleA', 'MaleB', 'FemaleA', 'FemaleB']
    # column 4
    regions = ['AMY', 'CPU', 'ERC', 'HCa',
            'HCp', 'HYP', 'NAC', 'PFC', 'VTA_SnR']

    out_bam_group_dirs: List[str] = [os.path.join(out_bam_dir, f"{r}_{s}_RNA")
                        for r in regions
                        for s in sexs]
    for d in out_bam_group_dirs:
        os.makedirs(d, exist_ok = True)
    nthread = 2
    cell_meta_fnm = os.path.join(proj_dir,
                                 "meta", "pairedtag.meta.v1.csv")
    sublib_id = 'A01'
    bam_fnm = os.path.join(proj_dir, "02.track",
                           "src/test/resource",
                           "ZW347_mm10_sorted_rmdup.bam")
    index_fnm = os.path.join(os.path.dirname(bam_fnm),
                             "ZW247.index")
    verbose = False
else:
    print("Under normal mode ...")
    verbose = False

# * set logger
logger = set_file_logger(logfnm, name = "run_macs2")
## works in Linux, but have some troubles in mac
sys.stdout = StreamToLogger(logger = logger, level = logging.INFO)
sys.stderr = StreamToLogger(logger = logger, level = logging.ERROR)

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
    [(row['barcode'] , os.path.join(out_bam_dir, "{r}_{s}_RNA".format(
        r = row['brainregion'],
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
logger.info(f"{len(barcodes)} from {sublib_id}.")
ubarcodes = list(set(barcodes))
logger.info(f"unique {len(ubarcodes)} from {sublib_id}.")

uubarcodes = [f"{sublib_id}:{u}" for u in ubarcodes]
uubarcodes = [u for u in uubarcodes if u in barcode2outfnm]
assert len(uubarcodes) == len(set(uubarcodes))
logger.info(f"{len(uubarcodes)} passed QC.")

if overwrite >= 1:
    logger.info("Overwrite existed bam files.")
else:
    logger.info("Only bam files not found will be generated.")
    uubarcodes = [u for u in uubarcodes
                  if not os.path.exists(barcode2outfnm[u])]
    if len(uubarcodes) < 1:
        logger.info("All the barcodes have the corresponding files.")
        sys.exit(status = 0)
    logger.info(f"{len(uubarcodes)} will be processed.")
    
# get new head
with pysam.AlignmentFile(
    bam_fnm, index_filename = index_fnm, mode = "rb") as abam:
    header = abam.header.to_dict()
    validKeys = [x for x in ["HD", "SQ", "RG"] if x in header.keys()]
    newhead = dict((k, header[k]) for k in validKeys)

# set up file handlers
logger.info(
    f"create {len(uubarcodes)} out bam filehandles for {sublib_id}.")
barcode2outfhs: Dict[str, Any] = OrderedDict(
    [(k,
      pysam.AlignmentFile(barcode2outfnm[k], 'wb', header = newhead))
     for k in uubarcodes]
)
logger.info(
    f"start to split bam for {sublib_id} with {len(uubarcodes)} barcodes.")
with pysam.AlignmentFile(bam_fnm, index_filename = index_fnm,
                         mode = 'rb') as abam:
    # for r in abam.fetch('chr1',0 , 195471971):
    for r in abam.fetch():
        write_oneread(r, sublib_id = sublib_id,
                      barcodes = set(uubarcodes),
                      outfhandles = barcode2outfhs,
                      logger = logger,
                      barcode_pattern = ':\w\w:\w\w:\w\w:',
                      verbose = verbose)
logger.info(f"end of splitting bam for {sublib_id}.")
for _, fh in barcode2outfhs.items():
    fh.close()
logger.info(f"close all the handles in {sublib_id}.")


