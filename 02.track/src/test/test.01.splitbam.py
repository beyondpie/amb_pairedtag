import os
import sys
import string
import random
import re

import pyprojroot
import sinto
from sinto import filterbarcodes
from sinto import utils
import pysam

# * configs
proj_dir = str(pyprojroot.here())
barcode_pattern = ':\w\w:\w\w:\w\w:'

# * load bam file
workdir = os.path.join(pyprojroot.here(), "02.track")
workdir = os.path.join("/Users/szu/git-recipes",
                       "mouseBrainAtlas", "amb_pairedtag",
                       "02.track")
out_dir = os.path.join(proj_dir, "02.track", "test_out")
logfnm = os.path.join(out_dir, "test.02.track.log")
bamdir = os.path.join(out_dir, "bams")
for d in [out_dir, bamdir]:
    os.makedirs(d, exist_ok = True)

hms = ['H3K27ac', 'H3K27me3', 'H3K4me1', 'H3K9me3']
# column 9
sexs = ['MaleA', 'MaleB', 'FemaleA', 'FemaleB']
# column 4
regions = ['AMY', 'CPU', 'ERC', 'HCa',
            'HCp', 'HYP', 'NAC', 'PFC', 'VTA_SnR']

out_bam_group_dirs: List[str] = [os.path.join(bamdir, f"{r}_{h}_{s}")
                    for r in regions
                    for h in hms
                    for s in sexs]
for d in out_bam_group_dirs:
    os.makedirs(d, exist_ok = True)
sublib_id = 'A01'
nthread = 2
cell_meta_fnm = os.path.join(proj_dir,
                                "meta", "pairedtag.meta.v1.csv")
bam_fnm = os.path.join(proj_dir, "02.track",
                        "src/test/resource",
                        "ZW347_mm10_sorted_rmdup.bam")

# * main
abam = pysam.AlignmentFile(
    os.path.join(workdir,
                 "src/test",
                 "resource", "ZW347_mm10_sorted_rmdup.bam"),
    mode = "rb"
)

ident = "".join(
    random.choice(
        string.ascii_uppercase + string.digits) for _ in range(6)
)
header = abam.header.to_dict()
validKeys = [x for x in ["HD", "SQ", "RG"] if x in header.keys()]
        newhead = dict((k, header[k]) for k in validKeys)

exp_reads = [i for i in abam.head(10)]
one_read = exp_reads[0]
re.search(barcode_pattern, one_read.query_name).group().strip(":")


# TODO:
# 1. test out bam is right (compare with samtools)
# 2. test out bam is OK to merge (if the header is OK)
