"""
Get pseudobulk-level RNA or DNA-related bam file based on
- region
- histone modification
- sex
- replication
- group
- OR Directly a group of barcodes
- omic: DNA or RNA

Currently, we only support
- one region or all regions
- one group or all groups

Dependencies:
- pysam
- pandas


Usage:

- to get bam from a specific region, histone modification, rep
  - python 03.mergebam.py -a PFC -i H3K27ac  -r A -o PFC_H3K27ac_repA.bam

- to get bam from a group of barcode in a file
  - The file has no column name, and each line is a barcode
    python 03.mergebam.py -b [your_barcode_file] -o [outname].bam

- to get bam from a group of barcode in a file for RNA omic
  - The file has no column name, and each line is a barcode
    python 03.mergebam.py -b [your_barcode_file] -o [outname].bam --omic RNA

- to get bam from a group of barcode in a file for RNA under tscc2 system
  - The file has no column name, and each line is a barcode
    python 03.mergebam.py -b [your_barcode_file] -o [outname].bam \
                          --omic RNA --system tscc2



Note:
- If having too many files open Error, run `ulimit -n 50000` firstly.
- The output is not sorted by coordinated.
  - To get the bigwig file:
    1. Use samtools sort and then index the sorted bam file.
    2. Use bamCoverage from deepTools to generate the bigwig file.
"""

import os
import tempfile
import argparse
import string
import random
from typing import List
import yaml

import pandas as pd
import pysam

regions: List[str] = ["AMY", "CPU", "ERC",
                      "HCa", "HCp", "HYP",
                      "NAC", "PFC", "VTA_SnR",
                      'NA']

rdm_str = ''.join(random.choices(
    string.ascii_uppercase + string.digits, k=6))
default_outfnm = os.path.join(".", f"mergedbam_{rdm_str}.bam")

parser = argparse.ArgumentParser(
    prog="mergebam",
    description="merge cell-level bam files."
)

parser.add_argument('-a', '--region',
                    choices=regions,
                    dest='region',
                    type=str,
                    default='NA',
                    help="mouse anatomic region")
parser.add_argument('-i', '--histone',
                    choices=['H3K27ac', 'H3K27me3',
                             'H3K9me3', 'H3K4me1',
                             'NA'],
                    dest='histone',
                    type=str,
                    default='NA',
                    help="histone modification type")
parser.add_argument('-s', '--sex',
                    choices=['Male', 'Female', 'NA'],
                    dest='sex',
                    default='NA',
                    type=str,
                    help="sex info")
parser.add_argument('-r', '--rep',
                    choices=['A', 'B', 'NA'],
                    dest='rep',
                    default='NA',
                    type=str,
                    help="replicate info")
parser.add_argument('-g', '--group',
                    dest="group",
                    type=str,
                    default='NA',
                    help="group of cells")
parser.add_argument('-c', '--groupcol',
                    dest='col',
                    type=str,
                    default='NA',
                    help="column for the group info")
parser.add_argument('-o', '--out',
                    dest='outpath',
                    type=str,
                    default=default_outfnm,
                    help="outfile path")
parser.add_argument('-b', '--barcode',
                    dest='barcode',
                    type=str,
                    default='NA',
                    help="file of barcodes for merging")
parser.add_argument('--omic',
                    choices=['RNA', 'DNA'],
                    dest="omic",
                    type=str, default='DNA',
                    help='DNA or RNA, default DNA')
# * functions
def check_args(args, metafnm: str) -> None:
    if args.barcode != 'NA' and (not os.path.exists(args.barcode)):
        raise FileNotFoundError(
            f"barcode file of {args.barcode} NOT exist.")
    meta = pd.read_csv(metafnm, header=0, sep=",")
    if args.col != 'NA' and (not args.col in meta.columns):
        raise RuntimeError(
            f"col of {args.col} NOT exist in the meta file.")
    if args.group != 'NA' and (not args.group in meta[args.col].unique()):
        raise RuntimeError(
            f"{args.group} not found in col of {args.col}."
        )
    outdir = os.path.dirname(args.outpath)
    if (len(outdir) > 0) and (not os.path.exists(outdir)):
        print(f"dir of {args.outpath} does not exist.")
        os.makedirs(outdir, exist_ok=True)
        print(f"create the dir {outdir}.")
    if (args.omic == "RNA") and (args.histone != 'NA'):
        raise RuntimeError(f"No Histones needed for RNA.")
    return None

def get_projd_from_git() -> str:
    curd = os.getcwd()
    while(not os.path.isdir(os.path.join(curd, ".git"))):
        curd = os.path.dirname(curd)
    return curd

def get_barcodes(args,
                 barcode_meta: pd.DataFrame) -> List[str]:
    if args.histone == 'NA':
        print(f"No filter on histone since it's {args.histone}.")
        m1 = barcode_meta
    else:
        print(f"filter histone: {args.histone}.")
        m1 = barcode_meta[barcode_meta['modality'] == args.histone]
    if args.sex == 'NA':
        print(f"No filter on sex since it's {args.sex}.")
        m2 = m1
    else:
        print(f"filter sex: {args.sex}.")
        m2 = m1[m1['sex'] == args.sex]
    if args.rep == 'NA':
        print(f"No filter on replicate since it's {args.rep}.")
        m3 = m2
    else:
        print(f"filter replicate: {args.rep}")
        m3 = m2[m2['rep'].str.endswith(args.rep)]
    if args.group == 'NA':
        print(f"No filter on group since it's {args.group}.")
        m4 = m3
    else:
        print(f"filter group: {args.group}.")
        m4 = m3[m3[args.col] == str(args.group)]
    if args.region == 'NA':
        print(f"No filter on region since it's {args.region}.")
        m5 = m4
    else:
        print(f"filter region: {args.region}.")
        m5 = m4[m4['brainregion'] == args.region]
    print(f"{len(m5)} barcodes found.")
    return m5['barcode'].tolist()


def read_barcodes(args, barcode_meta: pd.DataFrame) -> List[str]:
    print(f"load barcodes from file: {args.barcode}.")
    raw_barcodes: List[str] = pd.read_csv(
        args.barcode, header=None).iloc[:, 0].to_list()
    print(f"got {len(raw_barcodes)} barcodes.")
    all_barcodes = set(barcode_meta['barcode'])
    r = [b for b in raw_barcodes if b in all_barcodes]
    print(f"{len(r)} barcodes in our record.")
    return r


def get_bam_files(barcodes: List[str],
                  barcode_meta: pd.DataFrame,
                  bamdir: str,
                  data_type: str = 'DNA') -> List[str]:
    if len(barcodes) < 1:
        raise RuntimeError("No barcodes for bam files.")
    if data_type not in ['DNA', 'RNA']:
        raise RuntimeError(f"data_type {data_type} should be RNA or DNA.")
    m: pd.DataFrame = barcode_meta.loc[barcodes]
    if data_type == 'DNA':
        bamfiles = [os.path.join(bamdir,
                                 f"{r['brainregion']}_{r['modality']}_{r['rep']}",
                                 f"{r['barcode']}.srt.rmdup.bam")
                    for _, r in m.iterrows()]
    else:
        bamfiles = [os.path.join(bamdir,
                                 f"{r['brainregion']}_{r['rep']}_RNA",
                                 f"{r['barcode']}.srt.rmdup.bam")
                    for _, r in m.iterrows()]
    return bamfiles


def merge_bam_file(bam_files: List[str], outfnm: str) -> None:
    for i in bam_files:
        if not os.path.exists(i):
            raise FileNotFoundError(f"{i} NOT exist.")
    tmp_bamfnms = tempfile.NamedTemporaryFile(
        delete=False, mode='w')
    tmp_bamfnms.write('\n'.join(bam_files))
    tmp_bamfnms.close()
    print(f"put bam_files to temp {tmp_bamfnms.name}.")
    print(f"start merging {len(bam_files)} bams to {outfnm}...")
    pysam.merge("-f", "-o", outfnm, "-b", tmp_bamfnms.name)
    os.unlink(tmp_bamfnms.name)
    print(f"merging bam is done.")


# * main
args = parser.parse_args()
# set metafnm and bamdir
projd: str = get_projd_from_git()
config_yaml = os.path.join(projd, "config.yaml")
with open(config_yaml, 'r') as f:
    config = yaml.safe_load(f)
data_projd = os.path.join(config['disk'], config['project'])

meta_local_fnm = config['pt_scmeta_path']
DNA_dir = config['pt_scDNAbam_path']
RNA_dir = config['pt_scRNAbam_path']

metafnm = os.path.join(projs, meta_local_fnm)
if not os.path.exists(metafnm):
    raise FileNotFoundError(f"cannot find {metafnm}.")
bam_local_dir = (RNA_dir, DNA_dir)[args.omic == "DNA"]
bamdir = os.path.join(projs, bam_local_dir)
print(f"Under omic of {args.omic}, bamdir is set to {bamdir}.")
if not os.path.exists(bamdir):
    raise FileNotFoundError(f"cannot find {bamdir}.")
# check arguments
check_args(args, metafnm)
barcode_meta = pd.read_csv(metafnm,
                           header=0, sep=',',
                           index_col=False)
barcode_meta.set_index('barcode', inplace=True, drop=False)

if args.barcode != 'NA':
    barcodes = read_barcodes(args, barcode_meta)
else:
    barcodes = get_barcodes(args, barcode_meta)
bam_files = get_bam_files(barcodes, barcode_meta, bamdir,
                          data_type = args.omic)
merge_bam_file(bam_files, outfnm=args.outpath)
print("Done.")
