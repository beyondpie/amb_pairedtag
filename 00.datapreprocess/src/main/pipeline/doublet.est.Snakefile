envvars:
    "PATH"

import os
import pandas as pd
from typing import List, Dict, Tuple
import pyprojroot

# * configs
nthread = config["nthread"]
debug = config["debug"]
genome = config["genome"]

proj_dir = pyprojroot.here()
# FIXME: name has 882 unique values, while 889 rows in total
if debug < 1:
    meta_sublib = pd.read_csv(
        os.path.join(proj_dir, "meta", "amb_pairedtag_sublibrary.csv"),
        sep = ",",
        header = 0,
        names = None)
else:
    print("Debug mode.")
    meta_sublib = pd.read_csv(
        os.path.join(proj_dir, "00.datapreprocess",
                     "src", "test", "resource",
                     "amb_pairedtag_sublibrary.csv"),
        sep = ",",
        header = 0,
        names = None)

fastq_DNA: Dict[str, str] = {
    row['DNA'] : os.path.join(
        row['path'], "02.trimmed",
        f"{row['DNA']}_BC_cov_trimmed.fq.gz")
    for index, row in meta_sublib.iterrows()
}
sublib_DNA_ids = fastq_DNA.keys()

fastq_RNA: Dict[str, str] = {
    row['RNA'] : os.path.join(
        row['path'], "02.trimmed",
        f"{row['RNA']}_BC_cov_trimmed_trimmed.fq.gz")
    for index, row in meta_sublib.iterrows()
}
sublib_RNA_ids = fastq_RNA.keys()
# make sure files are exist.
for d in list(fastq_DNA.values()) + list(fastq_RNA.values()):
    if not os.path.exists(d):
        raise FileNotFoundError(d)
    # else:
    #     print(f"{d} exists.")

print(f"get {len(fastq_DNA)} sublibraries.")
if genome == "hg38":
    print("Under genome hg38, ignore DNA and RNA mapping from Exp14.")
# ignore exp14 since no Hela cell spike-in
    Ids_exp14: List[Tuple[str, str]] = [(row['DNA'], row['RNA'])
                                        for _, row in meta_sublib.iterrows()
                                        if row['exp'] == 'exp14']
    DNA_ids_exp14 = [t[0] for t in Ids_exp14]
    RNA_ids_exp14 = [t[1] for t in Ids_exp14]
    sublib_DNA_ids = [i for i in sublib_DNA_ids if i not in DNA_ids_exp14]
    print(f"{len(sublib_DNA_ids)} DNA sublibraries left.")
    sublib_RNA_ids = [i for i in sublib_RNA_ids if i not in RNA_ids_exp14]
    print(f"{len(sublib_RNA_ids)} RNA sublibraries left.")

# ignore ZW1428 since the fastq file has unknow issue
# when performing mapping: gzip: /projects/ps-renlab2/zhw063/99.MouseBrainPairedTag/08.MouseBrainExp8/02.MouseBrainPT_Exp8B_20230911/02.trimmed/ZW1428_BC_cov_trimmed_trimmed.fq.gz: unexpected end of file
# EXITING because of FATAL ERROR in reads input: quality string length is not equal to sequence length
# @LH00320:28:22C7N7LT3:7:2113:13750:14314:C1:CD:16:CCAGACGGGT
# GTATATAGATATACAGATGGATGGATGGATAGATGGACAAGTGGGTGTAACTAAAGTAATCCAGTCACAGA
# FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
# SOLUTION: fix your fastq file
# Dec 23 14:52:49 ...... FATAL ERROR, exiting

print("Remove RNA sublib: ZW1428 ... ")
sublib_RNA_ids = [i for i in sublib_RNA_ids if i not in ["ZW1428"]]
print(f"{len(sublib_RNA_ids)} RNA sublibraries left.")
    
work_dir = os.path.join(proj_dir,"00.datapreprocess")
outdir = os.path.join(work_dir,"out")
flagdir = os.path.join(outdir, "flag")
logdir = os.path.join(outdir, "log")
mapdir = os.path.join(outdir, "mapping")
mtxdir = os.path.join(outdir, "mtx")
mtxdir_DNA_hg38 = os.path.join(mtxdir, "mtx_DNA_hg38")
mtxdir_RNA_hg38 = os.path.join(mtxdir, "mtx_RNA_hg38")
mtxdir_DNA_mm10 = os.path.join(mtxdir, "mtx_DNA_mm10")
mtxdir_RNA_mm10 = os.path.join(mtxdir, "mtx_RNA_mm10")

for d in [outdir, flagdir, logdir,
          mapdir, mtxdir, mtxdir_DNA_hg38,
          mtxdir_RNA_hg38]:
    os.makedirs(d, exist_ok = True)

genome_ref_dir = "/projects/ps-renlab/y2xie/projects/genome_ref"
mm10_gtf = os.path.join(genome_ref_dir, "gencode.vM25.annotation.gtf")
mm10_STAR_index = os.path.join(genome_ref_dir, "mm10", "star")
# mm10_STAR = os.path.join("/home/y2xie/package/STAR-2.7.1a",
                         # "bin/Linux_x86_64", "STAR")
mm10_bowtie2_index = "/projects/ps-renlab/share/bowtie2_indexes/mm10"
mm10_5k = os.path.join(genome_ref_dir, "Paired-Tag",
                       "mm10/mm10.bin5k.txt")
mm10_rna = os.path.join(genome_ref_dir, "Paired-Tag", "mm10",
                        "mm10.PairedTag.txt")

hg38_gtf = os.path.join(genome_ref_dir, "gencode.vH35.annotation.gtf")
hg38_STAR_index = "/projects/ps-renlab/share/STAR_indices/hg38"
hg38_STAR = os.path.join("/home/zhw063/ps-renlab/software/STAR-2.5.3a/",
                         "bin/Linux_x86_64/STAR")
hg38_bowtie2_index = "/projects/ps-renlab/share/bowtie2_indexes/hg38"

hg38_5k = os.path.join(genome_ref_dir, "Paired-Tag",
                       "hg38/hg38.bin5k.txt")
hg38_rna = os.path.join("/projects/ps-renlab/zhw063/software",
                        "Paired-Tag-master",
                        "references", "hg38.RNA.txt")

PT48 = os.path.join("/projects/ps-renlab/zhw063/software",
                    "Paired-Tag-master",
                    "references", "PairedTag48_384")
mtxdir_RNA = mtxdir_RNA_hg38 if genome == "hg38" else mtxdir_RNA_mm10
mtxdir_DNA = mtxdir_DNA_hg38 if genome == "hg38" else mtxdir_DNA_mm10
bin_5k = hg38_5k if genome == "hg38" else mm10_5k
gene = hg38_rna if genome == "hg38" else mm10_rna

# * tools
bowtie2 = "/home/szu/miniforge3/bin/bowtie2"
reachtools = os.path.join("/projects/ps-renlab2/szu/softwares",
                          "Paired-Tag", "reachtools", "reachtools")
# mediator
conda_dir = "/home/szu/miniforge3"
samtools = os.path.join(conda_dir, "bin", "samtools")
bowtie2 = os.path.join(conda_dir, "bin", "bowtie2")

# my_STAR = "/home/szu/miniforge3/bin/STAR"
# hg38_STAR = os.path.join("/home/zhw063/ps-renlab/software/STAR-2.5.3a/",
#                          "bin/Linux_x86_64/STAR")
# STAR = hg38_STAR if genome == "hg38" else my_STAR
# now version 2.5.3a for hg38
STAR = "/home/szu/miniforge3/bin/STAR"

# * functions
def get_DNA_fastq(wildcards) -> str:
    return fastq_DNA[wildcards.s]
def get_RNA_fastq(wildcards) -> str:
    return fastq_RNA[wildcards.s]

def get_DNA_bin(wildcards) -> str:
    if wildcards.g == "mm10":
        return mm10_5k
    return hg38_5k

DNA_ref_index = hg38_bowtie2_index if genome == "hg38" else mm10_bowtie2_index
RNA_ref_index = hg38_STAR_index if genome == "hg38" else mm10_STAR_index

def get_mtx_DNA(wildcards) -> str:
    if wildcards.g == "mm10":
        return  mtxdir_DNA_mm10
    return mtxdir_DNA_hg38
def get_mtx_RNA(wildcards) -> str:
    if wildcards.g == "mm10":
        return  mtxdir_RNA_mm10
    return mtxdir_RNA_hg38


rule all:
    input:
        expand("{f}/{s}_{g}_DNA_mapping.done",
               f = flagdir, s = sublib_DNA_ids, g = genome),
        expand("{f}/{s}_{g}_RNA_mapping.done",
               f = flagdir, s = sublib_RNA_ids, g = genome),
        expand("{f}/{s}_{g}_DNA_mtx.done",
               f = flagdir, s = sublib_DNA_ids, g = genome),
        expand("{f}/{s}_{g}_RNA_mtx.done",
               f = flagdir, s = sublib_RNA_ids, g = genome)
rule DNA_mapping:
    input:
        fastq_gz = get_DNA_fastq
    output:
        bam = f"{mapdir}/{{s}}_{genome}_sorted.bam",
        flag = touch(f"{flagdir}/{{s}}_{genome}_DNA_mapping.done")
    params:
        ref_index = DNA_ref_index,
        sam = f"{mapdir}/{{s}}_{genome}.sam"
    threads: nthread
    log:
        f"{logdir}/{{s}}_{genome}_DNA_mapping.log"
    shell:
        """
        echo "run DNA bowtie2 ..."
        {bowtie2} -x {params.ref_index} -U {input.fastq_gz} \
                  --no-unal -p 16 -S {params.sam}
        echo "run samtools ..."
        {samtools} sort -@ 16 -T {mapdir} -o {output.bam} {params.sam}

        echo "remove sam file ..."
        rm -f {params.sam}
        """
        
rule RNA_mapping:
    input:
        fastq_gz = get_RNA_fastq
    output:
        bam = f"{mapdir}/{{s}}_{genome}_sorted.bam",
        flag = touch(f"{flagdir}/{{s}}_{genome}_RNA_mapping.done")
    params:
        ref_index = RNA_ref_index,
        outprefix = f"{mapdir}/{{s}}_{genome}_",
        bam = f"{mapdir}/{{s}}_{genome}_clean.bam"
    threads: nthread
    log:
        f"{logdir}/{{s}}_{genome}_RNA_mapping.log"
    shell:
        """
        echo "run RNA STAR ..."
        {STAR} --runThreadN {nthread} --genomeDir {params.ref_index} \
               --readFilesIn {input.fastq_gz} \
               --readFilesCommand zcat \
               --outFileNamePrefix {params.outprefix} \
               --outSAMtype BAM Unsorted --quantMode GeneCounts
        echo "run samtools to clean ..."
        {samtools} view -h -F 256 \
           {mapdir}/{wildcards.s}_{genome}_Aligned.out.bam -b \
           > {params.bam}
        echo "run samtools to sort ..."
        {samtools} sort -@ 16 -T {mapdir} -o {output.bam} {params.bam}
        
        echo "remove bam file ..."
        rm -f {params.bam}
        """

rule DNA_bam2mtx:
    input:
        bam = f"{mapdir}/{{s}}_{genome}_sorted.bam"
    output:
        flag = touch(f"{flagdir}/{{s}}_{genome}_DNA_mtx.done")
    params:
        bam = f"{mapdir}/{{s}}_{genome}_sorted_rmdup.bam",
        mtx2 = f"{mapdir}/{{s}}_{genome}_sorted_rmdup_mtx2"
    log:
        f"{logdir}/{{s}}_{genome}_DNA_mtx.log"
    shell:
        """
        echo "remove duplicates in bam file ..."
        {reachtools} rmdup2 {input.bam}
        echo "get gene matrix ..."
        {reachtools} bam2Mtx2 {params.bam} {bin_5k}
        echo "move mtx2 to mtx dir ..."
        mv {params.mtx2} {mtxdir_DNA}
        """

rule RNA_bam2mtx:
    input:
        bam = f"{mapdir}/{{s}}_{genome}_sorted.bam",
    output:
        flag = touch(f"{flagdir}/{{s}}_{genome}_RNA_mtx.done")
    params:
        bam = f"{mapdir}/{{s}}_{genome}_sorted_rmdup.bam",
        mtx2 = f"{mapdir}/{{s}}_{genome}_sorted_rmdup_mtx2"
    log:
        f"{logdir}/{{s}}_{genome}_RNA_mtx.log"
    shell:
        """
        echo "remove duplicates in bam file ..."
        {reachtools} rmdup2 {input.bam}
        echo "get gene matrix ..."
        {reachtools} bam2Mtx2 {params.bam} {gene}
        echo "move mtx2 to mtx dir ..."
        mv {params.mtx2} {mtxdir_RNA}
        """
