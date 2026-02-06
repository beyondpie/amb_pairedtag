.PHONY: ATAC_sa2DE

snakemake := /tscc/nfs/home/szu/miniforge3/bin/snakemake
ATAC_sa2DE: src/main/shell/sa2DE.Snakemake
	-mkdir -p build/$@
	cp $< build/$@/Snakefile
	cd build/$@ && ${snakemake} --snakefile Snakefile -R -c 20 \
    --keep-incomplete \
    --rerun-triggers mtime \
    --rerun-incomplete \
    --skip-script-cleanup \
    --keep-going \
    --jobs 100 \
    --profile=slurm
