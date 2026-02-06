snakemake := /home/szu/miniforge3/envs/sa2/bin/snakemake
#snakemake := /tscc/nfs/home/szu/miniforge3/bin/snakemake

.PHONY: test_doublet_est double_est run_zst
test_doublet_est: src/main/pipeline/doublet.est.Snakefile
	-rm -rf build/$@
	-mkdir -p build/$@
	cp $< build/$@/Snakefile
	cd build/$@ && ${snakemake} --snakefile Snakefile -R -c 10 \
    --config \
      nthread=4 \
      debug=1 \
      genome=hg38 \
    --keep-incomplete --rerun-triggers mtime \
    --skip-script-cleanup 


doublet_est: src/main/pipeline/doublet.est.Snakefile
	-mkdir -p build/$@
	cp $< build/$@/Snakefile
	cd build/$@ && ${snakemake} --snakefile Snakefile -R -c 50 \
    --config \
      nthread=4 \
      debug=0 \
      genome=hg38 \
    --keep-incomplete --rerun-triggers mtime \
    --skip-script-cleanup \
    --keep-going

run_zst: src/main/pipeline/fastq.gz2zst.Snakefile
	-mkdir -p build/$@
	cp $< build/$@/Snakefile
	cd build/$@ && ${snakemake} --snakefile Snakefile -R -c 101 \
    --keep-incomplete \
    --rerun-triggers mtime \
    --rerun-incomplete \
    --skip-script-cleanup \
    --keep-going \
    --jobs 5



# --profile=slurm
