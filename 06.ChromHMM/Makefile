snakemake := /tscc/nfs/home/szu/miniforge3/bin/snakemake

.PHONY: runChromHMM

runChromHMM: src/main/shell/runChromHMM.Snakefile
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

mapmCHGB: src/main/shell/mapDNAmCHonGenomicBins.Snakefile
	-mkdir -p build/$@
	cp $< build/$@/Snakefile
	cd build/$@ && ${snakemake} --snakefile Snakefile -R -c 2 \
    --keep-incomplete \
    --rerun-triggers mtime \
    --rerun-incomplete \
    --skip-script-cleanup \
    --keep-going \
    --jobs 500 \
    --profile=slurm

varofstate: src/main/shell/runVarOfState.Snakefile
	-mkdir -p build/$@
	cp $< build/$@/Snakefile
	cd build/$@ && ${snakemake} --snakefile Snakefile -R -c 15 \
    --keep-incomplete \
    --rerun-triggers mtime \
    --rerun-incomplete \
    --skip-script-cleanup \
    --keep-going \
    --jobs 200 --profile=slurm
