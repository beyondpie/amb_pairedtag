# mediator
snakemake := /home/szu/miniforge3/envs/sa2/bin/snakemake
# proj_dir := /projects/ps-renlab2/szu/projects/amb_pairedtag

# tscc2
proj_dir := /tscc/projects/ps-renlab2/szu/projects/amb_pairedtag
# imac
# snakemake := /Users/szu/mambaforge/envs/sa2/bin/snakemake
# proj_dir := /Users/szu/git-recipes/mouseBrainAtlas/amb_pairedtag

work_dir := ${proj_dir}/01.clustering
pipeline_dir := ${work_dir}/src/main/pipeline
resource_dir := ${work_dir}/src/main/resource
sa2ann_dir := ${proj_dir}/data/pairedtag_ann

.PHONY: L1_tscc2 L2_tscc2
L1_tscc2: ${pipeline_dir}/sa2.scRNAseq.Snakefile ${pipeline_dir}/config.yaml
	-mkdir -p build/$@
	cp $< build/$@/Snakefile
	cp $(word 2,$^) build/$@/config.yaml
	cd build/$@ && snakemake --snakefile Snakefile -R -c 22 \
     --use-conda \
     --config \
			b2gf=${resource_dir}/pt.barcode2group.L0.20240330.csv \
			annh5ad=${sa2ann_dir}/pt.RNAseq.sa2ann4L1.h5ad \
			allenk8=${proj_dir}/meta/AIT21_k8_markers.txt \
			outdir=out/tscc2/L1 \
			npc=50 \
			k_knn=50 \
			ncpu=8 \
			clustering_level=L0 \
			use_k8=1 \
			nfeature=5000 \
			nleiden=10\
			nsilht=6\
			r_start=0.1\
			r_end=2\
			nsample_silht=50000\
			nsample_umap=20000\
			nmin_cluster=50 \
		--keep-incomplete\
		--skip-script-cleanup \
		--rerun-triggers mtime \
		--rerun-incomplete \
		--keep-going \
    --jobs 100 \
		--profile=slurm

L2_tscc2: ${pipeline_dir}/sa2.scRNAseq.Snakefile ${pipeline_dir}/config.yaml
	-mkdir -p build/$@
	cp $< build/$@/Snakefile
	cp $(word 2,$^) build/$@/config.yaml
	cd build/$@ && snakemake --snakefile Snakefile -R -c 100\
     --use-conda \
     --config \
			b2gf=${resource_dir}/pt.barcode2group.L1.20240331.csv \
			annh5ad=${sa2ann_dir}/pt.RNAseq.sa2ann4L1.h5ad \
			allenk8=${proj_dir}/meta/AIT21_k8_markers.txt \
			outdir=out/tscc2/L2 \
			npc=30 \
			k_knn=50 \
			ncpu=2 \
			clustering_level=L1 \
			use_k8=0 \
			nfeature=3000 \
			nleiden=10\
			nsilht=6\
			r_start=0.1\
			r_end=2\
			nsample_silht=50000\
			nsample_umap=20000\
			nmin_cluster=50 \
		--keep-incomplete\
		--skip-script-cleanup \
		--rerun-triggers mtime \
		--rerun-incomplete \
		--keep-going \
    --jobs 100 \
		--profile=slurm


.PHONY: L3_tscc2
L3_tscc2: ${pipeline_dir}/sa2.scRNAseq.Snakefile ${pipeline_dir}/config.yaml
	-mkdir -p build/$@
	cp $< build/$@/Snakefile
	cp $(word 2,$^) build/$@/config.yaml
	cd build/$@ && snakemake --snakefile Snakefile -R -c 200\
     --use-conda \
     --config \
			b2gf=${resource_dir}/pt.barcode2group.L2.20240331.csv \
			annh5ad=${sa2ann_dir}/pt.RNAseq.sa2ann4L1.h5ad \
			allenk8=${proj_dir}/meta/AIT21_k8_markers.txt \
			outdir=out/tscc2/L3 \
			npc=30 \
			k_knn=30 \
			ncpu=1 \
			clustering_level=L1_2 \
			use_k8=0 \
			nfeature=2000 \
			nleiden=10\
			nsilht=6\
			r_start=0.1\
			r_end=2\
			nsample_silht=50000\
			nsample_umap=20000\
			nmin_cluster=50 \
		--keep-incomplete\
		--skip-script-cleanup \
		--rerun-triggers mtime \
		--rerun-incomplete \
		--keep-going \
    --jobs 190 \
		--profile=slurm

.PHONY: L4_tscc2
L4_tscc2: ${pipeline_dir}/sa2.scRNAseq.Snakefile ${pipeline_dir}/config.yaml
	-mkdir -p build/$@
	cp $< build/$@/Snakefile
	cp $(word 2,$^) build/$@/config.yaml
	cd build/$@ && snakemake --snakefile Snakefile -R -c 200\
     --use-conda \
     --config \
			b2gf=${resource_dir}/pt.barcode2group.L3.20240401.csv \
			annh5ad=${sa2ann_dir}/pt.RNAseq.sa2ann4L1.h5ad \
			allenk8=${proj_dir}/meta/AIT21_k8_markers.txt \
			outdir=out/tscc2/L4 \
			npc=30 \
			k_knn=30 \
			ncpu=1 \
			clustering_level=L1_2_3 \
			use_k8=0 \
			nfeature=2000 \
			nleiden=10\
			nsilht=6\
			r_start=0.1\
			r_end=2\
			nsample_silht=50000\
			nsample_umap=20000\
			nmin_cluster=50 \
      filtercl=${resource_dir}/L3noL4.20240401.txt \
		--keep-incomplete\
		--skip-script-cleanup \
		--rerun-triggers mtime \
		--rerun-incomplete \
		--keep-going \
    --jobs 190 \
		--profile=slurm


.PHONY: L5_tscc2
L5_tscc2: ${pipeline_dir}/sa2.scRNAseq.Snakefile ${pipeline_dir}/config.yaml
	-mkdir -p build/$@
	cp $< build/$@/Snakefile
	cp $(word 2,$^) build/$@/config.yaml
	cd build/$@ && snakemake --snakefile Snakefile -R -c 50\
     --use-conda \
     --config \
			b2gf=${resource_dir}/pt.barcode2group.L4.20240413.csv \
			annh5ad=${sa2ann_dir}/pt.RNAseq.sa2ann4L1.h5ad \
			allenk8=${proj_dir}/meta/AIT21_k8_markers.txt \
			outdir=out/tscc2/L5 \
			npc=30 \
			k_knn=30 \
			ncpu=1 \
			clustering_level=L1_2_3_4 \
			use_k8=0 \
			nfeature=2000 \
			nleiden=10\
			nsilht=6\
			r_start=0.1\
			r_end=2\
			nsample_silht=50000\
			nsample_umap=20000\
			nmin_cluster=50 \
      filtercl=${resource_dir}/L4noL5.20240413.txt \
		--keep-incomplete\
		--skip-script-cleanup \
		--rerun-triggers mtime \
		--rerun-incomplete \
		--keep-going \
    --jobs 50 \
		--profile=slurm

test: src/main/pipeline/sa2.scRNAseq.Snakefile src/main/pipeline/config.yaml
	-mkdir -p build/$@
	cp $< build/$@/Snakefile
	cp $(word 2,$^) build/$@/config.yaml
	cd build/$@ && ${snakemake} --snakefile Snakefile -R -c 2 \
     --use-conda \
     --config \
     b2gf=${work_dir}/src/test/resource/b2g_test.csv \
     annh5ad=${work_dir}/src/test/resource/scRNAseq_sa2_test.ann.h5ad \
     allenk8=${proj_dir}/meta/AIT21_k8_markers.txt \
     outdir=out/test \
     npc=30 \
     k_knn=50 \
     ncpu=1 \
     debug=0 \
     clustering_level=L2 \
     use_k8=0 \
     nfeature=3000 \
     nleiden=2\
     nsilht=2\
     r_start=0.1\
     r_end=0.2\
     nsample_silht=100\
     nsample_umap=100\
  --keep-incomplete \
  --skip-script-cleanup \
  --rerun-triggers mtime \
  --rerun-incomplete
