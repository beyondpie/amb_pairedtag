projroot := /projects/ps-renlab2/szu/projects/amb_pairedtag
workdir := ${projroot}/03.integration
outdir := ${projroot}/03.integration/out

datadir := ${projroot}/data
allend := ${datadir}/allen
allen_nneu_ds := ${allend}/allen.10xv3.pt.k8.ds.raw.h5ad

snmCd := ${datadir}/snmC_snm3C
snmC_nneu_ds := ${snmCd}/snmC.scaled.gmat.k8.ds.h5ad

snATACd := ${datadir}/snATAC
snATAC_nneu_ds := ${snATACd}/snATAC.ptregion.gmat.k8.ds.raw.h5ad

ptd := ${datadir}/pairedtag_ann
pt_nneu_ds := ${ptd}/pt.RNA.ann.k8.ds.raw.h5ad

conda_dir := /home/szu/miniforge
conda := ${conda_dir}/bin/mamba
python := ${conda_dir}/envs/sa2stable/bin/python
cca_script := ${workdir}/src/main/python/04.run.cca.integration.py


nneuds_allen_pt:
	-mkdir -p ${outdir}/$@
	${python} ${cca_script} \
       --ref ${allen_nneu_ds} --refnm allen_neu \
       --query ${pt_nneu_ds} --querynm pt_neu \
       --refnorm --querynorm \
       --reflavor seurat_v3 --refhv\
       --queryflavor seurat_v3 --queryhv\
       --nhv 7500 \
       --refPCA --queryPCA --npca 50 \
       --k_anchor 5 \
       --k_filter 200 \
       --k_score 30 \
       --ncca 50 \
       --maxfea 200 \
       --maxcell 100000 \
       --k_local 0 \
       --k_weight 100 \
       --anchorpath ${outdir}/$@ \
       --inth5ad ${outdir}/$@/intgn.h5ad \
       --logfnm ${outdir}/$@/log.txt \
       --flagfnm ${outdir}/$@/flag.done \
       --run_harmony \
       --run_umap \
       --umapfnm ${outdir}/$@/umap.pdf

nneuds_allen_snmC:
	-mkdir -p ${outdir}/$@
	${python} ${cca_script} \
       --ref ${allen_nneu_ds} --refnm allen_neu \
       --query ${snmC_nneu_ds} --querynm snmC_neu \
       --refnorm \
       --reflavor seurat_v3 --refhv\
       --queryflavor seurat --queryhv\
       --nhv 7500 \
       --refPCA --queryPCA --npca 50 \
       --k_anchor 5 \
       --k_filter 200 \
       --k_score 30 \
       --ncca 50 \
       --maxfea 200 \
       --maxcell 100000 \
       --k_local 0 \
       --k_weight 100 \
       --anchorpath ${outdir}/$@ \
       --inth5ad ${outdir}/$@/intgn.h5ad \
       --logfnm ${outdir}/$@/log.txt \
       --flagfnm ${outdir}/$@/flag.done \
       --run_harmony \
       --run_umap \
       --umapfnm ${outdir}/$@/umap.pdf

nneuds_allen_snATAC:
	-mkdir -p ${outdir}/$@
	${python} ${cca_script} \
       --ref ${allen_nneu_ds} --refnm allen_neu \
       --query ${snATAC_nneu_ds} --querynm snATAC_neu \
       --refnorm --querynorm\
       --reflavor seurat_v3 --refhv\
       --queryflavor seurat_v3 --queryhv\
       --nhv 7500 \
       --refPCA --queryPCA --npca 50 \
       --k_anchor 5 \
       --k_filter 200 \
       --k_score 30 \
       --ncca 50 \
       --maxfea 200 \
       --maxcell 100000 \
       --k_local 0 \
       --k_weight 100 \
       --anchorpath ${outdir}/$@ \
       --inth5ad ${outdir}/$@/intgn.h5ad \
       --logfnm ${outdir}/$@/log.txt \
       --flagfnm ${outdir}/$@/flag.done \
       --run_harmony \
       --run_umap \
       --umapfnm ${outdir}/$@/umap.pdf

.PHONY: all nneuds_allen_pt nneuds_allen_snmC nneuds_allen_snATAC

all : nneuds_allen_pt nneuds_allen_snmC nneuds_allen_snATAC
