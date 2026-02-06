projroot := /projects/ps-renlab2/szu/projects/amb_pairedtag
workdir := ${projroot}/03.integration
outdir := ${projroot}/03.integration/out


datadir := ${projroot}/data
allend := ${datadir}/allen
allen_nn := ${allend}/allen.10xv3.pt.nn.ds.h5ad
allen_neu := ${allend}/allen.10xv3.pt.neu.ds.h5ad

snmCd := ${datadir}/snmC_snm3C
snmC_nn := ${snmCd}/snmC.scaled.gmat.nn.ds.h5ad
snmC_neu := ${snmCd}/snmC.scaled.gmat.neu.ds.h5ad

snATACd := ${datadir}/snATAC
snATAC_nn := ${snATACd}/snATAC.gmat.nn.ds.h5ad
snATAC_neu := ${snATACd}/snATAC.gmat.neu.ds.h5ad

ptd := ${datadir}/pairedtag_ann
pt_nn := ${ptd}/pt.RNA.ann.nn.noimn.k8.ds.h5ad
pt_neu := ${ptd}/pt.RNA.ann.neu.withimn.k8.ds.h5ad

conda_dir := /home/szu/miniforge
conda := ${conda_dir}/bin/mamba
python := ${conda_dir}/envs/sa2stable/bin/python

cca_script := ${workdir}/src/main/python/04.run.cca.integration.py

.PHONY = int_allen_snmC_nn
int_allen_snmC_nn:
	-mkdir -p ${outdir}/$@
	${python} ${cca_script} \
       --ref ${allen_nn} --refnm allen_nn \
       --query ${snmC_nn} --querynm snmC_nn \
       --refnorm \
       --reflavor seurat_v3 --refhv\
       --queryflavor seurat --queryhv\
       --nhv 3000 \
       --refPCA --queryPCA --npca 50 \
       --k_anchor 5 \
       --k_filter 200 \
       --k_score 30 \
       --ncca 50 \
       --maxfea 200 \
       --maxcell 200000 \
       --k_local 0 \
       --k_weight 100 \
       --anchorpath ${outdir}/$@ \
       --inth5ad ${outdir}/$@/intgn.h5ad \
       --logfnm ${outdir}/$@/log.txt \
       --flagfnm ${outdir}/$@/flag.done \
       --run_harmony \
       --run_umap \
       --umapfnm ${outdir}/$@/umap.pdf


int_allen_snmC_neu:
	-mkdir -p ${outdir}/$@
	${python} ${cca_script} \
       --ref ${allen_neu} --refnm allen_neu \
       --query ${snmC_neu} --querynm snmC_neu \
       --refnorm \
       --reflavor seurat_v3 --refhv\
       --queryflavor seurat --queryhv\
       --nhv 4000 \
       --refPCA --queryPCA --npca 50 \
       --k_anchor 5 \
       --k_filter 200 \
       --k_score 30 \
       --ncca 50 \
       --maxfea 200 \
       --maxcell 200000 \
       --k_local 0 \
       --k_weight 100 \
       --anchorpath ${outdir}/$@ \
       --inth5ad ${outdir}/$@/intgn.h5ad \
       --logfnm ${outdir}/$@/log.txt \
       --flagfnm ${outdir}/$@/flag.done \
       --run_harmony \
       --run_umap \
       --umapfnm ${outdir}/$@/umap.pdf

int_allen_snATAC_nn:
	-mkdir -p ${outdir}/$@
	${python} ${cca_script} \
       --ref ${allen_nn} --refnm allen_nn \
       --query ${snATAC_nn} --querynm snATAC_nn \
       --refnorm --querynorm \
       --reflavor seurat_v3 --refhv\
       --queryflavor seurat_v3 --queryhv\
       --nhv 4000 \
       --refPCA --queryPCA --npca 50 \
       --k_anchor 5 \
       --k_filter 200 \
       --k_score 30 \
       --ncca 50 \
       --maxfea 200 \
       --maxcell 200000 \
       --k_local 0 \
       --k_weight 100 \
       --anchorpath ${outdir}/$@ \
       --inth5ad ${outdir}/$@/intgn.h5ad \
       --logfnm ${outdir}/$@/log.txt \
       --flagfnm ${outdir}/$@/flag.done \
       --run_harmony \
       --run_umap \
       --umapfnm ${outdir}/$@/umap.pdf


int_allen_snATAC_neu:
	-mkdir -p ${outdir}/$@
	${python} ${cca_script} \
       --ref ${allen_neu} --refnm allen_neu \
       --query ${snATAC_neu} --querynm snATAC_neu \
       --refnorm --querynorm\
       --reflavor seurat_v3 --refhv\
       --queryflavor seurat_v3 --queryhv\
       --nhv 4000 \
       --refPCA --queryPCA --npca 50 \
       --k_anchor 5 \
       --k_filter 200 \
       --k_score 30 \
       --ncca 50 \
       --maxfea 200 \
       --maxcell 200000 \
       --k_local 0 \
       --k_weight 100 \
       --anchorpath ${outdir}/$@ \
       --inth5ad ${outdir}/$@/intgn.h5ad \
       --logfnm ${outdir}/$@/log.txt \
       --flagfnm ${outdir}/$@/flag.done \
       --run_harmony \
       --run_umap \
       --umapfnm ${outdir}/$@/umap.pdf

int_allen_pt_nn:
	-mkdir -p ${outdir}/$@
	${python} ${cca_script} \
       --ref ${allen_nn} --refnm allen_nn \
       --query ${pt_nn} --querynm pt_nn \
       --refnorm --querynorm \
       --reflavor seurat_v3 --refhv\
       --queryflavor seurat_v3 --queryhv\
       --nhv 4000 \
       --refPCA --queryPCA --npca 50 \
       --k_anchor 5 \
       --k_filter 200 \
       --k_score 30 \
       --ncca 50 \
       --maxfea 200 \
       --maxcell 200000 \
       --k_local 0 \
       --k_weight 100 \
       --anchorpath ${outdir}/$@ \
       --inth5ad ${outdir}/$@/intgn.h5ad \
       --logfnm ${outdir}/$@/log.txt \
       --flagfnm ${outdir}/$@/flag.done \
       --run_harmony \
       --run_umap \
       --umapfnm ${outdir}/$@/umap.pdf

int_allen_pt_neu:
	-mkdir -p ${outdir}/$@
	${python} ${cca_script} \
       --ref ${allen_neu} --refnm allen_neu \
       --query ${pt_neu} --querynm pt_neu \
       --refnorm --querynorm \
       --reflavor seurat_v3 --refhv\
       --queryflavor seurat_v3 --queryhv\
       --nhv 4000 \
       --refPCA --queryPCA --npca 50 \
       --k_anchor 5 \
       --k_filter 200 \
       --k_score 30 \
       --ncca 50 \
       --maxfea 200 \
       --maxcell 200000 \
       --k_local 0 \
       --k_weight 100 \
       --anchorpath ${outdir}/$@ \
       --inth5ad ${outdir}/$@/intgn.h5ad \
       --logfnm ${outdir}/$@/log.txt \
       --flagfnm ${outdir}/$@/flag.done \
       --run_harmony \
       --run_umap \
       --umapfnm ${outdir}/$@/umap.pdf
