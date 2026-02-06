projroot := /projects/ps-renlab2/szu/projects/amb_pairedtag
workdir := ${projroot}/03.integration
outdir := ${projroot}/03.integration/out


datadir := ${projroot}/data
allend := ${datadir}/allen
allen_neu := ${allend}/allen.10xv3.pt.regions.neu.imn.ann.k8.rawcount.h5ad
allen_nn := ${allend}/allen.10xv3.pt.regions.nn.imn.ann.k8.rawcount.h5ad

snmCd := ${datadir}/snmC_snm3C
snmC_nn := ${snmCd}/snmC.scaled.gmat.k8.nn.h5ad
snmC_neu := ${snmCd}/snmC.scaled.gmat.k8.neu.h5ad

snATACd := ${datadir}/snATAC
snATAC_nn := ${snATACd}/snATAC.ptregion.gmat.nn.h5ad
snATAC_neu := ${snATACd}/snATAC.ptregion.gmat.neu.h5ad

ptd := ${datadir}/pairedtag_ann
pt_nn := ${ptd}/pt.RNA.ann.nn.noimn.k8.h5ad
pt_neu := ${ptd}/pt.RNA.ann.neu.withimn.k8.h5ad

conda_dir := /home/szu/miniforge
conda := ${conda_dir}/bin/mamba
python := ${conda_dir}/envs/sa2stable/bin/python

cca_script := ${workdir}/src/main/python/04.run.cca.integration.py


all_allen_pt_neu:
	-mkdir -p ${outdir}/$@
	${python} ${cca_script} \
       --ref ${allen_neu} --refnm allen_neu \
       --query ${pt_neu} --querynm pt_neu \
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

all_allen_snmC_neu:
	-mkdir -p ${outdir}/$@
	${python} ${cca_script} \
       --ref ${allen_neu} --refnm allen_neu \
       --query ${snmC_neu} --querynm snmC_neu \
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

all_allen_snATAC_neu:
	-mkdir -p ${outdir}/$@
	${python} ${cca_script} \
       --ref ${allen_neu} --refnm allen_neu \
       --query ${snATAC_neu} --querynm snATAC_neu \
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

