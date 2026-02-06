# projroot := /tscc/projects/ps-renlab2/szu/projects/amb_pairedtag
# conda_dir := /tscc/nfs/home/szu/miniforge3
projroot := /projects/ps-renlab2/szu/projects/amb_pairedtag
conda_dir := /home/szu/miniforge3
conda := ${conda_dir}/bin/mamba
tos5script := ${projroot}/extratools/singlecell/SnapATACAnnData2Seurat.R
tfscript := ${projroot}/extratools/singlecell/TransferLabel.R
sumtfscript := ${projroot}/03.integration/src/main/R/05.allen.transferlabel.summary.R
Rscript_bin := ${conda_dir}/envs/r/bin/Rscript
workdir := ${projroot}/03.integration
region_snakefile := ${workdir}/src/main/pipeline/transferlabel.by.region.Snakefile
snakemake := ${conda_dir}/bin/snakemake

.PHONY: sa2ann2s5 test_tf
sa2ann2s5: ${tos5script}
	${Rscript_bin} $< -f ${projroot}/01.clustering/out/scRNAseq_sa2_all.ann.h5ad \
                   -o ${projroot}/data/pairedtag_seurat/ptRNAdsL4 \
                   --conda ${conda} --condaenv sa2\
                   -m ${projroot}/meta/pt.barcode.meta.withL4.csv \
                   -d --dscol L1_2_3_4 --nds 100

.PHONY: tfneu_k8_cca_k50
tfneu_k8_cca_k50: ${tfscript}
	${Rscript_bin} $< -r ${projroot}/data/allen_seurat/allen.10xv3.pt.regions.neu.imn.k8.cl.ds100.rds \
                    -q ${projroot}/data/pairedtag_seurat/ptRNA.neu.k8.L5ds50.rds \
                    -o ${workdir}/out/$@ \
                    -f ${projroot}/meta/AIT21_k8_markers.txt \
                    --saveanchor -m cca -k 50 -p 50 -t cl --rerun --threads 3
sum_tfneu_k8_cca_k50: ${sumtfscript}
	${Rscript_bin} $< -a ${workdir}/out/tfneu_k8_cca_k50/tf.anchors.with-cca-kac50.rds \
     -q ${workdir}/out/tfneu_k8_cca_k50/query.with.tf-cca-kac50_on-cl.rds \
     --ametafnm ${projroot}/meta/AIT21_annotation_freeze_081523.tsv \
     --pmetafnm ${projroot}/meta/pt.barcode.meta.L5.csv \
     --outseufnm ${workdir}/out/tfneu_k8_cca_k50/sum_tfneu_k8_cca_k50.seurat.rds \
     --outumaprefix ${workdir}/out/tfneu_k8_cca_k50/sum_tfneu \
     --outcnssprefix ${workdir}/out/tfneu_k8_cca_k50/sum_tfneu \
     --rd cca.l2 --dmtrc euclidean --mdist 0.1 --kumap 15 \
     --lightcolor white --darkcolor red


tfneu_k8_pca_k50: ${tfscript}
	${Rscript_bin} $< -r ${projroot}/data/allen_seurat/allen.10xv3.pt.regions.neu.imn.k8.cl.ds100.rds \
                    -q ${projroot}/data/pairedtag_seurat/ptRNA.neu.k8.L5ds50.rds \
                    -o ${workdir}/out/$@ \
                    -f ${projroot}/meta/AIT21_k8_markers.txt \
                    --saveanchor -m pcaproject -k 50 -p 50 -t cl --rerun --threads 3

sum_tfneu_k8_pca_k50: ${sumtfscript}
	${Rscript_bin} $< -a ${workdir}/out/tfneu_k8_pca_k50/tf.anchors.with-pcaproject-kac50.rds \
     -q ${workdir}/out/tfneu_k8_pca_k50/query.with.tf-pcaproject-kac50_on-cl.rds \
     --ametafnm ${projroot}/meta/AIT21_annotation_freeze_081523.tsv \
     --acmetafnm ${projroot}/data/allen/allen.10xv3.cell.meta.csv \
     --pmetafnm ${projroot}/meta/pt.barcode.meta.L5.csv \
     --outseufnm ${workdir}/out/tfneu_k8_pca_k50/sum_tfneu_k8_pca_k50.seurat.rds \
     --outumaprefix ${workdir}/out/tfneu_k8_pca_k50/sum_tfneu \
     --outcnssprefix ${workdir}/out/tfneu_k8_pca_k50/sum_tfneu \
     --rd pcaproject.l2 --dmtrc euclidean --mdist 0.1 --kumap 15 \
     --lightcolor white --darkcolor red

tfneu_k8_cca_k5: ${tfscript}
	${Rscript_bin} $< -r ${projroot}/data/allen_seurat/allen.10xv3.pt.regions.neu.imn.k8.cl.ds100.rds \
                    -q ${projroot}/data/pairedtag_seurat/ptRNA.neu.k8.L5ds50.rds \
                    -o ${workdir}/out/$@ \
                    -f ${projroot}/meta/AIT21_k8_markers.txt \
                    --saveanchor -m cca -k 5 -p 50 -t cl --rerun --threads 3

sum_tfneu_k8_cca_k5: ${sumtfscript}
	${Rscript_bin} $< -a ${workdir}/out/tfneu_k8_cca_k5/tf.anchors.with-cca-kac5.rds \
     -q ${workdir}/out/tfneu_k8_cca_k5/query.with.tf-cca-kac5_on-cl.rds \
     --ametafnm ${projroot}/meta/AIT21_annotation_freeze_081523.tsv \
     --pmetafnm ${projroot}/meta/pt.barcode.meta.L5.csv \
     --outseufnm ${workdir}/out/tfneu_k8_cca_k5/sum_tfneu_k8_cca_k5.seurat.rds \
     --outumaprefix ${workdir}/out/tfneu_k8_cca_k5/sum_tfneu \
     --outcnssprefix ${workdir}/out/tfneu_k8_cca_k5/sum_tfneu \
     --rd cca.l2 --dmtrc euclidean --mdist 0.1 --kumap 15 \
     --lightcolor white --darkcolor red

tfnn_k8_cca_k5: ${tfscript}
	${Rscript_bin} $< -r ${projroot}/data/allen_seurat/allen.10xv3.pt.regions.nn.imn.k8.cl.ds1000.rds \
                    -q ${projroot}/data/pairedtag_seurat/ptRNA.nn.k8.L5ds30.rds \
                    -o ${workdir}/out/$@ \
                    -f ${projroot}/meta/AIT21_k8_markers.txt \
                    --saveanchor -m cca -k 5 -p 50 -t cl --rerun --threads 3

sum_tfnn_k8_cca_k5: ${sumtfscript}
	${Rscript_bin} $< -a ${workdir}/out/tfnn_k8_cca_k5/tf.anchors.with-cca-kac5.rds \
     -q ${workdir}/out/tfnn_k8_cca_k5/query.with.tf-cca-kac5_on-cl.rds \
     --ametafnm ${projroot}/meta/AIT21_annotation_freeze_081523.tsv \
     --pmetafnm ${projroot}/meta/pt.barcode.meta.L5.csv \
     --outseufnm ${workdir}/out/tfnn_k8_cca_k5/sum_tfnn_k8_cca_k5.seurat.rds \
     --outumaprefix ${workdir}/out/tfnn_k8_cca_k5/sum_tfnn \
     --outcnssprefix ${workdir}/out/tfnn_k8_cca_k5/sum_tfnn \
     --rd cca.l2 --dmtrc euclidean --mdist 0.1 --kumap 15 \
     --lightcolor white --darkcolor red

# * region-specific integration
tfneu_region: ${region_snakefile}
	-mkdir -p build/$@
	cp $< build/$@/Snakefile
	cd build/$@ && ${snakemake} --snakefile Snakefile -R -c 20 \
     --keep-incomplete --skip-script-cleanup \
     --rerun-triggers mtime --rerun-incomplete \
     --keep-going

test_seu_dir := ${projroot}/03.integration/src/test/resource
test_tf: ${tfscript}
	${Rscript_bin} $< -r ${test_seu_dir}/allenseu_ds.rds \
                 -q ${test_seu_dir}/ptseu_ds.rds \
                 -o ${test_seu_dir} \
                 -f ${projroot}/meta/AIT21_k8_markers.txt \
                 --saveanchor -m cca -k 20 -p 30 -t cl --rerun
