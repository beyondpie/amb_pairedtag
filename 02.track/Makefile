snakemake := /home/szu/miniforge3/bin/snakemake
proj_root := /projects/ps-renlab2/szu/projects/amb_pairedtag
work_dir := ${proj_root}/02.track

.PHONY: test_split_bam split_bam split_bam_v2 split_RNA_bam resplit_RNA_bam clean

test_split_bam : src/main/pipeline/splitbam.Snakefile
	-mkdir -p build/$@
	cp $< build/$@/Snakefile
	cd build/$@ && ${snakemake} --snakefile Snakefile -R -c 3 \
    --config \
     nthread=1\
     debug=1\
     samtools=/home/szu/miniforge3/bin/samtools \
     cellmeta=${proj_root}/meta/pairedtag.meta.v1.csv \
     sublib=${work_dir}/src/test/resource/amb_pairedtag_sublibrary_v2.csv \
     outdir=${work_dir}/test \
     overwrite=0 \
    --keep-incomplete --rerun-triggers mtime \
    --skip-script-cleanup

# Need to set:
# ulimit -n 50000
# instead of 1200
# to allow multiple open files
split_bam : src/main/pipeline/splitbam.Snakefile
	-mkdir -p build/$@
	cp $< build/$@/Snakefile
	cd build/$@ && ${snakemake} --snakefile Snakefile -R -c 20 \
    --config \
     nthread=2\
     debug=0\
     samtools=/home/szu/miniforge3/bin/samtools \
     cellmeta=${proj_root}/meta/pairedtag.meta.v1.csv \
     sublib=${proj_root}/meta/amb_pairedtag_sublibrary.csv \
     overwrite=0 \
     outdir=${work_dir}/out \
    --keep-incomplete --rerun-triggers mtime \
    --skip-script-cleanup \
    --keep-going

split_bam_v2 : src/main/pipeline/splitbam.Snakefile
	-mkdir -p build/$@
	cp $< build/$@/Snakefile
	cd build/$@ && ${snakemake} --snakefile Snakefile -R -c 30 \
    --config \
     nthread=2\
     debug=0\
     samtools=/home/szu/miniforge3/bin/samtools \
     cellmeta=${work_dir}/src/main/resource/cellMeta.20240321.csv \
     sublib=${proj_root}/meta/amb_pairedtag_sublibrary_v2.csv \
     overwrite=0 \
     outdir=${work_dir}/out \
    --keep-incomplete --rerun-triggers mtime \
    --skip-script-cleanup \
    --keep-going

# Need to set:
# ulimit -n 50000
# instead of 1200
# to allow multiple open files
split_RNA_bam : src/main/pipeline/splitbam.RNA.Snakefile
	-mkdir -p build/$@
	cp $< build/$@/Snakefile
	cd build/$@ && ${snakemake} --snakefile Snakefile -R -c 30 \
    --config \
     nthread=2\
     debug=0\
     samtools=/home/szu/miniforge3/bin/samtools \
     cellmeta=${work_dir}/src/main/resource/cellMeta.20240321.csv \
     sublib=${proj_root}/meta/amb_pairedtag_sublibrary.csv \
     overwrite=1 \
     outdir=${work_dir}/out_RNA \
    --keep-incomplete --rerun-triggers mtime \
    --skip-script-cleanup \
    --keep-going

# Need to set:
# ulimit -n 50000
# instead of 1200
# to allow multiple open files
resplit_RNA_bam : src/main/pipeline/resplitbam.RNA.Snakefile
	-mkdir -p build/$@
	cp $< build/$@/Snakefile
	cd build/$@ && ${snakemake} --snakefile Snakefile -R -c 10 \
    --config \
     nthread=2\
     debug=0\
     samtools=/home/szu/miniforge3/bin/samtools \
     cellmeta=${work_dir}/src/main/resource/cellMeta.20240321.csv \
     sublib=${proj_root}/meta/amb_pairedtag_sublibrary.csv \
     overwrite=1 \
     outdir=${work_dir}/out_RNA \
    --keep-incomplete --rerun-triggers mtime \
    --skip-script-cleanup \
    --keep-going

.PHONY: test_merge_DNA test_merge_RNA
merge_script := ${proj_root}/02.track/src/main/python/03.mergebam.py
test_barcodes := ${proj_root}/02.track/src/test/resource/test_barcodes.txt
test_outdir := ${proj_root}/02.track/test
test_merge_DNA:
	python ${merge_script} -b ${test_barcodes} \
         -o ${test_outdir}/tmp_DNA.bam --omic DNA --system mediator
test_merge_RNA:
	python ${merge_script} -b ${test_barcodes} --omic RNA\
         -o ${test_outdir}/tmp_RNA.bam --system mediator


clean:
	-rm -rf test_split_bam
	-rm -rf test_out
