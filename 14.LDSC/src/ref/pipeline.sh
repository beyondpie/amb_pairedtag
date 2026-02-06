### create annot files
cd ${root_dir}/subclass_regional_ldsc
conda activate ldsc

cwd=$(pwd)
plink=/tscc/projects/ps-renlab/yangli/resource/1000G_EUR_Phase3_plink
ldsc_home=~/softwares/ldsc/

function loadavg {
    while [ `cat /proc/loadavg | awk '{print int($1)}'` -gt 50 ]
    do
        sleep 120
        date
    done
}

for dir in $(ls ./); do
    echo $dir
    cd $dir
    
    mkdir -p sbatch_logs
    mkdir -p sbatch_scripts
    mkdir -p ld_score
    
    for file in peaks/*_hg19.bed; do
        prefix=$(basename $file _hg19.bed)
        echo $prefix
        mkdir -p ld_score/$prefix
        for j in {1..22}; do
            python ${ldsc_home}/make_annot.py \
                --bed-file $file \
                --bimfile ${plink}/1000G.EUR.QC.$j.bim \
                --annot-file ld_score/${prefix}/${prefix}.$j.annot.gz &
            sleep 3
            loadavg
        done
    done
    cd $cwd
done



### run ldsc
cd ${root_dir}/subclass_regional_ldsc
conda activate ldsc

ldsc_home=~/softwares/ldsc/
resources=/tscc/projects/ps-renlab/yangli/resource

cwd=$(pwd)
for dir in $(ls ./); do
    echo $dir
    cd $dir
    
    DIR=$(realpath ./)
    for file in peaks/*_hg19.bed; do
        prefix=$(basename $file _hg19.bed)
        echo $prefix
        cat >$DIR/sbatch_scripts/${prefix}.sbatch <<EOF
#! /bin/bash
#SBATCH -J ${dir}_${prefix}.ldscore
#SBATCH -A csd772
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -N 2
#SBATCH -c 8
#SBATCH -t 8:00:00
#SBATCH --mem 100G
#SBATCH -o ${DIR}/sbatch_logs/${prefix}.ldscore.out
#SBATCH -e ${DIR}/sbatch_logs/${prefix}.ldscore.err
#SBATCH --mail-user biy022@health.ucsd.edu
#SBATCH --mail-type FAIL

source ~/.bashrc
cd $DIR
conda activate ldsc

for i in {1..22}
do
    python ${ldsc_home}/ldsc.py \\
        --l2 \\
        --bfile ${resources}/1000G_EUR_Phase3_plink/1000G.EUR.QC.\$i \\
        --ld-wind-cm 1 \\
        --annot ld_score/${prefix}/${prefix}.\$i.annot.gz \\
        --thin-annot \\
        --out ld_score/${prefix}/${prefix}.\$i \\
        --print-snps ${resources}/1000G_Phase3_hapmap3/1000G_Phase3_hapmap3_print_snps.\$i.snp
done

conda deactivate
EOF
    done
    cd $cwd
done



### ldsc test
cd ${root_dir}/subclass_regional_ldsc
conda activate ldsc

sumstats=/tscc/projects/ps-renlab/yangli/resource/GWAStraits
resources=/tscc/projects/ps-renlab/yangli/resource
ldsc_home=~/softwares/ldsc
cts_name=types

cwd=$(pwd)
for dir in $(ls ./); do
    echo $dir
    cd $dir
    
    mkdir -p results
    mkdir -p sbatch_test_scripts
    mkdir -p sbatch_test_logs
    
    for file in peaks/*_hg19.bed; do
        prefix=$(basename $file _hg19.bed)
        if [[ $prefix == "union_set" ]]; then
            continue
        fi
        printf "${prefix}\tld_score/${prefix}/${prefix}.,ld_score/union_set/union_set.\n"
    done >${cts_name}.ldcts
    
    DIR=$(realpath ./)
    for file in ${sumstats}/*sumstats.gz; do
        i=$(basename $file .sumstats.gz)
        # echo $i
        cat >sbatch_test_scripts/$i.$cts_name.sbatch <<EOF
#! /bin/bash
#SBATCH -J ${i}.${cts_name}.ldsc
#SBATCH -A csd772
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -N 2
#SBATCH -c 8
#SBATCH -t 8:00:00
#SBATCH --mem 100G
#SBATCH -o $DIR/sbatch_test_logs/$i.$cts_name.ldsc.out
#SBATCH -e $DIR/sbatch_test_logs/$i.$cts_name.ldsc.err
#SBATCH --mail-user biy022@health.ucsd.edu
#SBATCH --mail-type FAIL

source ~/.bashrc
cd $DIR
conda activate ldsc

python ${ldsc_home}/ldsc.py \\
    --h2-cts $sumstats/$i.sumstats.gz \\
    --ref-ld-chr $resources/1000G_EUR_Phase3_baseline/baseline. \\
    --out results/${i}.$cts_name \\
    --ref-ld-chr-cts $cts_name.ldcts \\
    --w-ld-chr $resources/weights_hm3_no_hla/weights.
conda deactivate
EOF
    done
    
    cd $cwd
done
