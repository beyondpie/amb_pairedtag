envvars:
    "PATH"

import os
import pyprojroot

proj_dir = pyprojroot.here()
conda_dir = "/home/szu/miniforge3"
conda = f"{conda_dir}/bin/mamba"
Rscript_bin = f"{conda_dir}/envs/r/bin/Rscript"
tool_dir = f"{proj_dir}/extratools/singlecell"

work_dir = f"{proj_dir}/03.integration"
log_dir = f"{work_dir}/out/log"
out_8k_dir = f"{work_dir}/out/tfneu_8k_region_cca_k5"
out_vf_dir = f"{work_dir}/out/tfneu_vf_region_cca_k5"
for d in [log_dir, out_8k_dir, out_vf_dir]:
    os.makedirs(d, exist_ok = True)


tfscript = f"{tool_dir}/TransferLabel.R"
sumtfscript = f"{work_dir}/src/main/R/05.allen.transferlabel.summary.R"

regions = ["AMY", "CPU", "HYP", "HIP", "ERC",
           "NAC", "VTA", "PFC"]
allen_seu_dir = f"{proj_dir}/data/allen_seurat"
pt_seu_dir = f"{proj_dir}/data/pairedtag_seurat/neu_seu_region"

rule all:
    input:
        # tf_k8_flag = expand("{d}/{r}_8k_region.done", d = log_dir, r = regions),
        tf_vf_flag = expand("{d}/{r}_vf_region.done", d = log_dir, r = regions),
        # sumtf_k8_flag = expand("{d}/sumtf_{r}_8k_region.done", d = log_dir, r = regions),
        sumtf_vf_flag = expand("{d}/sumtf_{r}_vf_region.done", d = log_dir, r = regions)

rule tf_k8:
    input:
        ref = f"{allen_seu_dir}/allen.10xv3.{{r}}.neu.ds.seu.rds",
        query = f"{pt_seu_dir}/pt.neu.{{r}}.ds.seu.rds",
        k8fnm = f"{proj_dir}/meta/AIT21_k8_markers.txt"
    output:
        flag = touch(f"{log_dir}/{{r}}_8k_region.done")
    log: f"{log_dir}/{{r}}_8k_region.log"
    threads: 1
    shell:
        """
        {Rscript_bin} {tfscript} \
           -r {input.ref} -q {input.query} -f {input.k8fnm} \
           -o {out_8k_dir}/tf_{wildcards.r} \
           --saveanchor -m cca -k 5 -p 30 -t cl --rerun --threads 1
        """
rule tf_vf:
    input:
        ref = f"{allen_seu_dir}/allen.10xv3.{{r}}.neu.ds.seu.rds",
        query = f"{pt_seu_dir}/pt.neu.{{r}}.ds.seu.rds"
    output:
        flag = touch(f"{log_dir}/{{r}}_vf_region.done")
    log: f"{log_dir}/{{r}}_vf_region.log"
    threads: 1
    shell:
        """
        {Rscript_bin} {tfscript} \
           -r {input.ref} -q {input.query} --nvar 3000 --useref --usequery \
           -o {out_vf_dir}/tf_{wildcards.r} \
           --saveanchor -m cca -k 5 -p 30 -t cl --rerun --threads 1
        """
rule sumtf_k8:
    input:
        anchor = f"{out_8k_dir}/tf_{{r}}/tf.anchors.with-cca-kac5.rds",
        query = f"{out_8k_dir}/tf_{{r}}/query.with.tf-cca-kac5_on-cl.rds",
        almetafnm = f"{proj_dir}/meta/AIT21_annotation_freeze_081523.tsv",
        pmetafnm = f"{proj_dir}/meta/pairedtag.cell.meta.all.csv"
    output:
        flag = touch(f"{log_dir}/sumtf_{{r}}_8k_region.done"),
        outseufnm = f"{out_8k_dir}/sumtfneu_8k_{{r}}_cca_k5.seu.rds"
    log: f"{log_dir}/sumtf_{{r}}_8k_region.log"
    params:
        outumaprefix = f"{out_8k_dir}/sumtfneu_8k_{{r}}_cca_k5"
    shell:
        """
        {Rscript_bin} {sumtfscript} -a {input.anchor} -q {input.query} \
           --ametafnm {input.almetafnm} --pmetafnm {input.pmetafnm} \
           --outseufnm {output.outseufnm} \
           --outumaprefix {params.outumaprefix} \
           --outcnssprefix {params.outumaprefix} \
           --rd cca.l2 --dmtrc euclidean --mdist 0.1 --kumap 15 \
           --lightcolor white --darkcolor red
        """

rule sumtf_vf:
    input:
        anchor = f"{out_vf_dir}/tf_{{r}}/tf.anchors.with-cca-kac5.rds",
        query = f"{out_vf_dir}/tf_{{r}}/query.with.tf-cca-kac5_on-cl.rds",
        almetafnm = f"{proj_dir}/meta/AIT21_annotation_freeze_081523.tsv",
        pmetafnm = f"{proj_dir}/meta/pairedtag.cell.meta.all.csv"
    output:
        flag = touch(f"{log_dir}/sumtf_{{r}}_vf_region.done"),
        outseufnm = f"{out_vf_dir}/sumtfneu_vf_{{r}}_cca_k5.seu.rds"
    log: f"{log_dir}/sumtf_{{r}}_vf_region.log"
    params:
        outumaprefix = f"{out_vf_dir}/sumtfneu_vf_{{r}}_cca_k5"
    shell:
        """
        {Rscript_bin} {sumtfscript} -a {input.anchor} -q {input.query} \
           --ametafnm {input.almetafnm} --pmetafnm {input.pmetafnm} \
           --outseufnm {output.outseufnm} \
           --outumaprefix {params.outumaprefix} \
           --outcnssprefix {params.outumaprefix} \
           --rd cca.l2 --dmtrc euclidean --mdist 0.1 --kumap 15 \
           --lightcolor white --darkcolor red
        """
    
    
    
