envvars:
   "PATH"

import os
import pyprojroot

## states = [5, 10, 15, 18, 20, 23, 25]
states = [i for i in range(5, 25) if i != 18]

disk = "/tscc/projects/ps-renlab2/szu"
projd = f"{disk}/projects/amb_pairedtag"
chromHMMJar = f"{disk}/softwares/ChromHMM/ChromHMM.jar"
## bBedr = f"{projd}/06.ChromHMM/out/bBed_bypeak_b200_SPMq25"
## use raw peak binarization
bBedr = f"{projd}/06.ChromHMM/out/bBed_bypeak_b200"
java = "/tscc/nfs/home/szu/miniforge3/bin/java"
logd = f"{projd}/06.ChromHMM/log"
flagd = f"{projd}/06.ChromHMM/flag"
outd = f"{projd}/06.ChromHMM/out/"

def getOutd(wildcards):
    return f"{outd}/model_bypeak_b200_s{wildcards.s}"

rule all:
    input:
        expand("{f}/chromHMM_bpeak_s{s}.done", f = flagd,
               s = states)

rule runChromHMM:
    output:
        tag = touch(f"{flagd}/chromHMM_bpeak_s{{s}}.done")
    log:
        f"{logd}/chromHMM_bpeak_s{{s}}.log"
    threads: 20
    params:
        outd=getOutd 
    resources:
        mem = "500G",
        walltime = "24:00:00",
        queue = "condo",
        mail = "ae",
        email = "debug.pie@gmail.com",
        slurm = "qos=condo account=csd788 nodes=1"
    shell:
        """
        mkdir -p {params.outd}
        {java} -mx500000M -jar {chromHMMJar} LearnModel -p {threads} \
           {bBedr} {params.outd} {wildcards.s} mm10
        """
        
