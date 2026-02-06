envvars:
   "PATH"

import os

projd = "/tscc/projects/ps-renlab2/szu/projects/amb_pairedtag"
workd = f"{projd}/06.ChromHMM"
thejar = f"{projd}/target/scala-3.5.1/pairedtag-assembly-0.5.1.jar"
flagd = f"{projd}/out/flag"
logd = f"{projd}/out/log"


with open(f"{workd}/out/subclass.hasDNAme.txt", 'r') as f:
    scs = [ l.strip() for l in f.readlines()]

chrs = [f"chr{i}" for i in range(1, 20)]

# test
#scs = ["001_CLA_EPd_CTX_Car3_Glut"]
#chrs = ["chr1"]

rule all:
    input:
        expand("{f}/mapDNAmCHOnGenomicBins_{s}_{c}.done",
               f = flagd, s = scs, c = chrs),
        expand("{f}/getChromHMMStatemCH_{s}_{c}.done",
               f = flagd, s = scs, c = chrs)

rule mapDNAmCH:
    output:
        tag = touch(f"{flagd}/mapDNAmCHOnGenomicBins_{{s}}_{{c}}.done")
    log:
        f"{logd}/mapDNAmCHOnGenomicBins_{{s}}_{{c}}.log"
    threads: 1
    conda: "scala"
    resources:
        mem = "120G",
        walltime = "00:45:00",
        queue = "condo",
        mail = "ae",
        email = "debug.pie@gmail.com",
        slurm = "qos=condo account=csd788 nodes=1"
    shell:
        """
        java -cp {thejar} mapDNAMethOnGenomicBins {wildcards.s} {wildcards.c}
        """
        
rule getChromHMMmCH:
    input:
        f"{flagd}/mapDNAmCHOnGenomicBins_{{s}}_{{c}}.done"
    output:
        tag = touch(f"{flagd}/getChromHMMStatemCH_{{s}}_{{c}}.done")
    log:
        f"{logd}/mapDNAmCHOnGenomicBins_{{s}}_{{c}}.log"
    threads: 1
    conda: "scala"
    resources:
        mem = "20G",
        walltime = "00:30:00",
        queue = "condo",
        mail = "ae",
        email = "debug.pie@gmail.com",
        slurm = "qos=condo account=csd788 nodes=1"
    shell:
        """
        java -cp {thejar} GetChromHMMStateDNAMeth {wildcards.s} {wildcards.c}
        """
    
