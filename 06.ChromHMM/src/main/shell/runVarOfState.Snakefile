envvars:
   "PATH"

import os

projd = "/tscc/projects/ps-renlab2/szu/projects/amb_pairedtag"
workd = f"{projd}/06.ChromHMM"
thejar = f"{projd}/target/scala-3.6.3/pairedtag-assembly-0.5.1.jar"
flagd = f"{workd}/flag"
logd = f"{workd}/log"

scs = [r.replace("_18_dense.bed", "")
       for r in os.listdir(f"{workd}/out/updateDenseBed")]

# test
# scs = ["001_CLA_EPd_CTX_Car3_Glut"]
rule all:
    input:
        expand("{f}/varOfState_{s}.done", f = flagd, s = scs)
        
rule varOfState:
    output:
       tag = touch(f"{flagd}/varOfState_{{s}}.done")
    log:
       f"{logd}/varOfState_{{s}}.log"
    threads: 15
    conda: "scala"
    resources:
        mem = "30G",
        walltime = "02:00:00",
        queue = "condo",
        mail = "ae",
        email = "debug.pie@gmail.com",
        slurm = "qos=condo account=csd788 nodes=1"
    shell:
        """
        java -cp {thejar} CalStateChangesArosccSubclass {wildcards.s}
        """
    
    
