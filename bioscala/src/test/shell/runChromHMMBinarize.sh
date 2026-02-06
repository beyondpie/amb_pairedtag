#!/bin/bash

proj="/projects/ps-renlab2/szu/projects/amb_pairedtag"
rscd="${proj}/src/test/resource"
cellmarktable="${rscd}/cellmarkfiletable.tsv"
chromHMMjar="/projects/ps-renlab2/szu/softwares/ChromHMM/ChromHMM.jar"
chromSize="/projects/ps-renlab2/szu/softwares/ChromHMM/CHROMSIZES/mm10.txt"

java -mx50000M -jar $chromHMMjar BinarizeBed -b 200 -peaks $chromSize \
     $rscd $cellmarktable $rscd

