#!/usr/bin/bash

ls -la ~/softwares/homer/data/genomes/mm10/annotations/repeats |\
    awk '{print $9}' |\
    grep ann |\
    grep -v ? | grep -v Unknown | grep -v Other |\
    sed 's/\.ann\.txt$//' |\
    awk -F "|" -v OFS="," '{if (NF == 3) {print $2,$3,$1}}' |\
    sort \
     > ~/git-recipes/amb_pairedtag/meta/mm10TE_class2family2subfamily.csv
