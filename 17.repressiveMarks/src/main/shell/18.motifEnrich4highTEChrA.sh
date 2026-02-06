#!/bin/bash

projd="/Users/szu/git-recipes/amb_pairedtag"
HOMER="/Users/szu/softwares/homer/bin"
input="${projd}/17.repressiveMarks/out/drawFig5/teCRE4homer.bed"
motifpool="${projd}/meta/jaspar2022_homer.motif"
outd="${projd}/17.repressiveMarks/out/drawFig5/homer"

$HOMER/"findMotifsGenome.pl" $input mm10 $outd -size given \
      -mask -nomotif -p 4 -homer2
