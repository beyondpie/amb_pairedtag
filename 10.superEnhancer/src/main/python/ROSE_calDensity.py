"""
From ROSE_main to extract the density calculation procedure.

For ROSE programe, see below for further details.
1. https://bitbucket.org/young_computation/rose/src/master/
2. https://github.com/stjude/ROSE?tab=readme-ov-file
"""

import os
import sys
from tmpPypkg.ROSE_utils import mapCollection, gffToLocusCollection

# * main
stitchedGFFFile = sys.argv[1]
inputGFFFile = sys.argv[2]
bamFileList = [sys.argv[3]]
bamPeakGFF = sys.argv[4]
bamStitichedGFF = sys.argv[5]
outputFile1 = sys.argv[6]
stitchedGFFName = sys.argv[7]


print("BAM MAPPING COMPLETED NOW MAPPING DATA TO REGIONS")

# CALCULATE DENSITY BY REGION
stitchedCollection = gffToLocusCollection(stitchedGFFFile, window=500)
referenceCollection = gffToLocusCollection(inputGFFFile, window=500)
mappedFolder = f"{os.path.dirname(stitchedGFFFile)}/"
mapCollection(
    stitchedCollection,
    referenceCollection,
    bamFileList,
    mappedFolder,
    outputFile1,
    refName=stitchedGFFName,
)

print("Finish calDensity.")
