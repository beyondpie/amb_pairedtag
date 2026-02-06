"""
From ROSE_main to extract the stitching enhancer procedure.

For ROSE programe, see below for further details.
1. https://bitbucket.org/young_computation/rose/src/master/
2. https://github.com/stjude/ROSE?tab=readme-ov-file
"""

import sys
import tmpPypkg.ROSE_utils as ROSE_utils
from tmpPypkg.ROSE_utils import regionStitching


# * main
inputGFFFile = sys.argv[1]
stitchWindow = int(sys.argv[2])
tssWindow = int(sys.argv[3])
stitchedGFFFile = sys.argv[4]
annotFile = sys.argv[5]

# LOADING IN THE BOUND REGION REFERENCE COLLECTION
# print("LOADING IN GFF REGIONS")
# referenceCollection = ROSE_utils.gffToLocusCollection(inputGFFFile)

# NOW STITCH REGIONS
print("STITCHING REGIONS TOGETHER")
stitchedCollection, debugOutput = regionStitching(
    inputGFFFile, stitchWindow, tssWindow, annotFile, removeTSS=True
)

# NOW MAKE A STITCHED COLLECTION GFF
print("MAKING GFF FROM STITCHED COLLECTION")
stitchedGFF = ROSE_utils.locusCollectionToGFF(stitchedCollection)

# WRITE THE GFF TO DISK
print(("WRITING STITCHED GFF TO DISK AS %s" % (stitchedGFFFile)))
ROSE_utils.unParseTable(stitchedGFF, stitchedGFFFile, "\t")

print("Finish stitchCRE.")
