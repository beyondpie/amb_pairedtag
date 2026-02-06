#!/usr/bin/env python
"""
This is used for smoothing bigwigs not generating bigwigs.
"""
from deeptools import writeBedGraph_bam_and_bw
import os
from argparse import ArgumentParser
import numpy as np

def getType(fname):
    """
    function from deeptools bigwigAverage.py

    Tries to determine if a file is a wiggle file a bigWig file.
    Returns 'wiggle' if the file name ends with .wig, otherwise 'bigwig'
    """
    if fname.endswith(".wig") or fname.endswith(".wiggle"):
        return "wiggle"
    elif fname.lower().endswith(".bedgraph") or fname.endswith(".bdg"):
        return "bedgraph"
    else:
        return "bigwig"
    
def identity(tileCoverage, args):
    '''
    an identity function to map over the tileCoverage object
    '''
    # print('running function')
    # return [float(cov[0]) for cov in tileCoverage]
    return np.mean(tileCoverage)
        

def smooth_bigwig(bigwig, smooth_length, outpath=None, n_threads=1, as_bedgraph=False):
    if outpath == None:
        outpath = os.path.dirname(bigwig) + '/' + os.path.basename(bigwig)[0] + f'_smooth_{smooth_length}.bw'
    format = 'bigwig'
    if as_bedgraph:
        format = 'bedgraph'
    print(bigwig, smooth_length, outpath,getType(bigwig) )
    writeBedGraph_bam_and_bw.writeBedGraph([(bigwig, getType(bigwig))], numberOfProcessors=n_threads,
                                           outputFileName=outpath,
                                            fragmentLength=0,
                                            format=format,
                                            extendPairedEnds=False,
                                            func=identity, funcArgs=[0],
                                            smoothLength=smooth_length,
                                            # skipOverZero=False
                                            )

def parse_args():
    parser = ArgumentParser(description='A wrapper around deeptools to smooth biwgiwg tracks')
    parser.add_argument('--input', '-i',
                        help='the input file (can be bigwig, bedgraph or wiggle)',
                        required=True)
    parser.add_argument('--smooth_length', '-l',
                        help='The number of base pairs to smooth coverage. Coverage is computed as the moving '
                              'average of reads surround the given base',
                        required=True)
    parser.add_argument('--output', '-o',
                         help='a user defined output file. If not specficied files are saved in their input '
                              'directory with the final exstension _smooth_[smooth_length].bw',
                              default=None)
    parser.add_argument('--threads', '-t',
                        help='the number of threads to run the process on, default=1',
                        default=1
                        )
    parser.add_argument('--uncompressed', '-u',
                        help='whether to return outputs in uncompressed bedgraph format',
                        action='store_true')
    return parser


def main():
    args = parse_args().parse_args()
    smooth_bigwig(bigwig=args.input,
                  smooth_length=int(args.smooth_length), 
                  outpath=args.output,
                  n_threads=int(args.threads),
                  as_bedgraph=args.uncompressed
                  )
    
if __name__ == "__main__":
    main()
