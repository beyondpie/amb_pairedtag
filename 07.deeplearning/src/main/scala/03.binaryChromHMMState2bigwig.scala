import os._
import PairedTagChromHMM.chromHMMStates as states
import ChromHMM.loadDenseBed
import SimpleBigWig.BedGraphElement
import SimpleBigWig.BedGraph.{sortBedGraph, toBigWig, toBedGraph}
import GRange.GenomicRange
import SimpleBigWig.loadChromosomeData

object ChromHMMDenseBed2BigWig {
  val projd      = "/projects/ps-renlab2/szu/projects/amb_pairedtag"
  val densebedir = s"${projd}/data/chromHMM/denseBed"
  val outd  = s"${projd}/data/chromHMM/binarizedChromHMMStateAsBigWig"
  val bg2bw = "/home/szu/mambaforge/envs/deeptools/bin/bedGraphToBigWig"
  val chromSizef = s"${projd}/meta/mm10.chrom.sizes.lite"

  def main(args: Array[String]) = {
    val scs = os.list(
      os.Path(densebedir)).map(x => x.baseName.replace("_18_dense", ""))
    // test
    // val scs = List("335_BAM_NN")
    scs.foreach { sc =>
      {
        val densebed =
          loadDenseBed(f = s"${densebedir}/${sc}_18_dense.bed")._2
        val bedGraphMap =
          densebed
          .groupBy(_.name)
          .map((s, x) =>
            (s,
                x.map(y =>
                  new BedGraphElement(g = y.g.copy(), s = 1.0))))
        bedGraphMap.keys.foreach(state => {
          val ordbedGraph = sortBedGraph(bedGraphMap(state))
          val bgfnm       = s"${outd}/${sc}.${state}.bedgraph"
          val bwfnm       = s"${outd}/${sc}.${state}.bw"
          toBedGraph(x = ordbedGraph, outf = bgfnm)
          println(s"transform bedgraph to bigwig: ${state} of ${sc}.")
          toBigWig(bg2bw, chromSizef, bgfnm, bwfnm)
          // val a = loadChromosomeData(bwfnm, List("chr1"))
        })
      } // end of sc
    }   // end of scs
  }     // end of main
}       // end of object
