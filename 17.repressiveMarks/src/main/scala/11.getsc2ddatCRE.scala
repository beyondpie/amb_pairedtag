import os._
import bioscala.LightCoord.GenomeCoord.{GenomeCoord, GenomeCoords}
import bioscala.lightype.{DMat, Strings}
import RepressiveMarkOnTEOvlpCRE.loadSubclassesTrack
import RepressiveMarkOnTEOvlpCRE.getbwMat
import RepressiveMarkOnTEOvlpCRE.savebwMat
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*

object Getsc2DistalDARChrACREovlpTESignal {
  val projd    = os.Path("/Users/szu/git-recipes/amb_pairedtag")
  val workd    = projd / "17.repressiveMarks"
  val outd     = workd / "out" / "ovlpTE"
  val ptDNAbwd = projd / "data" / "ptDNAbw"
  val sc2TEd   = workd / "out" / "ovlpTE" / "sc2TE"

  val ptscs: Strings =
    os.list(sc2TEd)
      .map(_.baseName)
      .distinct
      .filter(x => x.length() > 1)
      .toVector
  // val hs = Vector("H3K9me3", "H3K27me3", "H3K27ac", "H3K4me1")
  val hs = Vector("H3K9me3", "H3K27ac", "H3K4me1")

  // ddat short for distal DAR ChrAs overlapped with TE
  val ddatCREfnm = outd / "distalDARChrACREs.ovlpwithTE.txt"
  val ddatCREs: GenomeCoords =
    os.read.lines
      .stream(ddatCREfnm)
      .map(x => x.strip())
      .filter(x => x.contains(":") && x.contains("-"))
      .map(x => {
        val chr_coords = x.split(":")
        val chr        = chr_coords(0)
        val coords     = chr_coords(1).split("-").map(i => i.toInt)
        (
          chrom = chr,
          coord = (startFrom = coords(0), endTo = coords(1)),
          strand = "."
        )
      }).toVector

  def main(args: Array[String]) = {
    // * get sc2ddatCRE signal matrix
    hs.par.foreach(h => {
      val sc2bw = loadSubclassesTrack(h)
      val mat = getbwMat(
        gs = ddatCREs, scs = ptscs, sc2bw = sc2bw
      )
      val outf = outd / s"sc2distalDARChrACREovlpTE.${h}"
      savebwMat(ddatCREs, ptscs, mat, outf)
      println(s"DDAT CREs on $h signals are done.")
    })

  } // end of main
}
