import os._
import SZUtils.{str2path, path2str, ifelse, writeStrings2File}
import MetaData.readSubclassMeta
import bioscala.LightCoord.GenomeCoord.{GenomeCoord, GenomeCoords}
import project.Peak.loadPairedTagSubclassPeak
import bioscala.LightCoord.GenomeCoord.isOverlapIgnoreStrand
import bioscala.LightCoord.HomerMouseGenomeAnnot.{
  Annot, readGenomeElement6Cols
}
import scala.collection.parallel.*
import scala.collection.parallel.CollectionConverters.*

object CalRepressiveRegionOnTE {
  val homerPath =
    os.Path("/projects/ps-renlab2/szu/softwares/homer")

  val repeatPath = homerPath / "data" / "genomes" /
    "mm10" / "annotations" / "repeats"

  val projd =
    os.Path("/projects/ps-renlab2/szu/projects/amb_pairedtag")

  val workd = projd / "17.repressiveMarks"
  val outd  = workd / "out"

  val ptscMeta = readSubclassMeta().filter(x => x.atac)
  val ptscs    = ptscMeta.map(x => x.name.PairedTagName)
  val mods     = List("H3K27me3", "H3K9me3")

  val TEgroups =
    List("LINE", "SINE", "LTR", "Satellite", "Simple_repeat")

  val head     = (TEgroups :+ "total").mkString(",")
  val autoChrs = 1.to(19).map(i => s"chr$i").toVector

  val TEs: Map[String, GenomeCoords] =
    TEgroups
      .map(te => {
        val fnm = repeatPath / s"$te.ann.txt"
        val r   = readGenomeElement6Cols(fnm, sep = "\t")
          .map(_.g)
          .filter(x => autoChrs.contains(x.chrom))
        (te, r)
      })
      .toMap

  // val chromHMMdensed = projd / "data" / "chromHMM" / "denseBed"

  def main(args: Array[String]) = {

    mods.foreach(mod => {
      val sc2peak: Vector[(String, GenomeCoords)] =
        ptscs.map(sc => {
          (sc, loadPairedTagSubclassPeak(sc, mod))
        })

      val sc2ovlpTERatio: Vector[(String, Vector[Double])] =
        sc2peak.par
          .map((sc, peaks) => {
            val r =
              TEgroups
                .map(te => {
                  isOverlapIgnoreStrand(peaks, TEs(te))
                    .filter(x => x)
                    .length
                    .toDouble
                })
                .toVector
            val rall = r :+ peaks.length.toDouble
            (sc, rall)
          })
          .toVector

      // save results
      writeStrings2File(
          content = sc2ovlpTERatio.map((sc, r) =>
            List(sc, r.mkString(",")).mkString(",")),
          to = outd / s"pt.${mod}.peakOvlpTE.ratio.csv",
          overwrite = true,
          head = s"subclass,$head"
      )
    })

  } // end of main

} // end of object
