import os._
import MetaData.readSubclassMeta
import SZUtils.{str2path, path2str}
import bioscala.LightCoord.GenomeCoord.GenomeCoords
import bioscala.LightCoord.BedGraph._
import org.broad.igv.bbfile.BBFileReader
import bioscala.LightCoord.GenomeCoord.genomeCoordIgnoreStrandOrd
import bioscala.LightCoord.GenomeCoord.GenomeCoord
import bioscala.LightCoord.GenomeCoord.mkStringGenomeCoord
import SZUtils.{ifelse, writeStrings2File}
import smile.math.MathEx
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*
import project.Peak.loadSubclassModUnifedPeak


def loadbw(sc: String, mod: String, d: os.Path): BBFileReader = {
  val suffix = ifelse(mod == "H3K9me3", "sm1000", "sm300")
  getBigWigReader(d / s"$sc.$mod.e100.bs100.$suffix.bw")
}

@main def GetWigValuesOnUnifiedPeak(): Unit = {

  val projd              = os.Path(MetaData.TSCCMeta.projd)
  val workd              = projd / "17.repressiveMarks"
  val outd               = workd / "out" / "allGenomicRange"
  val ptscMeta           = readSubclassMeta().filter(x => x.atac)
  val ptscs              = ptscMeta.map(x => x.name.PairedTagName)
  val mods               = Vector("H3K27me3", "H3K9me3")
  val emptyValue: Double = 0.0
  val ntop: Int          = 10000

  // 1. load subclass-specific UNIFIED peaks
  val upeakd = projd / "data" / "pairedtag_peak" /
    "merge_peak" / "subclassOnUnifiedPeak"

  // Use a specific mod at first
  // val mod  = "H3K27me3"
  val mod  = "H3K9me3"
  val sc2p =
    ptscs
      .map(sc => (sc, loadSubclassModUnifedPeak(sc, mod, upeakd)))
      .toMap

  // 2. map subclass-specific bigwigs on these peaks
  val bwd   = projd / "data" / "ptDNAbam" / "bigwig"
  val sc2bw =
    ptscs
      .map(sc => (sc, loadbw(sc, mod, bwd)))
      .toMap

  val sc2vd = outd / s"subclass2unifiedPeakWigValue_$mod"

  val sc2v: Map[String, Map[GenomeCoord, Double]] =
    sc2p.par
      .map((sc, regions) => {
        val s =
          regions
            .map(g => {
              val r = getWigValue(sc2bw(sc), g, emptyValue)
              if (r.length > 1) {
                val wv =
                  r.map(y => {
                    // RPKM: width is considered
                    val width =
                      y.g.coord.endTo.min(g.coord.endTo) -
                        y.g.coord.startFrom.max(g.coord.startFrom)
                    (width.toDouble, y.s)
                  })
                val tw = wv.map(_._1).sum.toDouble
                (g, wv.map(x => x._1 * x._2).sum / tw)
              } else {
                (g, r.head.s)
              }
            })
            .toMap
        println(
            s"$sc $mod broadPeak wigValue has been calculated.")
        (sc, s)
      })
      .toVector
      .toMap

  // save sc2v
  sc2v.par.foreach { (sc, gs) =>
    {
      val x =
        gs.toVector
          .sortBy(_._1)
          .map((g, s) =>
            mkStringGenomeCoord(g, sep = "\t",
                useStrand = false) + "\t" + s.toString)
      writeStrings2File(
          content = x,
          to = sc2vd / s"$sc.$mod.broadPeak.wigValue.bed",
          overwrite = true,
          head = ""
      )

      println(s"$sc $mod broadPeak wigValue has been written.")
    }
  } // end of foreach of sc2v

} // end of main
