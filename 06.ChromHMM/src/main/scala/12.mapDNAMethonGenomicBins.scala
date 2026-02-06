import os._
//import scala.collection.parallel._
//import scala.collection.parallel.CollectionConverters.*
import scala.collection.mutable.ListBuffer
import GRange.GenomicRange
import SimpleBigWig.BedGraphElement
import AllenMetaData.Allen.{subclassNamePairedTag, subclassNameDNAMeth}
import SZUtils.{writeStrings2File, readTable}
import SimpleBigWig.loadChromosomeData
import SimpleBigWig.calculateBigWigRegionStatistic

object mapDNAMethOnGenomicBins {
  val projd     = "/tscc/projects/ps-renlab2/szu/projects/amb_pairedtag"
  val metad     = s"${projd}/meta"
  val workd     = s"${projd}/06.ChromHMM/"
  val binSize   = 200
  val chrs      = 1.to(19).map(i => s"chr${i}").toList
  val outd      = s"${workd}/out/DNAMeth_binSize${binSize}"
  val mCGd      = s"${projd}/data/snmC_snm3C/mCG_bigwig"
  val mCHd      = s"${projd}/data/snmC_snm3C/mCH_bigwig"
  val overwrite = false

  val rawscs = readTable(fnm = s"${metad}/pt.subclass2count.csv",
      sep = ",", skipEmptyField = true, head = true).map(x => x(0))

  val ptscs = rawscs.map(sc => subclassNamePairedTag(sc))

  val mscs =
    rawscs
      .map(sc => subclassNameDNAMeth(sc))
      .map(x => x.replaceAll("^\\d+_", ""))

  val ptsc2msc: Map[String, String] = ptscs.zip(mscs).toMap

  val chr2size: Map[String, Int] =
    readTable(fnm = s"${metad}/mm10.chrom.sizes.lite", sep = "\t",
        head = false)
      .map(x => (x(0), x(1).toInt))
      .toMap

  def avgByCount(x: List[(Int, Double)], binSize: Int = 200): Double = {
    x.map(_._2).sum / x.length
  }

  def avgByRangeSize(x: List[(Int, Double)],
    binSize: Int = 200): Double = {
    x.map(x => x._1 * x._2).sum / binSize
  }

  /**
   * Get BigWig Signals on Genomic Bins for one chromosome The binId
   * start from 0.
   *
   * @param bw
   * @param binSize
   * @param chromSize
   * @param ignoreLastBinIfNotEnough
   * @param emptyValue
   * @param f
   * @return
   */
  def getBigWigSignalsOnGenomicBins(bw: List[BedGraphElement],
    binSize: Int = 200, chromSize: Int,
    ignoreLastBinIfNotEnough: Boolean = true,
    f: (List[(Int, Double)], Int) => Double): List[(Int, Double)] = {
    val nbins: Int =
      if ((!ignoreLastBinIfNotEnough) && (chromSize % binSize > 0)) {
        chromSize / binSize + 1
      } else {
        // floor Int
        chromSize / binSize
      }

    val r = new ListBuffer[(Int, Int, Double)]()
    bw.map(i => {
      i.g
        .getGenomicBinIndex(binSize)
        .toList
        .filter(j => j < nbins)
        .map(j => {
          val binStart = j * binSize
          val binEnd   = binStart + binSize
          val t1       = i.g.startFrom.max(binStart)
          val t2       = i.g.endTo.min(binEnd)
          r.append((j, t2 - t1, i.s))
        })
    })
    // update r
    r.toList
      .groupBy(x => x._1)
      .map((k, v) => {
        (k, f(v.map(x => (x._2, x._3)), binSize))
      })
      .toList
      .sortBy(_._1)
  }

  def main(args: Array[String]) = {
    val scs: List[String] =
      os.list(os.Path(s"${workd}/out/model_bypeak_b200_s18"))
        .map(x => x.toString)
        .filter(x => x.contains("_18_dense.bed"))
        .map(x => x.split("/").toList.last.replace("_18_dense.bed", ""))
        .toList

    // 142 subclasses
    // val sc2CGbwf: Map[String, String] =
    //   scs
    //     .map(sc => {
    //       (sc, s"${mCGd}/${ptsc2msc(sc)}.CGN.bw")
    //     })
    //     .filter((k, v) => os.exists(os.Path(v)))
    //     .toMap

    val sc2CHbwf: Map[String, String] =
      scs
        .map(sc => {
          (sc, s"${mCHd}/${ptsc2msc(sc)}.CHN.bw")
        })
        .filter((k, v) => os.exists(os.Path(v)))
        .toMap

    // mCH
    val mC  = "mCH"
    val sc  = args(0)
    val chr = args(1)
    if (sc2CHbwf.contains(sc)) {
      val bw = loadChromosomeData(sc2CHbwf(sc), List(chr))
      val outf =
        s"${outd}/${sc}_${chr}_${mC}_b${binSize}.avgByCount.csv"
      if (os.exists(os.Path(outf)) && (!overwrite)) {
        println(s"${outf} exists, and skip it.")
      } else {
        println(s"map ${mC} for ${sc} on ${chr} ... ")
        val r = getBigWigSignalsOnGenomicBins(bw(chr),
            binSize = binSize, chromSize = chr2size(chr),
            ignoreLastBinIfNotEnough = true, f = avgByCount)
        writeStrings2File(
            content = r.map((k, v) => f"${k},${v}%.4f"),
            to = outf,
            overwrite = true,
            head = s"binId,${mC}"
        )
      }
    }

    // scs.foreach(sc => {
    //   if (sc2CHbwf.contains(sc)) {
    //     // remove par cause too much memory
    //     chrs.foreach(chr => {
    //       val bw = loadChromosomeData(sc2CHbwf(sc), List(chr))
    //       val outf =
    //         s"${outd}/${sc}_${chr}_${mC}_b${binSize}.avgByCount.csv"
    //       println(s"map ${mC} for ${sc} on ${chr} ... ")
    //       val r = getBigWigSignalsOnGenomicBins(bw(chr),
    //           binSize = binSize, chromSize = chr2size(chr),
    //           ignoreLastBinIfNotEnough = true, f = avgByCount)
    //       writeListOfString2File(
    //           content = r.map((k, v) => f"${k},${v}%.4f"),
    //           to = outf,
    //           overwrite = true,
    //           head = s"binId,${mC}"
    //       )
    //     }) // end of chrs
    //   }
    // }) // end of scs
  } // end of main
}   // end of object
