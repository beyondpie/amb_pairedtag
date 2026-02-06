import os._
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*
import Peak.{Peak, Peaks}
import GRange.{GenomicRange, filterByRanges}
import MetaData.TSCCMeta
import MetaData.TSCCMeta.projd
import MetaData.GenomeMeta.{mm10fafnm, mm10BlacklistBed}
import MergePeak.readReproPeak
import MergePeak.{GenomeFilter, BlacklistFilter}
import MergePeak.{BestSPMerge2, BedToolMerge2}
import SZUtils.{writeMap2File, mergeListOfMap}
import SZUtils.{readTable, writeStrings2File}
import MetaData.TSCCMeta.modality
import AllenMetaData.sc2spNN

object MergePairedTagSubclassPeak {
  val peakd = s"${projd}/data/pairedtag_peak/subclass_peak"
  val blFilter = BlacklistFilter(bedfnm = mm10BlacklistBed)
  val genomeFilter = GenomeFilter(fanm = mm10fafnm)

  def getPeakfnmsNN(sc: String, m: String, mergeTag: String): List[String] = {
    sc2spNN(sc)
      .map(sp => s"${peakd}/${sp}-${m}.${mergeTag}.peak")
      .filter(f => os.exists(os.Path(f)))
  }

  def getSexPeakfnms(sc: String, m: String): List[String] = {
    List("Female", "Male")
      .map(s => s"${peakd}/${sc}-${m}-${s}.reproPeak")
      .filter(f => os.exists(os.Path(f)))
  }

  def mergescPeak(sc: String, m: String, force: Boolean = false,mergeNNsp: Boolean = false): Unit = {
    val mergeTag = if (m == "H3K27ac") {
      "BestSPM"
    } else { "bedtoolmerge" }

    val outfnm = s"${peakd}/${sc}-${m}.${mergeTag}.peak"
    val peakfnms: List[String] = if (mergeNNsp) {
      getPeakfnmsNN(sc, m, mergeTag)
    } else { getSexPeakfnms(sc, m) }

    if (os.exists(os.Path(outfnm)) && (!force)) {
      println(s"${outfnm} exists, and skip.")
    } else if (peakfnms.length < 1) {
      println(s"No peak files found for ${sc}-${m}.")
    } else {
      val rawPeaks: Map[String, Peaks] = peakfnms
        .map(f => readReproPeak(f, head = true))
        .flatten
        .groupBy(_.r.chrom)
        .map((chrom, peaks) =>
          (chrom, blFilter.filter(genomeFilter.filter(peaks)))
        )

      val mergedPeaks: Map[String, Peaks] = if (m == "H3K27ac") {
        BestSPMerge2.merge2(rawPeaks)
      } else {
        BedToolMerge2.merge2(rawPeaks)
      }
      // output
      val vv = mergedPeaks.map((chrom, peaks) =>
        (chrom, peaks.map(p => p.mkString("\t")))
      )
      writeMap2File(
        vv,
        out = os.Path(outfnm),
        ordedKeys = GenomicRange.chromOrd,
        head = Peak.heads.mkString("\t")
      )
    } // end of else
  } // end fn

  def main(args: Array[String]) = {
    val scs =
      os.list(os.Path(peakd))
        .map(x => x.baseName)
        .filter(x => x.contains("Male") || x.contains("Female"))
        .map(x => x.split("-")(0))
        .distinct
        .toList

    scs.foreach { sc =>
      modality.foreach { m =>
        {
          println(s"Merge Peaks for ${sc} ${m}.")
          mergescPeak(sc, m, force = false, mergeNNsp = false)
        }
      } // end of merge Peaks for modality
    } // end of merge Peaks for sc

    // for NN
    sc2spNN.par.foreach((sc, sps) => {
        modality.par.foreach { m =>
          {
            println(s"Merge Peaks for NN ${sc} ${m}.")
            mergescPeak(sc, m, force = true, mergeNNsp = true)
          }
        } // end of merge modality
    }) // end of merge NN sc

  } // end of main
}
