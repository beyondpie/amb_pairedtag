import os._
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*
import Peak.{Peak, Peaks}
import GRange.{GenomicRange, filterByRanges, mouseGenomicRangeOrd}
import SZUtils.{readTable, writeStrings2File, mergeListOfMap}
import MetaData.TSCCMeta
import MetaData.GenomeMeta.{mm10fafnm, mm10BlacklistBed}
import MergePeak.{GenomeFilter, BlacklistFilter}
import MergePeak.readReproPeak
import MergePeak.BedToolMerge2
import SZUtils.writeMap2File

object WriteAllCombinedPeaks {
  def main(args: Array[String]) = {
    val chromOrd: List[String] = GenomicRange.chromOrd
    val mods = TSCCMeta.modality
    val m2pc = TSCCMeta.mod2peakClass
    val sexes = ("Male", "Female")
    val bl = "blv2"
    val blFilter = BlacklistFilter(bedfnm = mm10BlacklistBed)

    val genomeFilter = GenomeFilter(fanm = mm10fafnm)

    val projd = TSCCMeta.projd
    val workd = s"${projd}/04.peakcalling"
    val outModPeak: String = s"${workd}/out/mergePeak"
    val outGroupOfModPeak: String = s"${workd}/out/mergePeakPerGroup"
    List("chrX", "chrY")

    // reproducible peak directory
    val rpd = s"${workd}/out/reproPeak_nolambda"
    lazy val groups =
      os.list(os.Path(rpd))
        .map(x => x.toString.split("/").last)
        .map(x => x.replace(".reproPeak", ""))

    lazy val mod2group: Map[String, List[(String, String)]] =
      groups
        .map(x => x.split("-"))
        .groupBy(_(1))
        .map((k, v) => (k, v.toList.map(x => (x(2), x(0)))))

    lazy val group2rp: Map[String, Map[String, Peaks]] =
      groups.par
        .map(g => {
          val peaks = readReproPeak(s"${rpd}/${g}.reproPeak")
          val p_gf = genomeFilter.filter(peaks)
          val p_bl = blFilter.filter(p_gf)
          val ffp =
            // X-inactiation
            if (
              (g.contains("H3K9me3") || g.contains("h3K27me3")) & g
                .contains("Female")
            ) {
              p_bl.filter(x => x.r.chrom != "chrX")
            } else {
              p_bl
            }
          val mffp = ffp
            .groupBy(_.r.chrom)
            .map((k, v) => (k, v.sortBy(x => (x.r.startFrom, x.r.endTo))))

          (g, mffp)
        })
        .toList
        .toMap

    lazy val mod2ListOfPeaks: Map[String, List[Map[String, Peaks]]] =
      mod2group.map((k, v) => (k, v.map(x => group2rp(s"${x._2}-$k-${x._1}"))))

    mod2group.keys.par
      .map(mod => {
        (
          mod,
          mergeListOfMap(mod2ListOfPeaks(mod), chromOrd).par
            .map((k, v) => (k, v.sortBy(_.r)))
            .map((k, v) => (k, v.map(x => x.mkString("\t"))))
            .toList
            .toMap
        )
      })
      .toList
      .foreach { (mod, allPeaks) =>
        {
          println(s"save peaks to mod: ${mod} ...")
          writeMap2File(
            allPeaks,
            out = os.Path(outModPeak) /
              s"${mod}.merge.all.${bl}.repoPeak",
            ordedKeys = chromOrd,
            head = ""
          )
        }
      }
  } // end of main
}
