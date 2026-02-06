import os._
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*
import GRange.GenomicRange
import MetaData.{TSCCMeta, PairedTagBarcodeMeta}
import SZUtils.TaskElement
import SZUtils.writeStrings2File

object MapHistones2ATACSummit {
  val workd = s"${TSCCMeta.projd}/data/ptHistoneCounts"
  val outd = s"${workd}/ATACPeak"
  val flagd = s"${workd}/flag"
  val logd = s"${workd}/log"
  val rawBamd = os.Path(TSCCMeta.ptclusterDNAbamd) / "bam"
  val atacSAF = (os.Path(
    CEMBATAC.wmbATACPeakfnm
  ) / os.up).toString + "/wmb.ATACpeak.annot.saf"
  class MapH2A(val sc: String, val h: String,
    val sex: String) extends TaskElement {
    val flagfnm = s"${flagd}/ATACPeak_${sc}_${h}_${sex}.done"
    val logfnm = s"${logd}/ATACPeak_${sc}_${h}_${sex}.log"
    val skip = true
    def runCore(): Unit = {
      println(s"start featureCounts on ${sc}-${h}-${sex}.")
      os.proc(
        "featureCounts",
        //"-M",
        "-O",
        "-F",
        "SAF",
        "-a",
        atacSAF,
        "-o",
        s"${outd}/${sc}_${h}_${sex}.featureCount.txt",
        s"${rawBamd}/${sc}/${h}-${sex}.srt.bam"
      ).call(check = false, stdout = os.Path(logfnm), stderr = os.Path(logfnm))
      println(s"finish featureCounts on ${sc}-${h}-${sex}.")
    }
  }

  def main(args: Array[String]) = {
    // 1. get ATAC-summit, and then extend +- 2
    val ATACPeaks: List[GenomicRange] = CEMBATAC.getATACPeaks
    val atacBiExt: List[GenomicRange] =
      ATACPeaks.map(p => p.biExtMid(extsize = 2500))
    val head = List("GeneID", "Chr", "Start", "End", "Strand")
    val atac2SAF: List[String] =
      atacBiExt
        .map(x => {
          val name = s"${x.chrom}:${x.startFrom}-${x.endTo}"
          List(name, x.chrom, x.startFrom.toString, x.endTo.toString, "*")
        })
        .map(x => x.mkString("\t"))
    writeStrings2File(
      content = head.mkString("\t") :: atac2SAF,
      to = atacSAF,
      overwrite = true
    )

    // 2. [optional] label ATAC-summit is proximal or distal
    // 3. for each subclass, map different histone reads on the region
    val rawBamd = os.Path(TSCCMeta.ptclusterDNAbamd) / "bam"
    val groups =
      os.list(rawBamd)
        .map(x => x.toString.split("/").last)
    val fullGroups =
      groups
        .map(g =>
          os.list(rawBamd / g)
            .map(d => d.toString.split("/").last)
            .filter(f =>
              f.contains("Male.srt.bam") || f.contains("Female.srt.bam")
            )
            .map(f => f.replace(".srt.bam", "").split("-").toList)
            .map(l => g :: l)
            .map(x => (x(0), x(1), x(2)))
            .toList
        )
        .toList
        .flatten
    val tasks = fullGroups.map((sc, h, sex) => new MapH2A(sc, h, sex))
    tasks.par.foreach{
      t => t.run()
    }
  }
}
