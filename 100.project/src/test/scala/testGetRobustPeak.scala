import os._
import scala.collection.immutable.HashMap
import MetaData.TSCCMeta
import Peak.{Peak, Peaks}
import SZUtils.readTable
import GRange.GenomicRange
import GRange.mouseGenomicRangeOrd
import project.Peak.GroupOfPeaks

val mod2minread = HashMap(
  ("H3K27ac", 10000),
  ("H3K27me3", 31622),
  ("H3K4me1", 31622),
  ("H3K9me3", 31622)
)

val mod2peakClass = HashMap(
  ("H3K27ac", "narrow"),
  ("H3K27me3", "broad"),
  ("H3K4me1", "broad"),
  ("H3K9me3", "broad")
)

lazy val groupStat = readTable(
  fnm = s"${TSCCMeta.projd}/04.peakcalling/src/main/resource/bam2nread.full.csv",
  sep = ",", head = false)

lazy val group2nread: Map[String, Int] = groupStat.map(
  x => {
    val suffix = x(1).split("\\.").head
    (s"${x(0)}-${suffix}", x(2).toInt)
  }
).toMap


object TestGetRobustPeak {
  val projd = TSCCMeta.projd
  val workd = s"${projd}/04.peakcalling"
  val from = s"${workd}/out/macs3_nolambda"
  val outd = s"${workd}/src/test/resource"

  val name = "001_CLA_EPd_CTX_Car3_Glut"
  val sex = "Male"
  val modality = "H3K27ac"

  val r = new GroupOfPeaks(
    from = from,
    name = name,
    sex = sex,
    modality = modality,
    neglog10qval = 2,
    outd = outd,
    flagd = outd,
    skipTask = true
  )
  // bedtools intersect -wo \
  // -a /projects/ps-renlab2/szu/projects/amb_pairedtag/04.peakcalling/out/macs3_nolambda/001_CLA_EPd_CTX_Car3_Glut/001_CLA_EPd_CTX_Car3_Glut-H3K27ac-Male_peaks.narrowPeak \
  // -b /projects/ps-renlab2/szu/projects/amb_pairedtag/04.peakcalling/out/macs3_nolambda/001_CLA_EPd_CTX_Car3_Glut/001_CLA_EPd_CTX_Car3_Glut-H3K27ac-MaleB_peaks.narrowPeak  \
  // -nonamecheck | cut -f 1-10 | sort -k1,1 -k2,2n | uniq | wc -l
  val treatfnm = r.fnms(0)
  lazy val tp = Peak.readBed(fnm = treatfnm.toString(), maxNeglog10qval = 1000, head = false)
  val repfnm = r.fnms(2)
  lazy val rp =
    Peak.readBed(fnm =repfnm.toString(), maxNeglog10qval = 1000, head = false)
      .toList
      .map(_.r)
      .groupBy(_.chrom)
      .map((k, v) => (k, v.sorted(using mouseGenomicRangeOrd).toIndexedSeq))
  tp.map(x => x.r.isInSortedRanges(rp)).map(x => if x then 1 else 0).sum

}


