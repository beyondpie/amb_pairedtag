/** Use deeptools multipleBamSummary to calculate the similarities between
  * different subclasses.
  */
import os._
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*

object RunMultiBamSummary {
  val projd = "/projects/ps-renlab2/szu/projects/amb_pairedtag"
  val bamd = s"${projd}/data/ptDNAbam/bam"
  val histonebedir = s"${projd}/data/chromHMM/pt_subclass"
  val atacbedir = s"${projd}/data/chromHMM/snATAC_subclass"

  val bamsoftlinkd = s"${projd}/06.ChromHMM/out/subclass_bam"
  val multBamSumd = s"${projd}/06.ChromHMM/out/multiBamSum"

  val bedtools = "/home/szu/mambaforge/bin/bedtools"

  val blv2 = s"${projd}/meta/mm10-blacklist.v2.bed"
  // use softlink to put all the subclass-level bams into one directory.
  // 1. get all the subclasses
  // val atac_scs = os.list(os.Path(atacbedir)).map(x => x.baseName.split("\\.")(0))
  def get_histone_scs(hm: String): List[String] = {
    os.list(os.Path(histonebedir))
      .map(x => x.baseName)
      .filter(_.contains(hm))
      .map(x => x.split("\\.")(0))
      .toList
  }
  val H3K27ac_scs = get_histone_scs("H3K27ac")
  val H3K27me3_scs = get_histone_scs("H3K27me3")
  val H3K4me1_scs = get_histone_scs("H3K4me1")
  val H3K9me3_scs = get_histone_scs("H3K9me3")

  // 2. link different modalities' bams to corresponding directories.
  // def addSoftLink2Bam(sc: String, hm: String): Unit = {
  //   val fromfnm = s"${bamd}/${sc}/${hm}.srt.bam"
  //   val tofnm = s"${bamsoftlinkd}/${hm}/${sc}.${hm}.srt.bam"
  //   if (os.exists(os.Path(tofnm))) {
  //     println(s"${tofnm} exists, skip it.")
  //   } else {
  //     println(s"create soft link for ${fromfnm} with ${tofnm}.")
  //     os.proc("ln", "-s", fromfnm, tofnm).call(check = false)
  //   }
  // }
  // H3K27ac_scs.foreach(sc => addSoftLink2Bam(sc, "H3K27ac"))
  // H3K27me3_scs.foreach(sc => addSoftLink2Bam(sc, "H3K27me3"))
  // H3K9me3_scs.foreach(sc => addSoftLink2Bam(sc, "H3K9me3"))
  // H3K4me1_scs.foreach(sc => addSoftLink2Bam(sc, "H3K4me1"))

  // 3.multiBamSummary
  def indexBam(bam: String, check:Boolean = true): Unit = {
    os.proc("samtools", "index", bam).call(check = check)
  }
  // 164 subclasses in total, -1 (tmp remove bad bam of 045 in H3K27ac)
  val scs =
    H3K27ac_scs
      .filter(x => H3K27me3_scs.contains(x))
      .filter(x => H3K4me1_scs.contains(x))
      .filter(x => H3K9me3_scs.contains(x))
      .filter(x => ! (x == "045_OB_STR_CTX_Inh_IMN"))
  // NOTE: 045_OB_STR_CTX_Inh_IMN has been fixed.
  // Need to be updated.

  // index bam files
  // val hm = "H3K27ac"
  // for (hm <- List("H3K4me1", "H3K9me3", "H3K27me3")) {
  //   scs.par.foreach(
  //     sc => {
  //       val bam = s"${bamsoftlinkd}/${hm}/${sc}.${hm}.srt.bam"
  //       println(s"index bam: ${bam}.")
  //       if (!os.exists(os.Path(s"${bam}.bai"))) {
  //         indexBam(bam, check = true)
  //       }
  //     }
  //   )
  // }
  def main(args: Array[String]) = {
  // 4. generate multiBamSummary for four histone modifications.
    for (hm <- List("H3K27ac", "H3K4me1", "H3K9me3", "H3K27me3")) {
      println(s"Start multiBamSum: ${hm}")
      val bams: List[String] = scs.map(
        sc => s"${bamsoftlinkd}/${hm}/${sc}.${hm}.srt.bam")

      val cmdstr: List[String] = List(
        "multiBamSummary", "bins",
        "--bamfiles") ::: bams ::: List(
        "-p", 42.toString,
          "--ignoreDuplicates", "-bl", blv2,
          "--binSize", 10000.toString,
          "--smartLabels",
          "-o", s"${multBamSumd}/${hm}/${hm}_binSize_10000.npz")
      val cmd = cmdstr.map(x => os.Shellable(List(x)))
      os.proc(cmd*).call(
        check = false,
        stderr = os.Path(s"${multBamSumd}/log/${hm}_binSize_10000.log"))
      println(s"Finish multiBamSum: ${hm}")
    } // end of for hm


  } // end of main
}  // end of object
