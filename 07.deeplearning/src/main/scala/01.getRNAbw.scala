/*
 NOTE:
 After generation the RNA bigwigs with strand specificity information,
 I replace fwstrand as reverse-strand, and rvstrand as forward-strand.
 2024-10-30
 */

import os._
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*

import SZUtils.TaskElement
import MetaData.TSCCMeta.projd
import MetaData.TSCCTools.MergeBamUsingSamtools
import SZUtils.writeStrings2File
import Genome.MouseGenome.mm10Size

object GetRNABigWig4DeepLearning {
  val blv2fnm = s"${projd}/meta/mm10-blacklist.v2.bed"
  val RNAbamd = s"${projd}/data/ptRNAbam"

  val groups: List[String] =
    os.list(os.Path(RNAbamd) / "bam")
      .map(p => p.baseName.split("-")(0))
      .map(x => x.split("\\.")(0))
      .toList

  val g2bamlist: Map[String, String] =
    groups
      .map(g => {
        (g, s"${RNAbamd}/blist/${g}.mergeSexBam.bamlist")
      })
      .toList
      .toMap

  val g2bam: Map[String, String] =
    groups
      .map(g => {
        (g, s"${RNAbamd}/bam/${g}.srt.bam")
      })
      .toList
      .toMap

  class Bam2Bigwig(
    val bamfnm: String,
    val bigwigfnm: String,
    val logfnm: String,
    val flagfnm: String,
    val params: List[String],
    val skip: Boolean = true
  ) extends TaskElement {
    val cmds: List[String] =
      List("bamCoverage", "-b", bamfnm) ::: params ::: List("-o", bigwigfnm)
    def runCore(): Unit = {
      println(f"Start to generate ${bigwigfnm} from ${bamfnm}.")
      // index file firstly
      os.proc("samtools", "index", bamfnm)
        .call(check = false, stdout = os.Path(logfnm), stderr = os.Path(logfnm))
      val cmd = cmds.map(x => os.Shellable(Seq(x)))
      os.proc(cmd*)
        .call(check = false, stdout = os.Path(logfnm), stderr = os.Path(logfnm))
      println(f"Finish generating ${bigwigfnm}.")
    }
  }

  @deprecated("bamCoverage now has smooth parameter.")
  class SmoothBigwig(
    val rawbwfnm: String,
    val sbwfnm: String,
    val logfnm: String,
    val flagfnm: String,
    val smoothSize: Int = 300,
    val ncpu: Int = 1,
    val skip: Boolean = true,
    val check: Boolean = false
  ) extends TaskElement {
    val python = "/projects/ps-renlab2/earmand/conda_envs/deeptools_pip/bin/python3"
    val python_script = s"${projd}/02.track/src/main/python/smooth_bigwig.py"
    val cmds: List[String] = List(
      python,
      python_script,
      "-i",
      rawbwfnm,
      "-l",
      smoothSize.toString,
      "-o",
      sbwfnm,
      "-t",
      ncpu.toString
    )
    def runCore(): Unit = {
      println(s"Start to smooth ${rawbwfnm}.")
      val cmd = cmds.map(x => os.Shellable(Seq(x)))
      os.proc(cmd*)
        .call(check = check, stdout = os.Path(logfnm), stderr = os.Path(logfnm))
      println(f"Genereated ${sbwfnm}.")
    }
  }

  def main(args: Array[String]) = {
    // 0. Update g2bamlist
    // g2bamlist.foreach((g, f) => {
    //   val r: List[String] = List("Female", "Male")
    //     .map(x => s"${RNAbamd}/bam/${g}-${x}.bam")
    //     .filter(x => os.exists(os.Path(x)))
    //   if (r.length > 0) {
    //     println(s"Updating ${f}")
    //     writeListOfString2File(r, to = f, overwrite = true)
    //   }
    // })

    // 1. Combine Sex-related Bams
    // val combineTasks = g2bamlist.keys.map(g =>
    //   new MergeBamUsingSamtools(
    //     samtools = "samtools",
    //     flagfnm = s"${RNAbamd}/flag/${g}.CombineSexBam.done",
    //     logfnm = s"${RNAbamd}/log/${g}.CombineSexBam.log",
    //     fnmofListOfBams = g2bamlist(g),
    //     toBamfile = g2bam(g),
    //     filterMAPQ = false,
    //     skip = true
    //   )
    // )
    // combineTasks.par.foreach { t =>
    //   t.run()
    // }

    // 2. Generate bigwig without strand-specific info

    // update g2bam
    // val g2bam: Map[String, String] =
    //   groups.distinct
    //     .map(g => {
    //       (g, s"${RNAbamd}/bam/${g}.srt.bam")
    //     })
    //     .filter((g, f) => os.exists(os.Path(f)))
    //     .toList
    //     .toMap

    // val normalizeUsing = "RPKM"
    // val tobwNoStrandTasks = g2bam.map((g, bamfnm) =>
    //   new Bam2Bigwig(
    //     bamfnm,
    //     bigwigfnm = s"${RNAbamd}/bigwig/${g}.${normalizeUsing}.bw",
    //     logfnm = s"${RNAbamd}/log/${g}.${normalizeUsing}.log",
    //     flagfnm = s"${RNAbamd}/flag/${g}.${normalizeUsing}.done",
    //     params = List(
    //       "-bl",
    //       blv2fnm,
    //       "--effectiveGenomeSize",
    //       mm10Size.toString,
    //       "--normalizeUsing",
    //       normalizeUsing
    //     )
    //   )
    // )
    // tobwNoStrandTasks.par.foreach { t =>
    //   t.run()
    // }

    // // 3. Generate bigwig with strand-specific info
    // val tobwForwardStrandTasks = g2bam.map((g, bamfnm) =>
    //   new Bam2Bigwig(
    //     bamfnm,
    //     bigwigfnm = s"${RNAbamd}/bigwig/${g}.fwstrand.${normalizeUsing}.bw",
    //     logfnm = s"${RNAbamd}/log/${g}.fwstrand.${normalizeUsing}.log",
    //     flagfnm = s"${RNAbamd}/flag/${g}.fwstrand.${normalizeUsing}.done",
    //     params = List(
    //       "-bl",
    //       blv2fnm,
    //       "--effectiveGenomeSize",
    //       mm10Size.toString,
    //       "--normalizeUsing",
    //       normalizeUsing,
    //       "--filterRNAstrand",
    //       "forward"
    //     )
    //   )
    // )

    // tobwForwardStrandTasks.par.foreach { t =>
    //   t.run()
    // }

    // val tobwReverseStrandStrandTasks = g2bam.map((g, bamfnm) =>
    //   new Bam2Bigwig(
    //     bamfnm,
    //     bigwigfnm = s"${RNAbamd}/bigwig/${g}.rvstrand.${normalizeUsing}.bw",
    //     logfnm = s"${RNAbamd}/log/${g}.rvstrand.${normalizeUsing}.log",
    //     flagfnm = s"${RNAbamd}/flag/${g}.rvstrand.${normalizeUsing}.done",
    //     params = List(
    //       "-bl",
    //       blv2fnm,
    //       "--effectiveGenomeSize",
    //       mm10Size.toString,
    //       "--normalizeUsing",
    //       normalizeUsing,
    //       "--filterRNAstrand",
    //       "reverse"
    //     )
    //   )
    // )

    // tobwReverseStrandStrandTasks.par.foreach { t =>
    //   t.run()
    // }

    // 4. Generate 300-soomthed Paired-Tag DNA bigwigs (RPKM)
    val allbwfiles = os.list(os.Path(RNAbamd) / "bigwig")
    val smoothSize: Int = 300
    val SmoothBigwigTasks = allbwfiles.map(f => {
      val p = f.baseName
      new SmoothBigwig(
        rawbwfnm = f.toString,
        sbwfnm = s"${RNAbamd}/bigwig/${p}.s${smoothSize}.bw",
        logfnm = s"${RNAbamd}/log/${p}.${smoothSize}.log",
        flagfnm = s"${RNAbamd}/flag/${p}.${smoothSize}.done",
        smoothSize = smoothSize,
        ncpu = 1,
        skip = true,
        check = true
      )
    })
    SmoothBigwigTasks.par.foreach{
      t => t.run()
    }

  }
}
