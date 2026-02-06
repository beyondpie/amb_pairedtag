package MetaData

import SZUtils.TaskElement

object TSCCTools {
  val home: String = "/tscc/nfs/home/szu"
  val conda_root: String =
    s"${home}/miniforge3"
  val samtools: String =
    s"${conda_root}/bin/samtools"
  val macs3: String = s"${conda_root}/bin/macs3"

  def mergeBamCore(
    from: String,
    to: String,
    filterMAPQ: Boolean = true,
    MAPQ: Int = 10
  ): os.CommandResult = {
    if (filterMAPQ) {
      println(s"merge bam with filter MAPQ under ${MAPQ}.")
      os.proc(samtools, "merge", "-", "-b", from)
        .pipeTo(os.proc(samtools, "view", "-b", "-q", MAPQ, "-o", to))
        .call(check = false)
    } else {
      println("merge bam without filterMAPQ.")
      os.proc(samtools, "merge", "-f", "-b", from, "-o", to)
        .call(check = false)
    }
  }

  def callPeak(
    macs3: String,
    bamfnm: String,
    name: String,
    shift: String = "-100",
    extsize: String = "200",
    genome: String = "mm",
    qvalue: String = "0.05",
    outdir: String,
    logfnm: String,
    broad: Boolean = false,
    broad_cutoff: String = "0.1",
    lambda: Boolean = true,
    savebdg: Boolean = true,
    cutoff_analysis: Boolean = false
  ) = {
    val t0: List[String] = List(
      macs3,
      "callpeak",
      "-t",
      bamfnm,
      "-n",
      name,
      "-f",
      "BAM",
      "--outdir",
      outdir,
      "-g",
      genome,
      "-q",
      qvalue,
      "--nomodel",
      "--shift",
      shift,
      "--extsize",
      extsize,
      "--keep-dup",
      "all"
    )
    val tag2param = Map(
      (
        "broad",
        (
          List("--broad", "--broad-cutoff", broad_cutoff),
          List("--call-summits")
        )
      ),
      ("lambda", (Nil, List("--nolambda"))),
      ("bdg", (List("-B", "--SPMR"), Nil)),
      ("cutoff", (List("--cutoff-analysis"), Nil))
    )
    val tag2use = Map(
      ("broad", broad),
      ("lambda", lambda),
      ("bdg", savebdg),
      ("cutoff", cutoff_analysis)
    )
    val append2cmd = tag2use
      .map((k, v) => if v then tag2param(k)._1 else tag2param(k)._2)
      .flatten
      .toList
    val t1 = t0 ::: append2cmd
    val cmd = t1.map(x => os.Shellable(List(x)))
    os.proc(cmd*)
      .call(check = false, stdout = os.Path(logfnm), stderr = os.Path(logfnm))
  }

  class MergeBamUsingSamtools(
    val samtools: String,
    val flagfnm: String,
    val logfnm: String,
    val fnmofListOfBams: String,
    val toBamfile: String,
    val filterMAPQ: Boolean = false,
    val MAPQ: Int = 10,
    val skip: Boolean = false
  ) extends TaskElement {

    /** Merge a list of bam files using samtools.
      *
      * Bams needs to be sorted firstly before merging. The generated bam file
      * will be sorted after merge. See ref:
      * http://www.htslib.org/doc/samtools-merge.html
      */
    def runCore(): Unit = {
      if (filterMAPQ) {
        println(s"merge bam with filter MAPQ under ${MAPQ}.")
        os.proc(samtools, "merge", "-", "-b", fnmofListOfBams)
          .pipeTo(os.proc(samtools, "view", "-b", "-q", MAPQ, "-o", toBamfile))
          .call(
            check = false,
            stdout = os.Path(logfnm),
            stderr = os.Path(logfnm)
          )
      } else {
        println("merge bam without filterMAPQ.")
        os.proc(samtools, "merge", "-f", "-b", fnmofListOfBams, "-o", toBamfile)
          .call(
            check = false,
            stdout = os.Path(logfnm),
            stderr = os.Path(logfnm)
          )
      }
    } // end of runCore
  } // end of class MergeBamUsingSamtools

  class Bam2Bed(
    val bedtools: String,
    val flagfnm: String,
    val logfnm: String,
    val bamfnm: String,
    val bedfnm: String,
    val skip: Boolean = true
  ) extends TaskElement {

    def runCore(): Unit = {
      os.proc(bedtools, "bamtobed", "-i", bamfnm)
        .call(check = false, stdout = os.Path(bedfnm), stderr = os.Path(logfnm))
    }
  } // end of Bam2Bed task

  class Bed2Bam(
    val inbed: String,
    val outbam: String,
    val genomefnm: String,
    val logfnm: String,
    val flagfnm: String,
    val skip: Boolean = true,
    val check: Boolean = false
  ) extends TaskElement {
    def runCore(): Unit = {
      val reader = if (inbed.endsWith("zst")) {
        "zstdcat"
      } else {
        "cat"
      }
      os.proc(reader, inbed)
        .pipeTo(os.proc("bedToBam", "-i", "stdin", "-g", genomefnm))
        .pipeTo(os.proc("samtools", "sort", "-o", outbam))
        .call(check = check, stdout = os.Path(logfnm), stderr = os.Path(logfnm))
    } // end of fn runCore
  } // end of Bed2Bam class

  class BamCoverage2BigWig(
    val inbam: String,
    val bw: String,
    val ext: Int = 100,
    val bs: Int = 20,
    val sm: Int = 60,
    val normal: String = "RPKM",
    val ncpu: Int = 1,
    val logfnm: String,
    val flagfnm: String,
    val skip: Boolean = true,
    val check: Boolean = false,
    val skipbw: Boolean = true
  ) extends TaskElement {
    def runCore(): Unit = {
      if (os.exists(os.Path(bw)) && skipbw) {
        println(s"${bw} already exists, and skip it.")
      } else {
        if (os.exists(os.Path(s"${inbam}.bai"))) {
          println(s"Inex file exists, skip indexing.")
        } else {
          println(s"Index ${inbam}.")
          os.proc("samtools", "index", inbam)
            .call(
              check = false,
              stdout = os.Path(logfnm),
              stderr = os.Path(logfnm)
            )
        }
        println(s"Generate bigwig from ${inbam} with smooth ${sm}.")
        os.proc(
          "bamCoverage",
          "-b",
          inbam,
          "-o",
          bw,
          "-p",
          ncpu.toString,
          "-e",
          ext.toString,
          "-bs",
          bs.toString,
          "--smoothLength",
          sm.toString,
          "--normalizeUsing",
          normal
        ).call(
          check = check,
          stdout = os.Path(logfnm),
          stderr = os.Path(logfnm)
        )
        println(s"Generated ${bw}.")
      }
    } // end of fn runCore
  } // end of BamCoverage2BigWig class
} // end of object TSCCTools
