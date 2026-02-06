import os._
import monocle.syntax.all._
import scala.collection.parallel._
import java.util.concurrent.ForkJoinPool
import scala.collection.parallel.CollectionConverters._
import SZUtils.readTable
import ZSTCompress.gz2zst
import MetaData.PairedTagSublib.{Sublib, SublibFASTQ, readSublib}


object OrgTrimmedFASTQs {
  val old = "/projects/ps-renlab"
  val d = "/projects/ps-renlab2"
  val zstd      = "/home/szu/miniforge3/bin/zstd"

  val projd     = s"${d}/szu/projects/amb_pairedtag"
  val sublibf   = s"${projd}/meta/amb_pairedtag_sublibrary.csv"
  val fastqd    = s"${d}/szu/shared_data/wmb_EP"
  val fastqDNAd = s"${fastqd}/fastqD"
  val fastqRNAd = s"${fastqd}/fastqR"

  val T         = 30
  val t         = 2

  val logd  = s"${projd}/out/log"
  val flagd = s"${projd}/out/flag"

  val old2newpath = (s"${old}/zhw063/renlab2",
      s"${d}/zhw063")

  def main(args: Array[String]) = {
    val suffixDNA = "BC_cov_trimmed.fq.gz"
    val suffixRNA = "BC_cov_trimmed_trimmed_trimmed.fq.gz"

    val subLibs = readSublib(sublibf)
      .map(x => x.focus(_.workd).modify(
        y => y.replace(old2newpath._1, old2newpath._2)))
    
    val fastqDNA: List[SublibFASTQ] = subLibs.map(x =>
      new SublibFASTQ(
        sublibid = x.sublibid,
        tag = x.idDNA,
        path = s"${x.workd}/02.trimmed/${x.idDNA}_${suffixDNA}"
      ))
    val fastqRNA: List[SublibFASTQ] = subLibs.map(x =>
      new SublibFASTQ(
        sublibid = x.sublibid,
        tag = x.idRNA,
        path = s"${x.workd}/02.trimmed/${x.idRNA}_${suffixRNA}"
      ))

    val fastqs = fastqDNA ++ fastqRNA

    val gz2zstTasks: List[gz2zst] = fastqs.map(x => {
      val prefix = s"${x.sublibid}_${x.tag}"
      new gz2zst(
        gzf = x.path,
        zstf = s"${fastqDNAd}/${prefix}_BC_cov_trimmed.fastq.zst",
        logfnm = s"${logd}/${prefix}_fastq.gz2zst.log",
        flagfnm = s"${flagd}/${prefix}_fastq.gz2zst.done",
        skip = true,
        zstd = zstd,
        T = 10
      )
    })
    val taskPar = gz2zstTasks.par
    taskPar.tasksupport = new ForkJoinTaskSupport(new ForkJoinPool(t))
    taskPar.foreach { t =>
      t.run()
    }
  } // end of main
} // end of orgTrimmedFASTQ
