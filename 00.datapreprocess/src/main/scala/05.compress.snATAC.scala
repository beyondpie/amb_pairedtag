import os._
import SZUtils.TaskElement
import ZSTCompress.gz2zst
import scala.collection.parallel._
import java.util.concurrent.ForkJoinPool
import scala.collection.parallel.CollectionConverters._

object CompressATACzst {
  // * tools
  // val zstd = "/home/szu/mambaforge/bin/zstd"
  val zstd = "/home/szu/miniforge3/bin/zstd"
  // multithreads for zstd
  val T = 10
  // multithreads for tasks
  val t = 2

  // * meta
  val atacd = "/projects/ps-renlab2/szu/projects/CEMBA_wmb_snATAC"
  val fromd = s"${atacd}/data/fastq"
  val projd = "/projects/ps-renlab2/szu/projects/amb_pairedtag"
  val logd  = s"${projd}/out/log"
  val flagd = s"${projd}/out/flag"
  val outd  = s"/projects/ps-renlab2/szu/shared_data/wmb_snATAC/fastq"

  def main(args: Array[String]) = {
    // * all the fastq files
    val fastqR1 =
      os.list(os.Path(fromd) / "R1")
        .map(x => x.toString)
        .map(x => (x.split("/").last.replace("\\.gz$", ""), x))

    val fastqR2 =
      os.list(os.Path(fromd) / "R2")
        .map(x => x.toString)
        .map(x => (x.split("/").last.replace("\\.gz$", ""), x))

    val fastqs = fastqR1 ++ fastqR2

    // * create command line task for each of them
    val fastqTasks = fastqs.map((k, v) =>
      new gz2zst(
          gzf = v,
          zstf = s"${outd}/${k}.zst",
          logfnm = s"${logd}/${k}.gz2zst.log",
          flagfnm = s"${flagd}/${k}.gz2zst.done",
          skip = true,
          zstd = zstd,
          T = T
      ))

    val parFastqTasks = fastqTasks.par
    parFastqTasks.tasksupport = new ForkJoinTaskSupport(
        new ForkJoinPool(t)
    )
    parFastqTasks.foreach(x => x.run())
    // test
    // val a = fastqTasks(0)
  }
}
