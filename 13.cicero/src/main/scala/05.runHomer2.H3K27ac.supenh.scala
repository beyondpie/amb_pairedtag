import os._
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*
import java.util.concurrent.ForkJoinPool
import SZUtils.{str2path, path2str}
import SZUtils.Extensions.*
import SZUtils.TaskElement
import MetaData.TSCCMeta
import Homer2.MotifFinderByHomer2

object H3K27acSupEnhBgGenome {
  val ncore: Int   = 4
  val projd        = TSCCMeta.projd
  val workd        = s"${projd}/13.cicero"
  val K27acSupEnhPeakd = s"${workd}/out/H3K27acOvlpSupEnhPeak"
  val outd         = s"${workd}/out/homer2_H3K27ac_SupEnh_bgGenome"
  val flagd        = s"${workd}/flag"
  val logd         = s"${workd}/log"

  def run() = {
    // * perform motif findings for
    val scs: List[String] =
      os.list(os.Path(K27acSupEnhPeakd))
        .map(x => x.baseName)
        .map(x => x.split("\\.")(0))
        .toList
        .distinct

    val inputfnms: Map[String, String] =
      scs
        .map(sc => {
          (sc, s"${K27acSupEnhPeakd}/${sc}.H3K27ac.ovlpsupenh.peak.bed")
        })
        .toMap

    val outds: Map[String, String] =
      scs
        .map(sc => {
          (sc, s"${outd}/${sc}.SupEnh")
        })
        .toMap

    outds.foreach((sc, d) => {
      if (!os.exists(d)) {
        os.makeDir(os.Path(d))
      }
    })

    val sct: List[MotifFinderByHomer2] = scs.map(sc => {
      val prefix = "H3K27ac.SupEnh.all.motif"
      new MotifFinderByHomer2(inputfnm = inputfnms(sc), bgfnm = None,
          flagfnm = s"${flagd}/${sc}.${prefix}.done",
          logfnm = s"${logd}/${sc}.${prefix}.log", name = sc,
          skip = true, outd = outds(sc), seqsize = 500, mask = "-mask",
          denovo = "-nomotif", hyperGenomic = " ", genome = "mm10",
          ncore = ncore)
    })
    val sctpar = sct.par
    sctpar.tasksupport = new ForkJoinTaskSupport(new ForkJoinPool(3))
    sctpar.foreach { t =>
      t.run()
    } // end of parall of motif finding
  } // end of main
}   // end of H3K27acSupEnhByGenome class

object H3K27acSupEnhCallMotif {
  def main(args: Array[String]) = {
    println("Run Homer2 default background.")
    H3K27acSupEnhBgGenome.run()
  }
}
