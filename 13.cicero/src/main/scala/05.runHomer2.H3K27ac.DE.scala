import os._
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*
import java.util.concurrent.ForkJoinPool
import SZUtils.{str2path, path2str}
import SZUtils.Extensions.*
import SZUtils.TaskElement
import MetaData.TSCCMeta
import Homer2.MotifFinderByHomer2

object H3K27acDEBgOthers {
  val ncore: Int   = 4
  val projd        = TSCCMeta.projd
  val workd        = s"${projd}/13.cicero"
  val K27acDEPeakd = s"${workd}/out/H3K27acDEPeak"
  val outd         = s"${workd}/out/homer2_H3K27ac_DE_bgOthers"
  val flagd        = s"${workd}/flag"
  val logd         = s"${workd}/log"

  def run() = {
    // * perform motif findings for
    val scs: List[String] =
      os.list(os.Path(K27acDEPeakd))
        .map(x => x.baseName)
        .map(x => x.split("\\.")(0))
        .toList
        .distinct

    val inputfnms: Map[String, String] =
      scs
        .map(sc => {
          (sc, s"${K27acDEPeakd}/${sc}.K27ac.DE.q0.05.bed")
        })
        .toMap

    val bgfnms: Map[String, String] =
      scs
        .map(sc => {
          (sc, s"${K27acDEPeakd}/${sc}.K27ac.DE.bg.bed")
        })
        .toMap

    val outds: Map[String, String] =
      scs
        .map(sc => {
          (sc, s"${outd}/${sc}.DE")
        })
        .toMap

    outds.foreach((sc, d) => {
      if (!os.exists(d)) {
        os.makeDir(os.Path(d))
      }
    })

    val sct: List[MotifFinderByHomer2] = scs.map(sc => {
      val prefix = "H3K27ac.DE.specific.motif"
      new MotifFinderByHomer2(inputfnm = inputfnms(sc),
          bgfnm = Some(bgfnms(sc)),
          flagfnm = s"${flagd}/${sc}.${prefix}.done",
          logfnm = s"${logd}/${sc}.${prefix}.log", name = sc,
          skip = false, outd = outds(sc), seqsize = 500, mask = "-mask",
          denovo = "-nomotif", hyperGenomic = " ", genome = "mm10",
          ncore = ncore)
    })
    val sctpar = sct.par
    sctpar.tasksupport = new ForkJoinTaskSupport(new ForkJoinPool(4))
    sctpar.foreach { t =>
      t.run()
    } // end of parall of motif finding

  } // end of main
}   // end of H3K27acDECiceroBgOthers

object H3K27acDEBgGenome {
  val ncore: Int   = 4
  val projd        = TSCCMeta.projd
  val workd        = s"${projd}/13.cicero"
  val K27acDEPeakd = s"${workd}/out/H3K27acDEPeak"
  val outd         = s"${workd}/out/homer2_H3K27ac_DE_bgGenome"
  val flagd        = s"${workd}/flag"
  val logd         = s"${workd}/log"

  def run() = {
    // * perform motif findings for
    val scs: List[String] =
      os.list(os.Path(K27acDEPeakd))
        .map(x => x.baseName)
        .map(x => x.split("\\.")(0))
        .toList
        .distinct

    val inputfnms: Map[String, String] =
      scs
        .map(sc => {
          (sc, s"${K27acDEPeakd}/${sc}.K27ac.DE.q0.05.bed")
        })
        .toMap

    val outds: Map[String, String] =
      scs
        .map(sc => {
          (sc, s"${outd}/${sc}.DE")
        })
        .toMap

    outds.foreach((sc, d) => {
      if (!os.exists(d)) {
        os.makeDir(os.Path(d))
      }
    })

    val sct: List[MotifFinderByHomer2] = scs.map(sc => {
      val prefix = "H3K27ac.DE.all.motif"
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
}   // end of H3K27acDECiceroBgGenome

object H3K27acDECallMotif {
  def main(args: Array[String]) = {
    // println("Run Homer2 Default background.")
    // H3K27acDEBgGenome.run()
    println("Run Homer2 Specific background.")
    H3K27acDEBgOthers.run()
  }
}
