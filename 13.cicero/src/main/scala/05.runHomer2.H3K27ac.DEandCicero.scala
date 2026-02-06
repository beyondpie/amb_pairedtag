import os._
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*
import java.util.concurrent.ForkJoinPool
import SZUtils.{str2path, path2str}
import SZUtils.Extensions.*
import SZUtils.TaskElement
import MetaData.TSCCMeta
import Homer2.MotifFinderByHomer2

object H3K27acDECiceroBgOthers {
  val ncore: Int         = 4
  val projd              = TSCCMeta.projd
  val workd              = s"${projd}/13.cicero"
  val K27acDECiceroPeakd = s"${workd}/out/H3K27acDECiceroPeak"
  val K27acDEPeakd       = s"${workd}/out/H3K27acDEPeak"
  val outd               = s"${workd}/out/homer2_H3K27ac_DEandCicero_bgOthers"
  val flagd              = s"${workd}/flag"
  val logd               = s"${workd}/log"

  def main(args: Array[String]) = {
    // * perform motif findings for
    val scs: List[String] =
      os.list(os.Path(K27acDECiceroPeakd))
        .map(x => x.baseName)
        .map(x => x.split("\\.")(0))
        .toList

    val inputfnms: Map[String, String] =
      scs
        .map(sc => {
          (sc, s"${K27acDECiceroPeakd}/${sc}.K27ac.DE.cicero.bed")
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
          (sc, s"${outd}/${sc}.DE.cicero")
        })
        .toMap

    outds.foreach((sc, d) => {
      if (!os.exists(d)) {
        os.makeDir(os.Path(d))
      }
    })

    val sct: List[MotifFinderByHomer2] = scs.map(sc => {
      val prefix = "H3K27ac.cicero.DE.specific.motif"
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
} // end of H3K27acDECiceroBgOthers

object H3K27acDECiceroBgGenome {
  val ncore: Int         = 4
  val projd              = TSCCMeta.projd
  val workd              = s"${projd}/13.cicero"
  val K27acDECiceroPeakd = s"${workd}/out/H3K27acDECiceroPeak"
  val K27acDEPeakd       = s"${workd}/out/H3K27acDEPeak"
  val outd               = s"${workd}/out/homer2_H3K27ac_DEandCicero_bgGenome"
  val flagd              = s"${workd}/flag"
  val logd               = s"${workd}/log"

  def main(args: Array[String]) = {
    // * perform motif findings for
    val scs: List[String] =
      os.list(os.Path(K27acDECiceroPeakd))
        .map(x => x.baseName)
        .map(x => x.split("\\.")(0))
        .toList

    val inputfnms: Map[String, String] =
      scs
        .map(sc => {
          (sc, s"${K27acDECiceroPeakd}/${sc}.K27ac.DE.cicero.bed")
        })
        .toMap

    val outds: Map[String, String] =
      scs
        .map(sc => {
          (sc, s"${outd}/${sc}.DE.cicero")
        })
        .toMap

    outds.foreach((sc, d) => {
      if (!os.exists(d)) {
        os.makeDir(os.Path(d))
      }
    })

    val sct: List[MotifFinderByHomer2] = scs.map(sc => {
      val prefix = "H3K27ac.cicero.DE.all.motif"
      new MotifFinderByHomer2(inputfnm = inputfnms(sc),
          bgfnm = None,
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
} // end of H3K27acDECiceroBgGenome
