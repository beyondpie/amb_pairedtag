import os._
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*
import MetaData.{PairedTagBarcodeMeta, toStringOfAllenIdLabel}
import MetaData.PairedTagBarcodeMeta.readPairedTagCellMetaFile
import MetaData.{TSCCMeta, TSCCTools}
import SZUtils.{writeStrings2File, readTable}
import MetaData.TSCCMeta.{ptclusterDNAbamd, ptscDNAbamd}

object CallPeak {
  val workd = s"${TSCCMeta.projd}/04.peakcalling"
  val rscd = s"${workd}/src/main/resource"
  val logd = s"${workd}/log/macs3_nolambda"
  val flagd = s"${workd}/flag/macs3_nolambda"
  val outdup = s"${workd}/out/macs3_nolambda"
  // val logd = s"${workd}/log/macs3_nolambda_test"
  // val flagd = s"${workd}/flag/macs3_nolambda_test"
  // val outdup = s"${workd}/out/macs3_nolambda_test"
  val macs3_encoder = "/home/szu/mambaforge/bin/macs3"

  class CallPeakElement(
    val cellCluster: String = "001_CLA_EPd_CTX_Car3_Glut",
    val isNeu: Boolean = true,
    val modality: String = "H3K27ac",
    val group: String = "Male",
    val qval: String = "0.05",
    val broad_cutoff: String = "0.1",
    val shift: String = "-100",
    val extsize: String = "200",
    val suffix: String = "srt",
    val lambda: Boolean = true,
    val savebdg: Boolean = true
  ) {
    val bamFile = ptclusterDNAbamd + "/bam/" +
      s"${cellCluster}/${modality}-${group}.${suffix}.bam"
    val outd = s"${outdup}/${cellCluster}"
    val name = s"${cellCluster}-${modality}-${group}"
    val ifbroad = if (modality == "H3K27ac") false else true
    val outFlag = s"${flagd}/${name}_macs3.done"
    val outLog = s"${logd}/${name}_macs3.log"
    def runmacs3(force: Boolean = false): Unit = {
      if (os.exists(os.Path(outFlag)) && !force) {
        println(s"${outFlag} exist, and will quit run macs3.")
      } else {
        println(s"run macs3 for ${name}.")
        TSCCTools.callPeak(
          macs3 = macs3_encoder,
          bamfnm = bamFile,
          name = name,
          shift = shift,
          extsize = extsize,
          genome = "mm",
          qvalue = qval,
          broad_cutoff = broad_cutoff,
          broad = ifbroad,
          outdir = outd,
          logfnm = outLog,
          lambda = lambda,
          savebdg = savebdg,
          cutoff_analysis = false
        )
        os.proc("touch", outFlag).call(check = false)
        println(s"run macs3 done for ${name}.")
      }
    }
    def isSuccess(): Boolean = {
      if (
        !os.exists(os.Path(outFlag)) ||
        !(os.exists(os.Path(outLog)))
      ) {
        false
      } else {
        os.read.lines
          .stream(os.Path(outLog))
          .slice(1, Int.MaxValue)
          .toList
          .last
          .contains("Done!")
      }
    }
  }

  def callPeak4SexLevel: Unit = {
    // generate CallPeakElement
    lazy val neuSexBamStat =
      MetaData.TSCCMeta.getSexBamStat(
        MetaData.TSCCMeta.NeuscSexRepBamStatFile,
        pattern = "A|B"
      )
    lazy val nnSexBamStat =
      MetaData.TSCCMeta.getSexBamStat(
        MetaData.TSCCMeta.NNspSexRepBamStatFile,
        pattern = "A|B"
      )
    lazy val sexBamStat = neuSexBamStat ++ nnSexBamStat
    lazy val sexCallPeakGroup = sexBamStat.map((k, v) =>
      (
        k,
        (
          CallPeakElement(
            cellCluster = k(0),
            isNeu = if (k(0).contains("NN")) false else true,
            modality = k(1),
            group = k(2),
            qval = "0.05",
            broad_cutoff = "0.1",
            shift = "-100",
            extsize = "200",
            lambda = false,
            // save bedgraph
            savebdg = true
          ),
          v
        )
      )
    )
    // generate the outdirs
    sexCallPeakGroup.foreach((k, v) => {
      val outd = v._1.outd
      if (!os.exists(os.Path(outd))) {
        os.makeDir(os.Path(outd))
      }
    })

    // run macs3 in parallel
    val forkJoinPool = new java.util.concurrent.ForkJoinPool(10)
    // Note: add one group
    lazy val sexCallPeakGroupPar = sexCallPeakGroup.take(1).par
    sexCallPeakGroupPar.tasksupport = new ForkJoinTaskSupport(forkJoinPool)
    sexCallPeakGroupPar.foreach { (k, v) =>
      v._1.runmacs3(force = false)
    }
    // // run macs3 in sequential
    // sexCallPeakGroup.foreach { (k, v) =>
    //   v._1.runmacs3(force = false)
    // }

    // check macs3 status
    lazy val jobStatus = sexCallPeakGroup.map((k, v) => v._1.isSuccess())
    jobStatus.map(x => if (x) 1 else 0).reduce((x, y) => x + y)
    // list of flag with failed jobs
    lazy val failedFlag =
      sexCallPeakGroup
        .filter((k, v) => !v._1.isSuccess())
        .map((k, v) => v._1.outFlag)
    failedFlag.foreach { x => os.remove(os.Path(x)) }
  }

  def callPeakAtSexRepLevel: Unit = {
    lazy val neuSexRepBamStat =
      MetaData.TSCCMeta.getSexBamStat(
        MetaData.TSCCMeta.NeuscSexRepBamStatFile,
        pattern = ""
      )
    lazy val nnSexRepBamStat =
      MetaData.TSCCMeta.getSexBamStat(
        MetaData.TSCCMeta.NNspSexRepBamStatFile,
        pattern = ""
      )
    lazy val sexRepBamStat = neuSexRepBamStat ++ nnSexRepBamStat
    lazy val sexRepCallPeakGroup = sexRepBamStat.map((k, v) =>
      (
        k,
        (
          CallPeakElement(
            cellCluster = k(0),
            isNeu = if (k(0).contains("NN")) false else true,
            modality = k(1),
            group = k(2),
            qval = "0.05",
            broad_cutoff = "0.1",
            shift = "-100",
            extsize = "200",
            lambda = false,
            savebdg = true
          ),
          v
        )
      )
    )
    sexRepCallPeakGroup.foreach((k, v) => {
      val outd = v._1.outd
      if (!os.exists(os.Path(outd))) {
        os.makeDir(os.Path(outd))
      }
    })
    // run macs3 in parallel
    val forkJoinPool = new java.util.concurrent.ForkJoinPool(5)
    lazy val sexRepCallPeakGroupPar = sexRepCallPeakGroup.par
    sexRepCallPeakGroupPar.tasksupport = new ForkJoinTaskSupport(forkJoinPool)
    sexRepCallPeakGroupPar.foreach { (k, v) =>
      v._1.runmacs3(force = false)
    }

    // check macs3 status
    lazy val jobStatus = sexRepCallPeakGroup.map((k, v) => v._1.isSuccess())
    jobStatus.map(x => if (x) 1 else 0).reduce((x, y) => x + y)
    // list of flag with failed jobs
    lazy val failedFlag =
      sexRepCallPeakGroup
        .filter((k, v) => !v._1.isSuccess())
        .map((k, v) => v._1.outFlag)
    failedFlag.foreach { x => os.remove(os.Path(x)) }
  }
  def callPeakAtSexShuffleRepLevel: Unit = {
    lazy val neuSexShufRep = readTable(
      fnm = s"${rscd}/Neu.subclass.SexShufRep.record.csv",
      sep = ",",
      head = false
    )
    lazy val nnSexShufRep = readTable(
      fnm = s"${rscd}/NN.supertype.SexShufRep.record.csv",
      sep = ",",
      head = false
    )

    lazy val sexRepBamStat = neuSexShufRep ++ nnSexShufRep
    lazy val sexRepCallPeakGroup = sexRepBamStat.map(k =>
      (
        k,
        (
          CallPeakElement(
            cellCluster = k(0),
            isNeu = if (k(0).contains("NN")) false else true,
            modality = k(1),
            group = k(2),
            qval = "0.05",
            broad_cutoff = "0.1",
            shift = "-100",
            extsize = "200",
            suffix = "shuffle.srt",
            lambda = false,
            savebdg = false
          )
        )
      )
    )
    // run macs3 in parallel
    val forkJoinPool = new java.util.concurrent.ForkJoinPool(2)
    lazy val sexRepCallPeakGroupPar = sexRepCallPeakGroup.par
    sexRepCallPeakGroupPar.tasksupport = new ForkJoinTaskSupport(forkJoinPool)
    sexRepCallPeakGroupPar.foreach { (k, v) =>
      v.runmacs3(force = false)
    }

    // // check macs3 status
    lazy val jobStatus = sexRepCallPeakGroup.map((k, v) => v.isSuccess())
    jobStatus.map(x => if (x) 1 else 0).reduce((x, y) => x + y)
    // list of flag with failed jobs
    lazy val failedFlag =
      sexRepCallPeakGroup
        .filter((k, v) => !v.isSuccess())
        .map((k, v) => v.outFlag)
    failedFlag.foreach { x => os.remove(os.Path(x)) }

  }

  def main(args: Array[String]) = {
    callPeak4SexLevel
    callPeakAtSexRepLevel
    callPeakAtSexShuffleRepLevel
  }
}
