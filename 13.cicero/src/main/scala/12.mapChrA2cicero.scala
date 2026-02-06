/**
 * Load previous cicero results, and check the overlaps with our Chr-A
 * enhancers. Ref codes: \- load subclass-specific ppdc
 * [CEMBA2]/src/main/R/09.sa2.subclass.specific.ppdc.R
 */

import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters._
import io._
import Bed.{BedPE10Element, BedElement4}
import Bed.BedPE10Element.writeBedPE
import MetaData.TSCCMeta.projd
import MetaData.readSubclassMetaFromChromHMM
import CEMBATAC.{ATACPeakd, loadAnnotscATAC}
import CEMBATAC.{globalppdc, loadscpdc}
import GRange.{mouseGenomicRangeOrd, mouseChrOrd, GenomicRange}
import PairedTagBigWig.{
  loadptRNAbw, BigWigScoreTask, ptDNAbwd, getptmodSuffix
}
import SimpleBigWig.loadChromosomeData
import SimpleBigWig.mapBigWigOnRegion
import SimpleBigWig.mapBigWigOnRegionOneChr
import SZUtils.writeStrings2File
import SZUtils.TaskElement
import SimpleBigWig.loadBigWigOneChr
import SimpleBigWig.BedGraphElement
import PairedTagBigWig.loadptRNAbwOneChr
import PairedTagBigWig.loadptDNAbw
import PairedTagBigWig.loadptDNAbwOneChr
import PairedTagBigWig.ptRNAbwd
import SZUtils.readTable
import MetaData.readSubclassMeta

@main def mapChrA2Cicero() = {
  val outd = s"$projd/13.cicero/out/ATACppec"
  if (!os.exists(os.Path(outd))) {
    os.makeDir(os.Path(outd))
  }
  val flagd             = s"$projd/13.cicero/flag"
  val ptscMeta          = readSubclassMeta()
  val scs = ptscMeta.map(x => x.name.PairedTagName)

  val ptsc2atacnm =
    ptscMeta
      .map(x =>
        (x.name.PairedTagName,
            x.name.ATACName.replaceAll("^\\d+_", "")))
      .toMap

  // // 1. load chrA-enhancer results per subclass
  // val sc2ATAC: Map[String, List[BedElement4]] =
  //   scs.map(sc => (sc, loadAnnotscATAC(sc))).toMap

  // val sc2enhs =
  //   scs
  //     .map(sc => (sc, sc2ATAC(sc).filter(x => x.name == "Chr-A")))
  //     .toMap
  // val allEnhs = sc2enhs.values.toList.flatten.distinct

  // // 2. per subclass, get ppdcs related with chrA-enhancer
  // val ppdc      = globalppdc
  // val sc2pdc    = scs.map(sc => (sc, loadscpdc(ptsc2atacnm(sc)))).toMap
  // val ppdcnames = ppdc.map(x => x.name).toSet
  // val sc2ppdc =
  //   sc2pdc
  //     .map((sc, pdcs) =>
  //       (sc, pdcs.filter(pd => ppdcnames.contains(pd.name))))
  //     .filter((sc, pdcs) => pdcs.length > 0)

  // // 3. per subclass, get ppdcs with positive correlated genes
  // //    related with enhancers
  // scs.par.foreach(sc => {
  //     val enhs = sc2enhs(sc).map(e => e.x).toSet
  //     val ppec =
  //       sc2ppdc(sc)
  //         .filter(pdc => enhs.contains(pdc.y))
  //         .sortBy(x => x.y)(mouseGenomicRangeOrd)
  //     writeBedPE(ppec, fnm = s"$outd/${sc}.ppec.bedpe")
  //   })

  // load sc2ppec
  val sc2ppec =
    scs
      .map(sc => {
        val ppec = BedPE10Element.readBedPE(
            fnm = s"$outd/${sc}.ppec.bedpe", head = false)
        (sc, ppec)
      })
      .toMap
  val allPosEnhs =
    sc2ppec.values.flatten
      .map(pdc => pdc.y)
      .toVector
      .sorted(using mouseGenomicRangeOrd)
      .distinct

  val ext: Int = 500
  val allPosEnhExts = allPosEnhs.map(x =>
    new GenomicRange(x.chrom, x.startFrom - 500, x.endTo + 500))
  // val EnhExt2Enh = allPosEnhExts.zip(allPosEnhs).toMap

  val chr2posEnh =
    allPosEnhs
      .groupBy(_.chrom)
      .map((chr, elms) =>
        (chr, elms.sorted(using mouseGenomicRangeOrd)))
      .toMap
  val chr2posExtEnh =
    allPosEnhExts
      .groupBy(_.chrom)
      .map((chr, elms) =>
        (chr, elms.sorted(using mouseGenomicRangeOrd)))
      .toMap

  // 4. get eRNA, H3K27ac, H3K4me1 signals on allPosEnhs
  // take > 13 hours for 151 subclasses and eRNA, H3K27ac, H3K4me1
  val chrs = 1.to(19).map(i => s"chr$i").toList

  val mods = List("H3K27ac", "H3K4me1", "H3K27me3")
  val hTasks = scs.flatMap(sc => {
    mods.map { mod =>
      {
        val prefix = s"${sc}.${mod}.${getptmodSuffix(mod)}"
        new BigWigScoreTask(
            bwf = s"$ptDNAbwd/${prefix}.bw",
            chr2CREs = chr2posExtEnh,
            outd = outd,
            flagd = flagd,
            prefix = prefix
        )
      }
    }
  })
  hTasks.par.foreach(t => {
    t.run()
    println(s"${t.prefix} is done.")
  })

  // [optional] check if they are significantly within
  // TAD boundaries.
  // [optional] check how many of them are located in
  // TE regions.

  // 5. merge all the singles as matrix
  def readSignals(fnm: String, leftShift: Int = 0,
    rightShift: Int = 0): Map[GenomicRange, Double] = {
    readTable(fnm = fnm, sep = ",", head = true)
      .map(x =>
        (GenomicRange(x(0), x(1).toInt + leftShift,
                x(2).toInt + rightShift), x(3).toDouble))
      .toMap
  }

  def mergeBigWigSignals_(
    sc: String): Vector[(GenomicRange, Vector[Double])] = {
    val eRNAMap = readSignals(s"$outd/${sc}.eRNA.csv")
    val tmpsuf  = "e100.bs100.sm300"
    val H3K4me1Map = readSignals(s"$outd/${sc}.H3K4me1.${tmpsuf}.csv",
        leftShift = 500, rightShift = -500)
    val H3K27acMap = readSignals(s"$outd/${sc}.H3K27ac.${tmpsuf}.csv",
        leftShift = 500, rightShift = -500)
    val H3K27me3Map = readSignals(s"$outd/${sc}.H3K27me3.${tmpsuf}.csv",
        leftShift = 500, rightShift = -500)
    eRNAMap
      .map((g, s) => {
        (g, Vector(s, H3K4me1Map(g), H3K27acMap(g), H3K27me3Map(g)))
      })
      .toVector
      .sortBy(_._1)(using mouseGenomicRangeOrd)
  }

  /**
   * Merge Signals from different subclasses for one modality.
   *
   * Based on the task above,
   *   the regions in each subclass are aligned and well sorted.
   *
   * @param fnms
   * @return
   * The genomic ranges are exactly the ones in the first file.
   * The vectors will follow the order of fnms
   */
  def mergeBigWgSignals(fnms: Vector[String], leftShift: Int = 0,
    rightShift: Int =
      0): (Vector[GenomicRange], Vector[Vector[Double]]) = {
    val grs = readTable(fnm = fnms(0), sep = ",", head = true)
      .map(x =>
        GenomicRange(x(0), x(1).toInt + leftShift,
            x(2).toInt + rightShift))
      .toVector
    val signal = fnms
      .map(f => {
        readTable(fnm = f, sep = ",", head = true)
          .map(x => x.last.toDouble)
          .toVector
      })
      .transpose
    (grs, signal)
  }

  def outSignal(g: Vector[GenomicRange], sc: Vector[String],
    m: Vector[Vector[Double]], outd: os.Path): Unit = {
    if (!os.exists(outd)) {
      os.makeDir(outd)
    }
    writeStrings2File(
        content = g.map(x => x.mkString(sep = "\t")).toList,
        to = (outd / "CRE.bed").toString,
        overwrite = true,
        head = ""
    )
    writeStrings2File(
        content = sc.toList,
        to = (outd / "subclass.csv").toString,
        overwrite = true,
        head = ""
    )
    writeStrings2File(
        content = m.map(x => x.mkString(sep = ",")).toList,
        to = (outd / "matCRE2subclass.csv").toString,
        overwrite = true,
        head = ""
    )
  }

  val outSignald = s"$projd/13.cicero/out"
  val eRNAfnms   = scs.map(sc => s"$outd/${sc}.eRNA.csv").toVector
  val eRNASignal = mergeBigWgSignals(eRNAfnms)
  outSignal(
      g = eRNASignal._1,
      sc = scs.toVector,
      m = eRNASignal._2,
      outd = os.Path(s"$outSignald/cicero_ppec_eRNA")
  )

  mods.foreach { h =>
    {
      println(s"prepare signal mat of $h")
      val fnms =
        scs.map(sc => s"$outd/${sc}.${h}.e100.bs100.sm300.csv").toVector
      val signals =
        mergeBigWgSignals(fnms, leftShift = 500, rightShift = -500)
      outSignal(
          g = signals._1,
          sc = scs.toVector,
          m = signals._2,
          outd = os.Path(s"$outSignald/cicero_ppec_$h")
      )
    }
  }

}
