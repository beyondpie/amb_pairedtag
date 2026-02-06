// The results based on peaks are weired.
// I am not sure if there are some bugs here.
// Currently, we are moving to genomic features, which,
// at least, makes more sense.
// commenter: szu
// time: 2025-06-19

import Math.*
import MetaData.readSubclassMeta
import SZUtils.{path2str, writeStrings2File}
import bioscala.LightCoord.BedGraph.*
import bioscala.LightCoord.BedGraph.LightBedGraphElement.mkString
import bioscala.LightCoord.GenomeCoord.{GenomeCoord, GenomeCoords}
import os.*
import project.Peak.loadAllModUnifiedPeak
import smile.math.MathEx
import smile.math.distance.EuclideanDistance

import scala.collection.parallel.*
import scala.collection.parallel.CollectionConverters.*
import scala.math.log1p

@main def Calsc2scSimUsingTopFeatures(): Unit = {
  val projd              = os.Path(MetaData.TSCCMeta.projd)
  val workd              = projd / "17.repressiveMarks"
  val outd               = workd / "out" / "allGenomicRange"
  val ptscMeta           = readSubclassMeta().filter(x => x.atac)
  val ptscs              = ptscMeta.map(x => x.name.PairedTagName)
  val mods               = Vector("H3K27me3", "H3K9me3")
  val emptyValue: Double = 0.0
  val topq: Double       = 0.05
  val lowq: Double       = 0.05

  // val mod: String = "H3K27me3"
  // val mod: String = "H3K9me3"
  val ntop: Int = 20000

  // functions
  def loadSubclassUnifiedPeakWigValue(
    mod: String): Map[String, Map[GenomeCoord, Double]] = {
    ptscs
      .map(sc => {
        val fnm =
          outd / s"subclass2unifiedPeakWigValue_$mod" / s"$sc.$mod.broadPeak.wigValue.bed"
        val r = os.read.lines
          .stream(fnm)
          .map(x => x.strip().split("\t"))
          .filter(x => x.nonEmpty)
          .map(x =>
            ((chrom = x.head,
                    coord = (startFrom = x(1).toInt,
                        endTo = x(2).toInt), strand = "."),
                x.last.toDouble))
          .toVector
          .toMap
        (sc, r)
      })
      .toMap
  }

  def getPeak2Value(sc2v: Map[String, Map[GenomeCoord, Double]],
    peak: Iterable[GenomeCoord], scs: Iterable[String],
    emptyValue: Double = 0.0): Vector[Vector[Double]] = {
    peak
      .map(p =>
        scs.map(sc => sc2v(sc).getOrElse(p, emptyValue)).toVector)
      .toVector
  }

  def getsc2Value(sc2v: Map[String, Map[GenomeCoord, Double]],
    peak: Iterable[GenomeCoord], scs: Iterable[String],
    emptyValue: Double = 0.0): Vector[Vector[Double]] = {
    scs
      .map(sc => {
        peak
          .map(p => {
            sc2v(sc).getOrElse(p, emptyValue)
          })
          .toVector
      })
      .toVector
  }

  def getTopUnifiedPeakWithScore(upeak: GenomeCoords,
    upeak2v: Vector[Vector[Double]], topq: Double = 0.05,
    lowq: Double = 0.05): LightBedGraph = {
    val npeak                  = upeak2v.size
    val upeakV: Vector[Double] = upeak2v.map(x => x.sum)
    upeakV.zipWithIndex
      .sortBy(x => -x._1)
      .zipWithIndex
      .filter((x, i) => i >= (npeak * topq).toInt)
      .filter((x, i) => i <= (npeak * (1 - lowq)).toInt)
      .map((x, _) => x._2)
      .map(i => (upeak(i), upeakV(i)))
  }

  def saveTopPeakWithScore(r: LightBedGraph,
    mod: String): Unit = {
    val outfnm = outd / s"topRankedUnifiedPeak_$mod.bed"
    writeStrings2File(
        content = r.map(x => x.mkString(sep = "\t")),
        to = outfnm,
        overwrite = true,
        head = ""
    )
  }

  def getsc2scSimMat(scs: Vector[String],
    sc2v: Vector[Vector[Double]],
    p: PairwiseClosenessOpt): Vector[Vector[Double]] = {
    val r = scs.indices
      .combinations(2)
      .map(xs => {
        val i = xs.head
        val j = xs.last
        val s = p.d(sc2v(i), sc2v(j))
        ((scs(i), scs(j)), s)
      })
      .toMap

    val sc2scMat: Vector[Vector[Double]] = scs.map(sc1 => {
      scs.map(sc2 => {
        if (sc1 == sc2) {
          p.toSelf
        } else if (r.contains((sc1, sc2))) {
          r((sc1, sc2))
        } else {
          r((sc2, sc1))
        }
      })
    })
    sc2scMat.map(x => p.s(x))
  }

  val simMap: Map[String, PairwiseClosenessOpt] = Map(
      "RawCor"            -> RawCor,
      "RawSpearmanCor"    -> RawSpearmanCor,
      "ReluCor"           -> ReluCor,
      "MinMaxCor"         -> MinMaxCor,
      "MinMaxRawEuc"      -> MinMaxRawEuclean,
      "MinMaxEucScaleLog" -> MinMaxEucleanScaleLog
  )

  List("H3K9me3", "H3K27me3").par.foreach(mod => {
    // 1. load subclass 2 wigvalue
    val sc2v: Map[String, Map[GenomeCoord, Double]] =
      loadSubclassUnifiedPeakWigValue(mod)

    // 2. load unified peaks.
    val allUnifiedPeakfnm = projd / "data" / "pairedtag_peak" /
      "merge_peak" / s"$mod.merged.all.blv2.me.peak"

    val allupeak = loadAllModUnifiedPeak(allUnifiedPeakfnm)

    // 3. get upeak'd wig values across all subclasses.
    val upeak2v = getPeak2Value(sc2v, allupeak, ptscs, emptyValue)

    // 4. select top features
    val topPeakWithScore =
      getTopUnifiedPeakWithScore(allupeak, upeak2v, topq, lowq)
    saveTopPeakWithScore(topPeakWithScore, mod)

    // 5. calulate subclass-level similarity
    val sc2selv = getsc2Value(sc2v, topPeakWithScore.map(_.g).take(ntop),
        ptscs, emptyValue)

    simMap.par.foreach((k, o) => {
      println(
          s"Start: sc2sc similarity mat for $mod with $k using $ntop features. ")
      val r = getsc2scSimMat(ptscs, sc2selv, o)
      writeStrings2File(
          content = r.map(x => x.mkString(",")),
          to = outd / s"sc2scSimMat.$mod.$k.$ntop.csv",
          overwrite = true,
          head = ptscs.mkString(sep = ",")
      )
      println(
          s"Done: sc2sc similarity mat for $mod with $k using $ntop features.")
    })

    // 6. intersect / join
    println(s"calculate Jaccard for $mod.")
    val sc2scJaccard = ptscs.map(sc1 => {
      ptscs.map(sc2 => {
        val sc1p = sc2v(sc1).keySet
        val sc2p = sc2v(sc2).keySet
        val nj = sc1p.intersect(sc2p).size
        val nu = sc1p.size + sc2p.size - nj
        nj.toDouble / nu.toDouble
      })
    })
    writeStrings2File(
      content = sc2scJaccard.map(x => x.mkString(",")),
      to = outd / s"sc2scSimMat.$mod.Jaccard.csv",
      overwrite = true,
      head = ptscs.mkString(sep = ",")
    )
    println(s"Done: calculate Jaccard for $mod.")
  })

} // end of main

trait PairwiseClosenessOpt {
  val toSelf: Double
  def d(x: Vector[Double], y: Vector[Double]): Double
  def s(x: Vector[Double]): Vector[Double]
}

object RawCor extends PairwiseClosenessOpt {
  override val toSelf: Double = 1.0

  override def d(x: Vector[Double], y: Vector[Double]): Double = {
    MathEx.cor(x.toArray, y.toArray)
  }

  override def s(x: Vector[Double]): Vector[Double] = x
}

object RawSpearmanCor extends PairwiseClosenessOpt {
  override val toSelf: Double = 1.0

  override def d(x: Vector[Double], y: Vector[Double]): Double = {
    MathEx.spearman(x.toArray, y.toArray)
  }

  override def s(x: Vector[Double]): Vector[Double] = x
}

object ReluCor extends PairwiseClosenessOpt {
  override val toSelf: Double = 1.0
  override def d(x: Vector[Double], y: Vector[Double]): Double = {
    MathEx.cor(x.toArray, y.toArray)
  }
  override def s(x: Vector[Double]): Vector[Double] = {
    x.map(y => relu(y))
  }
}

object MinMaxCor extends PairwiseClosenessOpt {
  override val toSelf: Double = 1.0
  override def d(x: Vector[Double], y: Vector[Double]): Double = {
    MathEx.cor(x.toArray, y.toArray)
  }

  override def s(x: Vector[Double]): Vector[Double] = {
    val m = x.min
    val M = toSelf
    minMaxScale(x, m = m, M = M)
  }
}

object MinMaxRawEuclean extends PairwiseClosenessOpt {
  override val toSelf: Double = 0.0
  override def d(x: Vector[Double], y: Vector[Double]): Double = {
    EuclideanDistance().d(x.toArray, y.toArray)
  }

  override def s(x: Vector[Double]): Vector[Double] = {
    val M = toSelf
    val m = x.max
    minMaxScale(x, m = m, M = M)
  }
}
object MinMaxEucleanScaleLog extends PairwiseClosenessOpt {
  override val toSelf: Double = 0.0
  override def d(x: Vector[Double], y: Vector[Double]): Double = {
    EuclideanDistance().d(scale(x.map(t => log1p(t))).toArray,
        scale(y.map(t => log1p(t))).toArray)
  }

  override def s(x: Vector[Double]): Vector[Double] = {
    val M = toSelf
    val m = x.max
    minMaxScale(x, m = m, M = M)
  }
}
