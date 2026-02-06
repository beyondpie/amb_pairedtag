package project.Peak

import os._
import Peak.{Peak, Peaks}
import scala.collection.immutable.HashMap
import SZUtils.{hereByGit, getFileSize, TaskElement}
import SZUtils.{readTable, writeStrings2File}
import GRange.{GenomicRange, mouseGenomicRangeOrd}
import MetaData.TSCCMeta
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*
import bioscala.LightCoord.GenomeCoord.GenomeCoords
import SZUtils.{str2path, path2str, ifelse}
import bioscala.LightCoord.GenomeCoord.genomeCoordIgnoreStrandOrd

def peakd: os.Path = {
  TSCCMeta.rawProjd / "data" / "pairedtag_peak" / "subclass_peak"
}

def groupStat = {
  readTable(
    fnm =
      s"${TSCCMeta.projd}/04.peakcalling/src/main/resource/bam2nread.full.csv",
    sep = ",", head = false
  )
}


def group2nread: Map[String, Int] = {
  groupStat.map(
    x => {
      val suffix = x(1).split("\\.").head
      (s"${x(0)}-${suffix}", x(2).toInt)
    }
  ).toMap
}

// for lambda
val mod2minpeak = HashMap(
  ("H3K27ac", 100),
  ("H3K27me3", 1000),
  // 10^2.8
  ("H3K4me1", 630),
  ("H3K9me3", 630)
)

// for nolambda
// 10^4.5: 31622
val mod2minread = HashMap(
  ("H3K27ac", 10000),
  ("H3K27me3", 31622),
  ("H3K4me1", 31622),
  ("H3K9me3", 31622)
)

val mod2peakClass = HashMap(
  ("H3K27ac", "narrow"),
  ("H3K27me3", "broad"),
  ("H3K4me1", "broad"),
  ("H3K9me3", "broad")
)


val ignoredGroup: Map[String, String] = HashMap(
  ("021 L4 RSP-ACA Glut", "021_L4_RSP_ACA_Glut")
)

def ignoreGroupByPeak1(modality: String, p: Peaks) =
  p.length < mod2minpeak(modality)

def ignoreGroupByPeak2(peakfnm: String, modality: String) =
  val group = peakfnm.split("/")
    .last.split("_").toList.init.mkString("_")
  if group2nread.isDefinedAt(group) then
    group2nread(group) < mod2minread(modality)
  else false

class GroupOfPeaks(
  from: String,
  name: String,
  sex: String,
  modality: String,
  neglog10qval: Float = 2,
  outd: String,
  flagd: String,
  skipTask: Boolean = true
) extends TaskElement {
  val peakClass = mod2peakClass(modality)
  val prefix = List(name, modality, sex).mkString("-")
  val flagfnm = s"${flagd}/${prefix}.done"
  val outfnm = s"${outd}/${prefix}.reproPeak"
  val skip = skipTask
  val fnms: List[os.Path] =
    List(sex, s"${sex}A", s"${sex}B", s"${sex}-shufA", s"${sex}-shufB")
      .map(x => s"${name}-${modality}-${x}_peaks.${peakClass}Peak")
      .map(x => os.Path(from) / name / x)
      .filter(os.exists)
  def getTreatPeaks: Option[Peaks] =
    if (
      (fnms.length > 0) &&
      os.exists(fnms(0)) &&
      (getFileSize(fnms(0).toString) > 2)
    ) {
      val tmp = Peak
        .readBed(fnm = fnms(0).toString)
        .filter(_.neglog10qval >= neglog10qval)
      if tmp.length < 1 then None
      else Option(tmp.sortBy(x => (x.r.chrom, x.r.startFrom, x.r.endTo)))
    } else {
      None
    }

  def getReplicatePeaks: List[Map[String, List[GenomicRange]]] =
    if (fnms.length < 2) then Nil
    else
      fnms.tail.map(x => {
        if getFileSize(x.toString) < 2 then Map()
        else {
          val tmp = Peak.readBed(x.toString).toList
          //if ignoreGroupByPeak1(modality, tmp) then Map()
          if ignoreGroupByPeak2(x.toString, modality) then Map()
          else
            tmp
              .filter(_.neglog10qval >= neglog10qval)
              .map(_.r)
              .groupBy(_.chrom)
              .map((k, v) => (k, v.sorted.toList))
        }
      })

  def getReproPeak: Peaks =
    val optionTreatPeaks = getTreatPeaks
    if optionTreatPeaks.isEmpty then
      println(s"${prefix} has no peaks from treatPeaks.")
      Nil
    else
      val treatPeaks = optionTreatPeaks.get
      println(s"${prefix} has ${treatPeaks.length} peaks.")
      val repliPeaks = getReplicatePeaks
      if (repliPeaks.isEmpty) {
        treatPeaks
      } else {
        val tmp = treatPeaks
          .map(x =>
            (
              x,
              // use parallel
              repliPeaks.par
                .map(l => {
                  if l.isEmpty | x.r.isInSortedRanges(l) then 1 else 0
                })
                .sum
            )
          )
          .filter(_._2 >= (repliPeaks.length - 1))
        if tmp.length < 1 then Nil
        else tmp.map(_._1)
      }
  def runCore(): Unit = {
    println(s"Run reproPeak on: ${prefix}")
    val r = getReproPeak
    if r == Nil then
      println(s"No peaks from ${prefix} after getReproPeak.")
    else
      println(s"write file to ${outfnm}.")
      Peaks.write(r, fnm = outfnm, sep = "\t", overwrite = true)
  }
}


// TODO: double check the unified peaks.
def loadSubclassModUnifedPeak(sc: String, mod: String,
                              d: os.Path): GenomeCoords = {
  os.read.lines
    .stream(d / s"$sc.$mod.upeak.bed")
    .map(x => x.strip().split("\t"))
    .filter(x => x.length >= 3)
    .map(x =>
      (chrom = x.head, coord = (x(1).toInt, x(2).toInt),
        strand = "."))
    .toVector
}

def loadAllModUnifiedPeak(fnm: String): GenomeCoords = {
    os.read.lines
      .stream(fnm)
      .drop(1)
      .map(x => x.strip().split("\t"))
      .filter(_.nonEmpty)
      .filter(x =>
        !Set("chrx", "chry").contains(x.head.toLowerCase))
      .map(x =>
        (chrom = x.head,
          coord = (
            startFrom = x(1).toInt,
            endTo = x(2).toInt
          ), strand = "."))
      .toVector
      .sorted
}

/**
  * Returns subclass-enriched peaks (not unified) of
  * histone modificatiosn from the Paired-Tag data.
  *
  * @param sc PairedTag subclass name
  * @param mod One of four histone modifications
  * @param peakd os.Path of the directory, default
  * [[project.Peak.peakd()]].
  */
def loadPairedTagSubclassPeak(sc: String, mod: String,
  peakd: os.Path = peakd): GenomeCoords = {
  val fnm: os.Path = ifelse(
      mod == "H3K27ac",
      peakd / s"${sc}-${mod}.BestSPM.peak",
      peakd / s"${sc}-${mod}.bedtoolmerge.peak"
  )
  os.read.lines
    .stream(fnm)
    .drop(1)
    .map(x => x.strip().split("\t"))
    .filter(x => x.nonEmpty)
    .filter(x => x.length >= 3)
    .map(x =>
      (chrom = x.head, coord = (x(1).toInt, x(2).toInt),
          strand = "."))
    .toVector
}

