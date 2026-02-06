import GRange.GenomicRange
import bioscala.LightCoord._
import SZUtils.{writeStrings2File, ifelse, str2path, path2str}
import os.*
import MetaData.PairedTagCluster
import bioscala.LightCoord.GenomeCoord.genomeCoordIgnoreStrandOrd
import bioscala.LightCoord.GenomeCoord.{GenomeCoord, GenomeCoords}
import bioscala.LightCoord.GenomeCoord.isOverlapIgnoreStrand
import bioscala.LightCoord.GenomeCoord.{
  mkStringGenomeCoord, toStringGenomeCoord
}
import bioscala.LightCoord.Bed.LightBed
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*

case class CRE(g: GenomeCoord, isDAR: String, subclass: String,
  chromHMMState: String, isDistal: String, cnsrv: String) {
  def mkString(sep: String = "\t"): String = {
    val x = mkStringGenomeCoord(g, sep, useStrand = false)
    List(x, isDAR, isDistal, chromHMMState, subclass, cnsrv)
      .mkString(sep)
  }
}

object CRE {
  def CREHeader(sep: String = "\t"): String = {
    List(
        "chrom",
        "startFrom",
        "endTo",
        "isDAR",
        "isDistal",
        "chromHMMState",
        "subclass",
        "convervation"
    ).mkString(sep)
  }
}

def annotCREWithChromHMMandDAR_(
  p: LightBed,
  a: LightBed,
  d: Set[String],
  sc: String,
  CAdiv: Set[String],
  CAcnsrv: Set[String]
): Vector[CRE] = {
  val isProximal: Vector[Boolean] =
    isOverlapIgnoreStrand(query = a.map(_.g),
        subject = p.map(_.g))
  a.zipWithIndex.map((x, i) => {
    val gstr = toStringGenomeCoord(x.g, useStrand = false)
    new CRE(
        g = x.g,
        isDAR = ifelse(d.contains(toStringGenomeCoord(x.g)),
            "DAR", "common"),
        subclass = sc,
        chromHMMState = x.name,
        isDistal = ifelse(isProximal(i), "proximal", "distal"),
        cnsrv = ifelse(CAdiv.contains(gstr), "CAdiv",
            ifelse(CAcnsrv.contains(gstr), "CAcnsrv", "mouse"))
    )
  })
}

@main def AnnotCRE_ChromHMMStateDARConserv(): Unit = {
  val projd: String = MetaData.TSCCMeta.projd
  val workd         = projd + "/12.DE"
  val DARfnm = projd + "/16.celloracle/out/" + "DAR.log2fd0.2.csv"
  val promoterfnm = "/projects/ps-renlab2/szu".+(
      "/genome/mouseGenePromoter.GRCm38-M25.promoter.1kbAwayTSS.tsv")
  val CREChromHMMAnnotd =
    projd + "/06.ChromHMM/out/CREAnnot250429"
  val CAcnsrvfnm =
    projd + "/meta/mba.hba.CAcons.orthInHg38.rmdup.bed"
  val CAdivfnm = projd + "/meta/mba.hba.CAdiver.orthInHg38.bed"

  val ptscMeta: Vector[PairedTagCluster] =
    MetaData.readSubclassMeta().filter(_.atac)

  val ptscs: Vector[String] = ptscMeta.map(_.name.PairedTagName)

  val outd = projd + "/data/chromHMM/CRE"
  if (!os.exists(os.Path(outd))) {
    os.makeDir(os.Path(outd))
  }

  val promoter: LightBed =
    os.read.lines
      .stream(os.Path(promoterfnm))
      .map(x => x.strip().split("\t"))
      .map(x =>
        (g = (chrom = x.head,
                coord = (startFrom = x(1).toInt,
                    endTo = x(2).toInt), strand = x(4)),
            name = x(3), score = 0.0))
      .toVector

  val DAR: Map[String, Set[String]] =
    os.read.lines
      .stream(os.Path(DARfnm))
      .drop(1)
      .map(x => x.strip().split(","))
      .map(x => (x.head, x.last))
      .toVector
      .groupMap(_._2)(_._1)
      .map((k, v) => (k, v.toSet))

  val CAdiv: Set[String] = os.read.lines
    .stream(CAdivfnm)
    .map(x => x.strip().split("\t"))
    .filter(x => x.nonEmpty)
    .filter(x => !Vector("chrX", "chrY").contains(x(0)))
    .map(_.last)
    .toSet

  val CAcnsrv: Set[String] = os.read.lines
    .stream(CAcnsrvfnm)
    .map(x => x.strip().split("\t"))
    .filter(x => x.nonEmpty)
    .filter(x => !Vector("chrX", "chrY").contains(x(0)))
    .map(_.last)
    .toSet

  ptscs.par.foreach(sc => {
    println(s"annotCRE for $sc.")
    val annot: LightBed =
      os.read.lines
        .stream(os.Path(
                CREChromHMMAnnotd + s"/$sc.CRE.annotBySummit.bed"))
        .map(x => x.strip().split("\t"))
        .map(x =>
          (g = (chrom = x.head,
                  coord = (startFrom = x(1).toInt,
                      endTo = x(2).toInt), strand = "."),
              name = x.last, score = 0.0))
        .toVector
    val r =
      annotCREWithChromHMMandDAR_(promoter, annot, DAR(sc), sc,
          CAdiv, CAcnsrv)
    val outfnm = outd + s"/$sc.CRE.ChromHMM.DAR.distal.annot.tsv"
    writeStrings2File(
        content = r.map(x => x.mkString("\t")),
        to = outfnm,
        overwrite = true,
        head = CRE.CREHeader(sep = "\t")
    )
  })
}
