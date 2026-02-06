package MetaData

/**
  * The module is designed specifically to the mouse brain PariedTag project.
  * Without any considerations of re-uses in other projects or in general.
  * 
  * Add note @ 2025-02-04.
  */

import os._
import scala.collection.immutable.HashMap
import GRange.{GenomicRange, mouseGenomicRangeOrd}
import SZUtils.{TaskElement, toStringOption, readTable}

case class BrainRegion(name: String, marjor: String)
case class Doublet(score: Float, prob: Float, isdlt: Boolean)
case class Experiment(exp: String, sublib: String, modularity: String)
case class Sample(sex: String, rep: String, region: BrainRegion)
case class MappingQuality(
  nRNA: Int,
  nDNA: Int,
  prcntMT: Float,
  prcntMouseReads: Float
)
case class Quality(mapping: MappingQuality, dlt: Doublet)
case class UMAP(x: Option[Double], y: Option[Double])
case class Cluster(id: Int, umap: UMAP)
case class IterCluster(
  l1: Cluster,
  l2: Cluster,
  l3: Cluster,
  l4: Cluster,
  l5: Cluster
)
case class PerformedAnnot(
  cl: Option[String],
  sp: Option[String]
)

case class FinalAnnot(
  sp: Option[String],
  sc: Option[String],
  c: Option[String]
)

case class PairedTagBarcodeMeta(
  barcode: String,
  sample: Sample,
  quality: Quality,
  cluster: IterCluster,
  exp: Experiment,
  isNeuL1: Boolean,
  l5r: Option[String],
  ptcl: String,
  regionUMAP: UMAP,
  gannot: PerformedAnnot,
  rannot: PerformedAnnot,
  annot: FinalAnnot,
  annotQuality: String
)


object PairedTagBarcodeMeta {
  def apply(csvline: String): PairedTagBarcodeMeta = {
    val t = csvline.split(",")
    new PairedTagBarcodeMeta(
      t(0),
      Sample(t(1), t(2), BrainRegion(t(3), t(4))),
      Quality(
        MappingQuality(t(5).toInt, t(6).toInt, t(7).toFloat, t(8).toFloat),
        Doublet(t(9).toFloat, t(10).toFloat, t(11).toBoolean)
      ),
      IterCluster(
        Cluster(t(12).toInt, UMAP(t(13).toDoubleOption, t(14).toDoubleOption)),
        Cluster(t(15).toInt, UMAP(t(16).toDoubleOption, t(17).toDoubleOption)),
        Cluster(t(18).toInt, UMAP(t(19).toDoubleOption, t(20).toDoubleOption)),
        Cluster(t(21).toInt, UMAP(t(22).toDoubleOption, t(23).toDoubleOption)),
        Cluster(t(24).toInt, UMAP(t(25).toDoubleOption, t(26).toDoubleOption))
      ),
      Experiment(t(27), t(28), t(29)),
      t(30).toBoolean,
      toStringOption(t(31)),
      t(32),
      UMAP(t(33).toDoubleOption, t(34).toDoubleOption),
      PerformedAnnot(toStringOption(t(35)), toStringOption(t(36))),
      PerformedAnnot(toStringOption(t(37)), toStringOption(t(38))),
      FinalAnnot(
        toStringOption(t(39)),
        toStringOption(t(40)),
        toStringOption(t(41))
      ),
      t(42)
    )
  }
  def readPairedTagCellMetaFile(
    path: String = TSCCMeta.ptCellMetaFile
  ): List[PairedTagBarcodeMeta] = {
    os.read.lines
      .stream(os.Path(path))
      .slice(1, Int.MaxValue)
      .toList
      .map(x => PairedTagBarcodeMeta(x))
  }
}

def isNeu(b: PairedTagBarcodeMeta): Boolean = {
  if b.isNeuL1 then true
  else if b.annot.sc.get.contains("IMN") then true
  else false
}

def toStringOfAllenIdLabel(s: String): String = {
  s.replace("/", "_").replaceAll(" +", "_").replaceAll("-", "_")
}



object PairedTagSublib{
  case class Sublib(expid: String, sublibid: String, workd: String,
    idDNA: String, idRNA: String, name: String, bamDNA: String,
    bamRNA: String)

  case class SublibFASTQ(sublibid: String, tag: String, path: String)
  case class SublibRawFASTQ(sublibid: String, tag: String, R1: String, R2: String)

  def readSublib(f:String): List[Sublib] = {
    readTable(f, ",", head = true)
      .map(x =>
        new Sublib(x(0), x(1), x(2), x(3),
            x(4), x(5), x(6), x(7)))
  }
}
