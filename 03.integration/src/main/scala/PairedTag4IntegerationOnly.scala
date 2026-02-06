package PairedTag

import scala.collection.immutable.HashMap
import scala.collection.immutable.Seq
import scala.math._
import MouseRegion.RegionMapping
import os._

case class UMAP(umap1: Option[Double], umap2: Option[Double])
case class Cluster(id: Int, umap: UMAP)
case class IterCluster(
    l1: Cluster,
    l2: Cluster,
    l3: Cluster,
    l4: Cluster,
    l5: Cluster
)
case class Annotation(
    l5r: String,
    cl: String,
    spIdLabel: String,
    scIdLabel: String,
    classIdLabel: String
)
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

case class PairedTagCell(
    barcode: String,
    sample: Sample,
    quality: Quality,
    cluster: IterCluster,
    annot: Annotation,
    exp: Experiment,
    isNeuFromL1: Boolean
)

object ReadPairedTagCellMeta {
  def readPairedTagCellMetaLine(s: String, pattern: String): PairedTagCell = {
    val x = s.split(pattern)
    PairedTagCell(
      x(0),
      Sample(
        x(11),
        x(12),
        BrainRegion(x(7), RegionMapping.PairedTagRegion2MajorRegion(x(7)))
      ),
      Quality(
        MappingQuality(x(4).toInt, x(5).toInt,
          if (x(6).length < 1) 0.toFloat else x(6).toFloat,
          if (x(14).length < 1) 0.toFloat else x(14).toFloat),
        Doublet(
          x(1).toFloat,
          x(2).toFloat,
          if (x(16) == "FALSE") false else true
        )
      ),
      IterCluster(
        Cluster(x(19).toInt, UMAP(Some(x(17).toDouble), Some(x(18).toDouble))),
        Cluster(x(22).toInt, UMAP(Some(x(20).toDouble), Some(x(21).toDouble))),
        Cluster(x(26).toInt, UMAP(Some(x(24).toDouble), Some(x(25).toDouble))),
        Cluster(
          x(30).toInt,
          UMAP(
            if (x(28).length < 1) None else Some(x(28).toDouble),
            if (x(29).length < 1) None else Some(x(29).toDouble)
          )
        ),
        Cluster(
          x(34).toInt,
          UMAP(
            if (x(32).length < 1) None else Some(x(32).toDouble),
            if (x(33).length < 1) None else Some(x(33).toDouble)
          )
        )
      ),
      Annotation(x(39), x(36), x(37), x(38),
        if(x(39) == "LQ_tf_by_region") "" else x(41)),
      Experiment(x(13), x(10), x(8)),
      if (x(40) == "FALSE") false else true
    )
  }

  def readPairedTagCellMetaFile(
      fnm: os.Path,
      header: Boolean = true,
      startLine: Int = 1,
      endLine: Int = Int.MaxValue,
      pattern: String = ","
  ): Seq[Option[PairedTagCell]] = {
    if (!os.exists(fnm)) {
      Seq(None)
    } else {
      val stream = if (header) {
        os.read.lines.stream(fnm).slice(max(1, startLine), endLine)
      } else {
        os.read.lines.stream(fnm).slice(max(0, startLine), endLine)
      }
      stream.map(x => Some(readPairedTagCellMetaLine(x, pattern))).toSeq
    }
  }
}
