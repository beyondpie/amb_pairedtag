import os._
import MetaData.TSCCMeta.{
  projd,
  metad,
  AllenClusterMetaFile,
  ptMetaWithGlobalTFile,
  ptMetaWithRegionTFile,
  ptRegionUMAPFile,
  ptPerRegionCheckFile,
  ptRegionTFRawBarcode2ClFile
}
import MetaData.PairedTagQCMeta.{unCoEmbedLQFiles, unCoEmbedHQFiles}
import PairedTag.{BrainRegion, Doublet, Experiment}
import PairedTag.{MappingQuality, Quality, Sample}
import MouseRegion.RegionMapping.PairedTagRegion2MajorRegion as r2mr

case class UMAP(x: Double, y: Double)
case class Cluster(id: Int, umap: Option[UMAP])
case class IterCluster(
    l1: Cluster,
    l2: Cluster,
    l3: Cluster,
    l4: Cluster,
    l5: Cluster
)

// g: global; r: region; f: final
// for g and r: cl is single-cell level
// for f: cl is L5 / L5r level
case class PerformedAnnot(
    cl: Option[String],
    sp: String
)

case class SimpleBarcodeMeta(
    barcode: String,
    sample: Sample,
    quality: Quality,
    cluster: IterCluster,
    exp: Experiment,
    isNeuL1: Boolean,
    l5r: Option[String]
)

// here annot is final supertype
case class PairedTagBarcodeMeta(
    barcodeMeta: SimpleBarcodeMeta,
    ptcl: String,
    regionUMAP: Option[UMAP],
    gannot: PerformedAnnot,
    rannot: Option[PerformedAnnot],
    annot: Option[String],
    annotQuality: String
)

lazy val allenspMeta: Map[String, (String, String)] =
  os.read.lines
    .stream(os.Path(AllenClusterMetaFile))
    .slice(1, Int.MaxValue)
    .toSeq
    .map(x => x.replace("\"", ""))
    .map(x => x.split("\t"))
    .map(x => x(5) -> (x(8), x(11)))
    .toMap

def showOption[A](x: Option[A]) = x match {
  case Some(s) => s.toString
  case None    => "NA"
}

def showOptionUMAP(x: Option[UMAP]) = x match {
  case Some(UMAP(x, y)) => s"${x.toString},${y.toString}"
  case None             => "NA,NA"
}

def showOptionRAnnot(x: Option[PerformedAnnot]) = x match {
  case Some(PerformedAnnot(Some(x), y)) =>
    s"$x,$y"
  case Some(PerformedAnnot(None, y)) =>
    s"NA,$y"
  case _ => s"NA,NA"
}

def showOptionFAnnot(x: Option[String]) = x match {
  case Some(x) =>
    s"$x,${allenspMeta(x)._1},${allenspMeta(x)._2}"
  case None => "NA,NA,NA"
}

object updatePairedtagMeta {
  def main(args: Array[String]) = {
    lazy val allenClMeta: Map[String, (String, String)] =
      os.read.lines
        .stream(os.Path(AllenClusterMetaFile))
        .slice(1, Int.MaxValue)
        .toSeq
        .map(x => x.replace("\"", ""))
        .map(x => x.split("\t"))
        .map(x => x(5) -> (x(8), x(11)))
        .toMap

    lazy val barcodeRegionUMAP: Map[String, UMAP] =
      os.read.lines
        .stream(os.Path(ptRegionUMAPFile))
        .slice(1, Int.MaxValue)
        .map(x => x.split(","))
        .map(x => x(2) -> UMAP(x(0).toDouble, x(1).toDouble))
        .toSeq
        .toMap

    lazy val barcodeLQClassUMAPCheck: Set[String] =
      os.read.lines
        .stream(os.Path(ptPerRegionCheckFile))
        .slice(1, Int.MaxValue)
        .map(x => x.replace("\"", ""))
        .filter(_.contains("remove"))
        .map(x => x.split(",").head)
        .toSeq
        .toSet

    lazy val barcodeLQL5rNum: Set[String] =
      os.read.lines
        .stream(os.Path(ptMetaWithRegionTFile))
        .slice(1, Int.MaxValue)
        .filter(_.contains("LQ_tf_by_region"))
        .map(x => x.split(",").head)
        .toSet

    lazy val L5rLQTFScore: Set[String] =
      unCoEmbedLQFiles
        .map(x => os.read.lines.stream(os.Path(x)))
        .map(x => x.slice(1, Int.MaxValue))
        .map(x => x.toList)
        .flatten
        .map(x => x.split(",").head)
        .toSet

    lazy val L5r2spHQUnCoEmbed: Map[String, String] =
      unCoEmbedLQFiles
        .map(x => os.read.lines.stream(os.Path(x)))
        .map(x => x.slice(1, Int.MaxValue))
        .map(x => x.toList)
        .flatten
        .map(_.split(","))
        .map(x => x.head -> x(3))
        .toMap

    lazy val barcode2ga: Map[String, PerformedAnnot] =
      os.read.lines
        .stream(os.Path(ptMetaWithGlobalTFile))
        .slice(1, Int.MaxValue)
        .toSeq
        .map(x => x.split(","))
        .map(x =>
          x.takeRight(3)(0) match
            case "" => (x.head, PerformedAnnot(None, x.takeRight(3)(1)))
            case _ =>
              (
                x.head,
                PerformedAnnot(
                  Some(x.takeRight(3)(0)),
                  x.takeRight(3)(1)
                )
              )
        )
        .map(x => x._1 -> x._2)
        .toMap
    // per-region transfer label raw barcode 2 cl
    // single-cell level (cells after downsampling)
    lazy val rtfBarcode2cl: Map[String, String] =
      os.read.lines
        .stream(os.Path(ptRegionTFRawBarcode2ClFile))
        .slice(1, Int.MaxValue)
        .toSeq
        .map(x => x.split(","))
        .map(x => x(0) -> x(1))
        .toMap

    lazy val barcode2ra: Map[String, Option[PerformedAnnot]] =
      os.read.lines
        .stream(os.Path(ptMetaWithRegionTFile))
        .slice(1, Int.MaxValue)
        .toSeq
        .map(x => x.split(","))
        .map(x =>
          if (x.takeRight(2).head == "LQ_tf_by_region")
            x(0) -> None
          else if (x.takeRight(2).head == "FALSE")
            x(0) -> None
          else if (rtfBarcode2cl.contains(x(0)))
            x(0) -> Some(
              PerformedAnnot(
                Some(x.takeRight(6).head), x.takeRight(5).head))
          else
            x(0) -> Some(PerformedAnnot(None, x.takeRight(5).head))
        )
        .toMap

    // only filtered by LQ_tf_by_region
    lazy val barcode2raFilter: Map[String, PerformedAnnot] =
      os.read.lines
        .stream(os.Path(ptMetaWithRegionTFile))
        .slice(1, Int.MaxValue)
        .filter(!_.contains("LQ_tf_by_region"))
        .toSeq
        .map(x => x.split(","))
        .map(x =>
          x(0) -> PerformedAnnot(
            Some(x.takeRight(6).head), x.takeRight(5).head)
        )
        .toMap

    lazy val barcode2Meta: List[SimpleBarcodeMeta] =
      os.read.lines
        .stream(os.Path(ptMetaWithRegionTFile))
        .slice(1, Int.MaxValue)
        .toList
        .map(x => x.split(","))
        .map(x =>
          SimpleBarcodeMeta(
            x(0),
            Sample(x(11), x(12), BrainRegion(x(7), r2mr(x(7)))),
            Quality(
              MappingQuality(
                x(4).toInt,
                x(5).toInt,
                if (x(6).length < 1) 0.toFloat else x(6).toFloat,
                if (x(14).length < 1) 0.toFloat else x(14).toFloat
              ),
              Doublet(
                x(1).toFloat,
                x(2).toFloat,
                if (x(16) == "FALSE") false else true
              )
            ),
            IterCluster(
              Cluster(x(19).toInt, Some(UMAP(x(17).toDouble, x(18).toDouble))),
              Cluster(x(22).toInt, Some(UMAP(x(20).toDouble, x(21).toDouble))),
              Cluster(x(26).toInt, Some(UMAP(x(24).toDouble, x(25).toDouble))),
              Cluster(
                x(30).toInt,
                if (x(28).length < 1) None
                else Some(UMAP(x(28).toDouble, x(29).toDouble))
              ),
              Cluster(
                x(34).toInt,
                if (x(32).length < 1) None
                else Some(UMAP(x(32).toDouble, x(33).toDouble))
              )
            ),
            Experiment(x(13), x(10), x(8)),
            if (x(40) == "FALSE") false else true,
            if (x(39) == "LQ_tf_by_region") None
            else Some(x(39))
          )
        )
    lazy val barcode2FullMeta: List[PairedTagBarcodeMeta] =
      barcode2Meta
        .map(x =>
          PairedTagBarcodeMeta(
            x,
            List(
              x.cluster.l1.id,
              x.cluster.l2.id,
              x.cluster.l3.id,
              x.cluster.l4.id,
              x.cluster.l5.id
            ).mkString("-"),
            barcodeRegionUMAP.get(x.barcode),
            barcode2ga(x.barcode),
            barcode2ra(x.barcode),
            if (x.l5r.isEmpty)
              None
            else if (L5r2spHQUnCoEmbed.contains(x.l5r.get))
              Some(barcode2ga(x.barcode).sp)
            else if (barcode2raFilter.contains(x.barcode))
              Some(barcode2raFilter(x.barcode).sp)
            else
              None,
            if (barcodeLQClassUMAPCheck.contains(x.barcode))
              "LQ_UMAP_disperse_at_class_per_region"
            else if (barcodeLQL5rNum.contains(x.barcode))
              "LQ_L5r_low_number"
            else if (x.l5r.isDefined && L5rLQTFScore.contains(x.l5r.get))
              "LQ_L5r_low_transferlabel_score"
            else
              "Good"
          )
        )
    // output
    val outptCellMetaFile = os.Path(metad) /
    "pairedtag.cell.meta.all.240626.csv"
    if (os.exists(outptCellMetaFile))
      println(s"${outptCellMetaFile} exists. Remove it.")
      os.remove(outptCellMetaFile)

    val head = List(
      "barcode",
      "sample.sex",
      "sample.rep",
      "sample.region.name",
      "sample.region.major",
      "mapping.nRNA",
      "mapping.nDNA",
      "maping.prcntMT",
      "mappiing.prcntMouseReads",
      "doublet.score",
      "doublet.prob",
      "doublet.isdlt",
      "cluster.l1.id",
      "cluster.l1.umap.x",
      "cluster.l1.umap.y",
      "cluster.l2.id",
      "cluster.l2.umap.x",
      "cluster.l2.umap.y",
      "cluster.l3.id",
      "cluster.l3.umap.x",
      "cluster.l3.umap.y",
      "cluster.l4.id",
      "cluster.l4.umap.x",
      "cluster.l4.umap.y",
      "cluster.l5.id",
      "cluster.l5.umap.x",
      "cluster.l5.umap.y",
      "exp",
      "exp.sublib",
      "exp.modularity",
      "isNeuL1",
      "l5r",
      "pairedtagCluster",
      "regionUMAP.x",
      "regionUMAP.y",
      "gannot.cl",
      "gannot.sp",
      "rannot.cl",
      "rannot.sp",
      "annot.sp",
      "annot.sc",
      "annot.c",
      "annotQuality"
    ).mkString(",")

    os.write(outptCellMetaFile, s"${head}\n")
    lazy val result = barcode2FullMeta
      .map(x =>
        List(
          x.barcodeMeta.barcode,
          x.barcodeMeta.sample.sex,
          x.barcodeMeta.sample.rep,
          x.barcodeMeta.sample.region.name,
          x.barcodeMeta.sample.region.marjor,
          x.barcodeMeta.quality.mapping.nRNA,
          x.barcodeMeta.quality.mapping.nDNA,
          x.barcodeMeta.quality.mapping.prcntMT,
          x.barcodeMeta.quality.mapping.prcntMouseReads,
          x.barcodeMeta.quality.dlt.score,
          x.barcodeMeta.quality.dlt.prob,
          x.barcodeMeta.quality.dlt.isdlt,
          x.barcodeMeta.cluster.l1.id,
          showOptionUMAP(x.barcodeMeta.cluster.l1.umap),
          x.barcodeMeta.cluster.l2.id,
          showOptionUMAP(x.barcodeMeta.cluster.l2.umap),
          x.barcodeMeta.cluster.l3.id,
          showOptionUMAP(x.barcodeMeta.cluster.l3.umap),
          x.barcodeMeta.cluster.l4.id,
          showOptionUMAP(x.barcodeMeta.cluster.l4.umap),
          x.barcodeMeta.cluster.l5.id,
          showOptionUMAP(x.barcodeMeta.cluster.l5.umap),
          x.barcodeMeta.exp.exp,
          x.barcodeMeta.exp.sublib,
          x.barcodeMeta.exp.modularity,
          x.barcodeMeta.isNeuL1,
          showOption(x.barcodeMeta.l5r),
          x.ptcl,
          showOptionUMAP(x.regionUMAP),
          showOption(x.gannot.cl),
          x.gannot.sp,
          showOptionRAnnot(x.rannot),
          showOptionFAnnot(x.annot),
          x.annotQuality
        )
      ).map(
        x => x.mkString(",")
      ).mkString("\n")
    os.write.append(outptCellMetaFile, result)
  }
}
