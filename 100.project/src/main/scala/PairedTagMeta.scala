package MetaData
import scala.collection.immutable.HashMap
import GRange.{GenomicRange, mouseGenomicRangeOrd}
import os._

object TSCCMeta {
  val modality = List("H3K4me1", "H3K9me3", "H3K27ac", "H3K27me3")
  val sex      = List("Male", "Female")
  val rep      = List("MaleA", "MaleB", "FemaleA", "FemaleB")

  val projd: String =
    "/tscc/projects/ps-renlab2/szu/projects/amb_pairedtag"

  val rawProjd =
    os.Path("/projects/ps-renlab2/szu/projects/amb_pairedtag")

  val metad: String  = s"${projd}/meta"
  val tfoutd: String =
    s"${projd}/03.integration/out/tfneu_vf_region_cca_k5"
  val ptscDNAbamd: String =
    s"${projd}/02.track/out/cell_bams"
  val ptscRNAbamd: String =
    s"${projd}/02.track/out_RNA/cell_bams"

  val ptclusterDNAbamd: String =
    s"${projd}/data/ptDNAbam"
  val ptclusterRNAbamd: String =
    s"${projd}/data/ptRNAbam"

  val AllenClusterMetaFile: String =
    s"${metad}/AIT21_annotation_freeze_081523.tsv"

  val ptPerRegionCheckFile: String =
    s"${tfoutd}/amb_PT_allcellmetadata.ZW.20240614.csv"

  val ptRegionUMAPFile: String =
    s"${tfoutd}/pt.umap.per.region.csv"
  val ptRawMetaFile: String =
    s"${metad}/pairedtag.cell.meta.all.csv"

  val ptMetaWithGlobalTFile: String =
    s"${metad}/pairedtag.cell.meta.all.with.init.tf.csv"
  val ptMetaWithRegionTFile: String =
    s"${metad}/pairedtag.cell.meta.all.with.tfv3.csv"
  val ptCellMetaFile: String =
    s"${metad}/pairedtag.cell.meta.all.240626.csv"

  val ptRegionTFRawBarcode2ClFile: String =
    s"${tfoutd}/pt.raw.barcode2cl.tfbyregion.csv"

  val ptscMetafnm =
    s"$projd/meta/PairedTagSubclassMetaFromCellMetaChromHMM.csv"

  val NeuscSexRepBamStatFile =
    s"${projd}/04.peakcalling/src/main/resource/Neu.subclass.SexRep.stat.csv"
  val NNspSexRepBamStatFile =
    s"${projd}/04.peakcalling/src/main/resource/NN.supertype.SexRep.stat.csv"

  val mod2peakClass = HashMap(
      ("H3K27ac", "narrow"),
      ("H3K27me3", "broad"),
      ("H3K4me1", "broad"),
      ("H3K9me3", "broad")
  )

  def getBarcode2Bam(
    cellMeta: List[PairedTagBarcodeMeta]
  ): Map[String, String] = {
    cellMeta
      .map(x =>
        x.barcode -> {
          List(
              ptscDNAbamd,
              List(x.sample.region.name, x.exp.modularity,
                  x.sample.rep)
                .mkString("_"),
              s"${x.barcode}.srt.rmdup.bam"
          ).mkString("/")
        })
      .toMap
  }

  def getRawBamStat(bamStatFile: String) = {
    os.read.lines
      .stream(
          os.Path(bamStatFile)
      )
      .slice(1, Int.MaxValue)
      .toList
      .map(x => x.split(","))
      .groupBy(x => (x(0), x(1)))
      .map((k, v) =>
        (
            k,
            (
                v.map(_(3).toInt).fold(0)((x, y) => x + y),
                v.map(_(4).toInt).fold(0)((x, y) => x + y)
            )
        ))
      .map((k, v) =>
        (k, (v(0), v(1), v(1).toFloat / v(0).toFloat)))
  }

  def readSexRepBamStat(statfnm: String) = {
    os.read.lines
      .stream(os.Path(statfnm))
      .slice(1, 9999)
      .toList
      .map(x => x.split(","))
      .map(x =>
        List(
            toStringOfAllenIdLabel(x(0)),
            x(1),
            x(2),
            x(3),
            x(4)
        ))
  }

  def getSexBamStat(statfnm: String, pattern: String = "A|B") = {
    val sexreg = pattern.r
    val t      = readSexRepBamStat(statfnm)
    t.map(x =>
      List(
          x(0),
          x(1),
          sexreg.replaceAllIn(x(2), ""),
          x(3),
          x(4)
      ))
      .groupMapReduce(x => x.slice(0, 3))(x =>
        (x(3).toInt, x(4).toInt))((x, y) =>
        (x._1 + y._1, x._2 + y._2))
      .map((k, v) =>
        (k, (v._1, v._2, v._2.toFloat / v._1.toFloat)))
  }

} // end of TSCCMeta
