import scala.collection.immutable.HashMap
import os._
import PairedTag.ReadPairedTagCellMeta
import PairedTag.PairedTagCell
import MouseRegion.RegionMapping.PairedTagRegion2MajorRegion

object EnrichL5 {
  def main(args: Array[String]) = {

    val projd = os.pwd
    val metaFile = projd / "meta" / "pairedtag.cell.meta.all.with.tfv3.csv"
    val workd = projd / "03.integration"
    val uncoed = workd / "out" / "tfneu_vf_region_cca_k5" / "un_coembed"

    lazy val mCells: Seq[Option[PairedTagCell]] =
      ReadPairedTagCellMeta.readPairedTagCellMetaFile(fnm = metaFile)
    lazy val barcode2LQ: Map[String, String] =
      mCells.map(_.get).filter(x => x.annot.l5r.contains("LQ")).
        map(x => (x.barcode, x.annot.l5r)).
        toMap

    lazy val barcode2L5r: Map[String, String] = mCells.map(_.get).
      map(x => x.barcode -> x.annot.l5r).
      toMap
    lazy val L5rtoSize: Map[String, Int] = mCells.map(_.get).
      map(x => (x.barcode, x.annot.l5r)).
      groupMapReduce(_._2)(x => 1)((x, y) => x + y)

    val fnmZW = workd / "out" / "tfneu_vf_region_cca_k5" /
      "amb_PT_allcellmetadata.ZW.20240614.csv"

    lazy val barcode2IfKeep: Map[String, String] = os.read.lines.stream(fnmZW).
      slice(1, Int.MaxValue).
      toSeq.
      map(_.replace("\"", "")).
      map(x => x.split(",").head -> x.split(",").last).toMap

    val ptGlobalTFile = projd / "meta"/
      "pairedtag.cell.meta.all.with.init.tf.csv"

    lazy val barcode2annot = os.read.lines.stream(ptGlobalTFile).
      slice(1, Int.MaxValue).
      toSeq.
      map(x => x.split(",").head -> x.split(",").takeRight(2)).
      toMap

    lazy val L5toAnnot = os.read.lines.stream(ptGlobalTFile).
      slice(1, Int.MaxValue).
      toSeq.
      map(x => x.split(",").takeRight(4)).
      map(x => (x(0), x(2))).
      groupBy(_(0)).
      toList.
      map(x => s"L5-${x._1}" -> x._2(0)._2).
      toMap


    def enrich(barcodes: List[String], cutoff: Float = 0.8) = {
      barcodes.filter(x => barcode2IfKeep.contains(x)).
        filter(x => barcode2IfKeep(x) == "keep").
        map(x => (barcode2L5r(x), 1)).
        groupBy(_._1).
        toList.
        map(x => (x._1, x._2.length,
          x._2.length.toFloat / L5rtoSize(x._1).toFloat
        )).
        filter(x => x._3 > cutoff)
    }


    def outf(barcodes: List[String], fnm: os.Path, cutoff: Float = 0.75)  = {
      val u = barcodes.map(x => (x, barcode2L5r(x))).
        groupMapReduce(_._2)(x => 1)((x, y) => x + y).toMap.
        map((k, v) => (k, v, L5rtoSize(k))).
        map(x => (x._1, x._2, x._3, x._2.toFloat / x._3.toFloat))
      val v = u.filter(_._4 > cutoff).toList
      os.write(fnm, s"# ${barcodes.length} barcodes: ${u.size} L5r \n")
      v.foreach {
        x => os.write.append(fnm, s"${x._1},${x._2},${x._3},${x._4}\n")
      }
      u
    }


    // // * get enriched L5r.
    // // ** AMY region
    // val r = "AMY"
    // lazy val b = os.read.lines(
    //   uncoed / s"${r}.barcode.txt").toList
    // val outfnm = uncoed / s"${r}.to.be.deleted.txt"
    // val u = outf(barcodes, outfnm, cutoff = 0.8)
    
    // // ** CPU region
    // // [KEEP] barcode 1 
    val r:String = "CPU"
    // lazy val barcodes = os.read.lines(
    //     uncoed / s"${r}.barcode.1.txt").toList
    // // val outfnm = uncoed / s"${r}.enrich.L5r.1.txt"
    // // val u = outf(barcodes, outfnm, cutoff = 0.8)
    // val u = enrich(barcodes, cutoff = 0.8)
    // val outfnm2 = uncoed / s"${r}.enrich.L5r.with.annot.txt"
    // val uu = u.map(x => (x._1, x._2, L5toAnnot(x._1.split(":")(0))))
    // uu.foreach {
    //   x => os.write.append(outfnm2, s"${x.toList.mkString(",")}\n")
    // }

    // // barcode 2
    // lazy val barcodes = os.read.lines(
    //   uncoed / s"${r}.barcode.2.txt").toList
    // val outfnm = uncoed / s"${r}.to.be.deleted.2.txt"
    // // no enriched L5r
    // val u = outf(barcodes, outfnm, cutoff = 0.8)
    // // barcode 3
    lazy val b = os.read.lines(
      uncoed / s"${r}.barcode.3.txt").toList
    val outfnm = uncoed / s"${r}.to.be.deleted.3.txt"
    // no enriched L5r
    val u = outf(b, outfnm, cutoff = 0.8)

    // // ** HYP
    // val r = "HYP"
    // lazy val barcodes = os.read.lines(
    //    uncoed / s"${r}.barcode.txt").toList
    // val outfnm = uncoed / s"${r}.to.be.deleted.txt"

    // // ** HIP
    // val r = "HIP"
    // lazy val barcodes = os.read.lines(
    //   uncoed / s"${r}.barcode.1.txt").toList
    // val outfnm = uncoed / s"${r}.to.be.deleted.1.txt"
    // val u = outf(barcodes, outfnm, cutoff = 0.8)

    // lazy val barcodes = os.read.lines(
    //   uncoed / s"${r}.barcode.2.txt").toList
    // val outfnm = uncoed / s"${r}.to.be.deleted.2.txt"
    // val u = outf(barcodes, outfnm, cutoff = 0.8)
   
    // lazy val barcodes = os.read.lines(
    //   uncoed / s"${r}.barcode.3.txt").toList
    // val outfnm = uncoed / s"${r}.to.be.deleted.3.txt"
    // val u = outf(barcodes, outfnm, cutoff = 0.8)

    // ** ERC
    // val r = "ERC"
    // lazy val barcodes = os.read.lines(
    //    uncoed / s"${r}.barcode.txt").toList
    // // val outfnm = uncoed / s"${r}.enrich.L5r.txt"
    // // lazy val u = outf(barcodes, outfnm, cutoff = 0.8)
    // val u = enrich(barcodes, cutoff = 0.8)
    // val outfnm2 = uncoed / s"${r}.enrich.L5r.with.annot.txt"
    // val uu = u.map(x => (x._1, x._2, L5toAnnot(x._1.split(":")(0))))
    // uu.foreach {
    //   x => os.write.append(outfnm2, s"${x.toList.mkString(",")}\n")
    // }

    // // ** NAC
    // val r = "NAC"
    // lazy val barcodes = os.read.lines(
    //   uncoed / s"${r}.barcode.txt").toList
    // val outfnm = uncoed / s"${r}.to.be.deleted.txt"
    // lazy val u = outf(barcodes, outfnm, cutoff = 0.8)

    // // ** VTA
    // val r = "VTA"
    // lazy val barcodes = os.read.lines(
    //   uncoed / s"${r}.barcode.1.txt").toList
    // val outfnm = uncoed / s"${r}.to.be.deleted.1.txt"
    // lazy val u = outf(barcodes, outfnm, cutoff = 0.8)

    // lazy val barcodes = os.read.lines(
    //   uncoed / s"${r}.barcode.2.txt").toList
    // val outfnm = uncoed / s"${r}.to.be.deleted.2.txt"
    // lazy val u = outf(barcodes, outfnm, cutoff = 0.8)

    // // ** PFC
    // val r = "PFC"

    // lazy val barcodes = os.read.lines(
    //    uncoed / s"${r}.barcode.1.txt").toList
    // val outfnm = uncoed / s"${r}.to.be.deleted.1.txt"
    // lazy val u = outf(barcodes, outfnm, cutoff = 0.8)
 
    // lazy val barcodes = os.read.lines(
    //   uncoed / s"${r}.barcode.2.txt").toList
    // val outfnm = uncoed / s"${r}.to.be.deleted.2.txt"
    // lazy val u = outf(barcodes, outfnm, cutoff = 0.8)

    // lazy val barcodes = os.read.lines(
    //   uncoed / s"${r}.barcode.3.txt").toList
    // // no enrich
    // val u = enrich(barcodes, cutoff = 0.8)
    // val outfnm2 = uncoed / s"${r}.enrich.L5r.with.annot.txt"
    // val uu = u.map(x => (x._1, x._2, L5toAnnot(x._1.split(":")(0))))
    // uu.foreach {
    //   x => os.write.append(outfnm2, s"${x.toList.mkString(",")}\n")
    // }
  }
}
