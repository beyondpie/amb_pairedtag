// Provide data for consensus matrix used in Supplementary Figure.
// Here we focus on the L4 or L5-level clustering we have made, and their
// conresponding to subclass (neuronal cells) or supertype (non-neuronals)
// annotations.

// x-axis: our clusters
// y-axis: Allen's clusters (sc or sp) we have finally.

import os._
import MetaData.PairedTagBarcodeMeta
import SZUtils.writeStrings2File

object PrepareFinalConsensusMat {

  case class ConsensusMat(rownames: Vector[String], colnames: Vector[String], x: Vector[Vector[Double]]) {
    def write(d: String): Unit = {
      if (!os.exists(os.Path(d))) {
        os.makeDir(os.Path(d))
      }
      writeStrings2File(content = rownames.toList,
          to = s"${d}/rownames.txt", overwrite = true, head = "")
      writeStrings2File(content = colnames.toList,
          to = s"${d}/colnames.txt", overwrite = true, head = "")
      writeStrings2File(
          content =
            x.map(y => y.map(yy => yy.toString()).mkString(",")).toList,
          to = s"${d}/mat.csv",
          overwrite = true,
          head = ""
      )
    }
  }

  def getClusterName(x: PairedTagBarcodeMeta, clevel: String = "L3"): String = clevel match {
    case "L3" =>
      List(x.cluster.l1.id, x.cluster.l2.id, x.cluster.l3.id).mkString(
          "-")
    case "L4" =>
      List(x.cluster.l1.id, x.cluster.l2.id, x.cluster.l3.id,
          x.cluster.l4.id).mkString("-")
    case _ => x.annot.c.get
  }

  def getAnnotName(x: PairedTagBarcodeMeta, annotLevel: String = "sc"): String = annotLevel match {
    case "sc" => x.annot.sc.get
    case "sp" => x.annot.sp.get
    case _    => x.annot.c.get
  }

  def getConsensusMat(clevel: String = "L3", cellMeta: List[PairedTagBarcodeMeta], annotlevel: String = "sc"): ConsensusMat = {

    // Allen BICCN names for whole mouse brain
    // sorted
    val annotNames = cellMeta
      .map(x => getAnnotName(x, annotlevel))
      .distinct
      .toVector
      .sortBy(x => x.split(" ").head.toInt)

    val annot2Index: Map[String, Int] = annotNames.zipWithIndex.toMap

    val r: Map[String, Map[String, Int]] = cellMeta
      .groupBy { x => getClusterName(x, clevel) }
      .map((cl, xs) =>
        (cl,
            xs.groupBy { x => getAnnotName(x, annotlevel) }
              .map((annot, xx) => (annot, xx.length))))

    // cluster names ordered by annotNames with highest frequency
    val clusterNames = r
      .map((cl, annots) => {
        val annotId = annots
          .map((annot, n) => (annot2Index(annot), n))
          .maxBy(_._2)
          ._1
        (cl, annotId)
      })
      .toVector
      .sortBy(_._2)
      .map(_._1)

    val x = clusterNames.map(rnm => {
      val t      = r(rnm)
      val v      = annotNames.map(cnm => t.getOrElse(cnm, 0).toDouble)
      val ntotal = v.sum.toDouble
      v.map(x => x / ntotal).toVector
    })

    ConsensusMat(clusterNames, annotNames, x)
  }

  def combineBlocks(x: ConsensusMat, y: ConsensusMat): ConsensusMat = {
    val extx = x.x.map(v => v ++ Vector.fill(y.colnames.length)(0.0))
    val exty = y.x.map(v => Vector.fill(x.colnames.length)(0.0) ++ v)
    val X    = extx ++ exty
    val colnames = x.colnames ++ y.colnames
    val rownames = x.rownames ++ y.rownames
    ConsensusMat(rownames, colnames, X)
  }

  def main(args: Array[String]) = {
    val projd = "/tscc/projects/ps-renlab2/szu/projects/amb_pairedtag"
    val workd = s"${projd}/03.integration"
    val outd  = s"${workd}/out"

    val cellMeta: List[PairedTagBarcodeMeta] =
      PairedTagBarcodeMeta
        .readPairedTagCellMetaFile()
        .filter(x => x.annotQuality == "Good")

    val neuscs: List[String] =
      cellMeta
        .filter(x => x.isNeuL1)
        .filter(x => x.annot.sc.isDefined)
        .map(x => x.annot.sc.get)
        .distinct

    val nnsps: List[String] =
      cellMeta
        .filter(x => !x.isNeuL1)
        .filter(x => x.annot.sp.isDefined)
        .map(x => x.annot.sp.get)
        .distinct

    for (i <- List("L3", "L4", "L5")) {
      val neuConMat =
        getConsensusMat(i, cellMeta.filter(_.isNeuL1), "sc")
      neuConMat.write(s"${outd}/consensusMat_neu_$i")

      val nnConMat =
        getConsensusMat(i, cellMeta.filter(!_.isNeuL1), "sp")
      nnConMat.write(s"${outd}/consensusMat_nn_$i")

      val conMat = combineBlocks(neuConMat, nnConMat)
      conMat.write(s"$outd/consensusMat_all_$i")
    }
  } // end of main
}   // end of object
