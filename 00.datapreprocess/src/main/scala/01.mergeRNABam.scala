import os._
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*

import SZUtils.TaskElement
import MetaData.{TSCCMeta, PairedTagBarcodeMeta}
import MetaData.toStringOfAllenIdLabel
import SZUtils.writeStrings2File

object MergeRNABam {
  val myworkd = s"${TSCCMeta.projd}/data/ptRNAbam"
  val myflagd = s"${myworkd}/flag"
  val logd = s"${myworkd}/log"
  val outd = s"${myworkd}/bam"
  val blistd = s"${myworkd}/blist"
  val scRNAbamd = TSCCMeta.ptscRNAbamd
  val scMetafnm = TSCCMeta.ptCellMetaFile
  val scMeta: List[PairedTagBarcodeMeta] =
    PairedTagBarcodeMeta
      .readPairedTagCellMetaFile()
      .filter(x => x.annotQuality == "Good")

  extension (a: PairedTagBarcodeMeta) {
    def scTostr: String = {
      if (a.annot.sc.isDefined) {
        toStringOfAllenIdLabel(a.annot.sc.get)
      } else {
        ""
      }
    }
    def toRNABam: String = {
      val r = a.sample.region.name
      val rep = a.sample.rep
      val t = s"${r}_${rep}_RNA"
      s"${scRNAbamd}/$t/${a.barcode}.srt.rmdup.bam"
    }
  }

  val scSex2Barcode: Map[(String, String), List[String]] =
    scMeta
      .groupBy(x => (x.scTostr, x.sample.sex))
      .map((k, p) => (k, p.map(x => x.toRNABam)))

  val scRep2Barcode: Map[(String, String), List[String]] =
    scMeta
      .groupBy(x => (x.scTostr, x.sample.rep))
      .map((k, p) => (k, p.map(x => x.toRNABam)))

  val scG2Barcode =
    (scSex2Barcode.toSeq ++ scRep2Barcode.toSeq)
      .groupBy(_._1)
      .view.mapValues(_.map(_._2).toList.flatten)

  val scG2blstfnm =
    scG2Barcode.keys
      .map((sc, rep) => ((sc, rep), s"${blistd}/${sc}-${rep}.barcodes.txt"))
      .toMap

  val scG2outfnm =
    scG2Barcode.keys
      .map((sc, rep) => ((sc, rep), s"${outd}/${sc}-${rep}.bam"))
      .toMap

  class MergeRNABamTask(
    val group: String,
    val rep: String = "Male"
  ) extends TaskElement {
    val flagfnm = s"${myflagd}/${group}-${rep}.mergeBam.done"
    val logfnm = s"${logd}/${group}-${rep}.mergeRNABam.log"
    val skip = false
    def runCore(): Unit = {
      println(s"runCore merge: ${group}-${rep}")
      val blstfnm = scG2blstfnm((group, rep))
      val outfnm = scG2outfnm((group, rep))
      os.proc("samtools", "merge", "-f", "-o", outfnm, "-b", blstfnm)
        .call(check = false, stdout = os.Path(logfnm), stderr = os.Path(logfnm))
    }
  }

  def main(args: Array[String]) = {
    // 1. generate the barcode list fnm for each subclass group
    // scG2blstfnm.foreach { (k, blstfnm) =>
    //   {
    //     writeListOfString2File(
    //       content = scG2Barcode(k),
    //       to = blstfnm,
    //       overwrite = true
    //     )
    //   }
    // }
    // 2. start in parallel the mergebam script
    val tasks = scG2blstfnm.keys.map((sc, rep) => new MergeRNABamTask(sc, rep))
    tasks.par.foreach { t =>
      {
        t.run()
        println(s"done for ${t.group}-${t.rep}.")
      }
    }

    // val checkSorted = scG2outfnm.map(
    //   (k, f) => {
    //     val head = os.proc("samtools", "view", "-H", f)
    //          .call(check = false)
    //     (k, head)
    //   }
    // )

  }
}
