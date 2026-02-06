import os._
import scala.collection.parallel.CollectionConverters.*
import scala.collection.parallel._

import MetaData.TSCCMeta.{projd, metad, rep, modality, ptclusterDNAbamd}
import MetaData.{PairedTagBarcodeMeta, toStringOfAllenIdLabel, isNeu}
import MetaData.TSCCTools.samtools as samtools
import MetaData.PairedTagBarcodeMeta.readPairedTagCellMetaFile as readPairedTagMeta

def statReads(bamfnm: String): Int = {
  os.proc(samtools, "view -F 0x904 -c".split(" "),bamfnm)
  .call(check = false).out.text().strip.toInt
}

object StatMergedBam {
  def main(args: Array[String]) = {
    lazy val cellMeta: List[PairedTagBarcodeMeta] = readPairedTagMeta().
      filter(x => x.annotQuality == "Good")
    lazy val neuRep2NCell =
      cellMeta
        .filter(x => isNeu(x))
        .groupMapReduce(
          x => (x.annot.sc.get, x.exp.modularity, x.sample.rep))(x => 1)(
          (x, y) => x + y)
    lazy val neuRep2Bamfile =
      neuRep2NCell.keys
        .map(x => x ->
          List(ptclusterDNAbamd, "bam", toStringOfAllenIdLabel(x._1),
            s"${x._2}-${x._3}.srt.bam").mkString("/"))
        .toMap

    val forkJoinPool = new java.util.concurrent.ForkJoinPool(7)

    lazy val neuRep2nReadsPar = neuRep2Bamfile.par
    neuRep2nReadsPar.tasksupport = new ForkJoinTaskSupport(forkJoinPool)

    lazy val neuRep2nReads =
      neuRep2nReadsPar.map((k, v) => (k, statReads(v))).toMap

    lazy val neuRepStat =
      neuRep2NCell
        .map((k, v) => (k,  (v, neuRep2nReads(k)))).
        map((k,v) => (k, (v._1, v._2, v._2.toFloat / v._1.toFloat))).
        toList.sortBy((k, v) => v._1)
    val outfnm = os.Path(projd) / "04.peakcalling" / "src" / "main"/
    "resource" / "Neu.subclass.SexRep.stat.csv"
    val head = List("subclass", "modality",
      "sexrep", "ncell", "nreads", "ratio").mkString(",")
    os.write(outfnm, s"${head}\n")
    neuRepStat.map((k,v) => List(k._1, k._2, k._3, v._1, v._2, v._3)).
      map(x => x.mkString(",")).foreach {
      x =>  os.write.append(outfnm, s"${x}\n")
    }
  }
}





