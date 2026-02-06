import os._
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*
import GRange.{GenomicRange, genomeCoordStartOrd, mouseGenomicRangeOrd}
import ChromHMM.{ChromHMMState, PeakChromHMMStateAnnot}
import PairedTagChromHMM.PeakChromHMMStateSum
import SZUtils.{writeStrings2File, readTable}
import Bed.BedElement4
import MetaData.TSCCMeta
import ChromHMM.readAnnotFromDenseBed
import CEMBATAC.loadscATAC 


object AnnotCREbySummitUsingChromHMM {

  // TODO: add this to tool
  def annotSummit(s: Vector[(String, Int)],
    a: Vector[BedElement4]): Vector[String] = {
    val c2a = a.groupBy(_.x.chrom)
    s.map(x => {
      val t = c2a(x._1).sortBy(_.x)(using genomeCoordStartOrd)
      t.find(p => (p.x.startFrom <= x._2) && (p.x.endTo >= x._2)).get.name
    })
  }

  def main(args: Array[String]) = {
    val projd = TSCCMeta.projd
    val densebedir = projd + "/06.ChromHMM/out/updateDenseBed"
    val autoChroms = 1.to(19).map(i => s"chr${i}")
    val outd = projd + "/06.ChromHMM/out/CREAnnot250429"
    if(!os.exists(os.Path(outd))) {
      os.makeDir(os.Path(outd))
    }

    val ptscMeta = readTable(fnm = TSCCMeta.ptscMetafnm, sep = ",", head = true)
    val scs: List[String] =
      ptscMeta
        .filter(x => x(10).toInt > 0)
        .map(x => x(0))

    scs.par.foreach(sc => {
      val chromAnnot = readAnnotFromDenseBed(
        fnm = s"${densebedir}/${sc}_18_dense.bed")
      val CREs: Vector[GenomicRange] = loadscATAC(sc).map(_.x).toVector
      val s: Vector[(String, Int)] = CREs.map(
        // in scala ( Int + Int ) / 2 will be floor (Int again)
        x => (x.chrom, (x.startFrom + x.endTo) / 2 ))
      val sa = annotSummit(s, chromAnnot)
      val linesOfCREAnnot = CREs.zip(sa).map((x, a) => {
        x.mkString("\t") + "\t" + a
      })
      writeStrings2File(content = linesOfCREAnnot,
        to = s"$outd/$sc.CRE.annotBySummit.bed",
        overwrite = true,
        head = ""
      )
      CREs.zip(sa)
        .groupBy((x, a) => a)
        .foreach{(annot, y) => {
          val x = y.map(_._1)
          val o = x.sorted(using mouseGenomicRangeOrd)
            .map(x => x.mkString(sep = "\t"))
          writeStrings2File(
            content = o,
            to = s"$outd/$sc.CRE.$annot.annotBySummit.bed",
            overwrite = true,
            head = ""
          )
        }}
      println(s"Annot CREs for ${sc}: done.")
    }) // end of scs foreach
  }
}
