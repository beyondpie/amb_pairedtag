import os._
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*
import GRange.{GenomicRange, genomeCoordStartOrd}
import ChromHMM.{ChromHMMState, PeakChromHMMStateAnnot}
import PairedTagChromHMM.PeakChromHMMStateSum
import SZUtils.{writeStrings2File, readTable}
import Bed.BedElement4
import MetaData.TSCCMeta
import ChromHMM.readAnnotFromDenseBed

// TODO: move GFF3 related functions to tools
object AnnotTSSbyChromHMM {

  /** Genes from Genetic Feature Fromat Version 3 (GFF3) format.
    *
    * Check object GFF3's doc for GFF3 format in details.
    */
  case class GFF3GeneElement(
    g: GenomicRange, geneId: String, geneType: String, geneName: String,
    strand: String) {
    def getTSS: (String, Int) = strand match {
      case "+" => (g.chrom, g.startFrom)
      case "-" => (g.chrom, g.endTo)
    }
  }

  /** Read Genetic Feature Fromat Version 3 (GFF3) format.
    *
    * See Wikipedia for deatails. GFF3, tab sep with 9 fields.
    *  - start is 1-base offset, included. While bed is 0-base offset,
    *    included.
    *  - end is 1-base offset, included while bed is not.
    */
  object GFF3{
    def readGFF3Gene(fnm: String): Vector[GFF3GeneElement] = {
       os.read.lines(os.Path(fnm))
        .filter(x => !x.matches("^#.*"))
        .map(x => x.split("\t"))
        .filter(x => x(2) == "gene")
        .toVector
        .map(y => {
          val d = y
            .last
            .split(";")
            .map(x => {
              val z = x.split("=")
              (z(0), z(1))
            })
            .toMap
          // organize the result
          new GFF3GeneElement(
            g = new GenomicRange(chrom = y(0), startFrom = y(3).toInt,
              endTo = y(4).toInt),
            geneId = d("gene_id"),
            geneType = d("gene_type"),
            geneName = d("gene_name"),
            strand = y(6)
          )
        })
    } // end of readGFF3Gene
  }

  def annotTSS(tss: Vector[(String, Int)], chromAnnot: Vector[BedElement4]): Vector[String] = {
    val c2a = chromAnnot.groupBy(_.x.chrom)
    tss.map(x => {
      val t = c2a(x._1).sortBy(_.x)(using genomeCoordStartOrd)
      t.find(p => (p.x.startFrom <= x._2) && (p.x.endTo >= x._2)).get.name
    })
  }

  def main(args: Array[String]) = {
    val projd = TSCCMeta.projd
    val densebedir = projd + "/06.ChromHMM/out/updateDenseBed"
    val autoChroms = 1.to(19).map(i => s"chr${i}")
    val outd = projd + "/06.ChromHMM/out/TSSAnnot"

    // 1. load TSS as Genomic Ranges
    val mouseGFF3fnm = "/projects/ps-renlab2/szu/genome/" +
      "gencode.vM25.annotation.gff3"

    // only genes on autosomes are considered.
    val genes = GFF3.readGFF3Gene(mouseGFF3fnm).filter(
      x => autoChroms.contains(x.g.chrom)
   )

    val tss = genes.map(_.getTSS)
    // 2. annotte TSS per subclass and save results
    // - load subclasses
    val ptscMeta = readTable(fnm = TSCCMeta.ptscMetafnm,
      sep=",", head = true)
    val scs: List[String] =
      ptscMeta
        .filter(x => x(10).toInt > 0)
        .map(x => x(0))
    // - for each subclass get the annotation.
    scs.par.foreach(sc => {
      val chromAnnot = readAnnotFromDenseBed(
        fnm = s"${densebedir}/${sc}_18_dense.bed")
      val tssAnnot = annotTSS(tss, chromAnnot)
      val linesOfTSSAnnot = tssAnnot.zip(genes).zip(tss).map(
        (x, z) => List(
          z._1,
          z._2,
          z._2 + 1,
          x._1,
          x._2.geneName,
          x._2.geneId,
          x._2.geneType
        ).mkString("\t")
      ).toList
      writeStrings2File(content = linesOfTSSAnnot,
        to = s"${outd}/${sc}.TSSAnnot.tsv",
        overwrite = true,
        head = List("chrom", "startFrom", "endTo", "ChromAnnot",
          "geneName", "geneId", "geneType").mkString("\t")
      )
      println(s"Annot TSSs for ${sc}: done.")
    }) // end of scs foreach
  } // end of main
}
