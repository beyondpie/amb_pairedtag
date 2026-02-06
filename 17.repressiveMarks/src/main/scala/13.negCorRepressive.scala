import io._
import MetaData.readSubclassMeta
import bioscala.LightCoord.GenomeCoord.GenomeCoord
import bioscala.LightCoord.GenomeCoord.GenomeCoords
import bioscala.LightCoord.GenCode.getGFF3Genes
import bioscala.lightype.Strings
import bioscala.LightCoord.BedGraph.getBigWigReader
import org.broad.igv.bbfile.BBFileReader
import SZUtils.{ifelse, str2path, path2str, writeStrings2File}
import bioscala.LightCoord.BedGraph.getbwOneRegion
import bioscala.LightCoord.Bed.LightBed
import bioscala.LightCoord.Bed.LightBedElement
import bioscala.LightCoord.Coord
import bioscala.LightCoord.GenCode.GFF3
import bioscala.LightCoord.GenCode.readGeneFromGenCodeGFF
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*

object NegCorRepGene {

  val projd      = os.Path("/Users/szu/git-recipes/amb_pairedtag")
  val workd      = projd / "17.repressiveMarks"
  val outd       = workd / "out" / "negCorRepGene"
  val ptDNAbwd   = projd / "data" / "ptDNAbw"
  val ptRNAbwd   = projd / "data" / "ptRNAbw"
  val gencodefnm = projd / "meta" / "gencode.vM25.annotation.gff3"

  val autoChrs: Vector[String] =
    1.to(19).map(i => s"chr$i").toVector

  val ptscMeta = readSubclassMeta(
      s"$projd/meta/PairedTagSubclassMetaFromCellMetaChromHMM.csv")
  val nMin: Int      = 50
  // 136 subclasses
  val ptscs: Strings =
    ptscMeta
      .filter(x =>
        x.h3K27me3.female >= nMin && x.h3K27me3.male >= nMin)
      .filter(x =>
        x.h3K9me3.female >= nMin && x.h3K9me3.male >= nMin)
      .map(x => x.name.PairedTagName)


  def getGeneBody(x: LightBedElement, upstream: Int = 2000): LightBedElement = {
    val g = if (x.g.strand == "-") {
      x.g.toTuple
        .copy(_2 = x.g.coord.toTuple
          .copy(_2 = x.g.coord.endTo + upstream)
          .asInstanceOf[Coord])
        .asInstanceOf[GenomeCoord]
    } else {
      x.g.toTuple
        .copy(_2 = x.g.coord.toTuple
          .copy(_1 = (x.g.coord.startFrom - upstream).max(0))
          .asInstanceOf[Coord])
        .asInstanceOf[GenomeCoord]
    }
    x.toTuple.copy(_1 = g).asInstanceOf[LightBedElement]
  }

  def main(args: Array[String]) = {

    val gff3: GFF3 =
      readGeneFromGenCodeGFF(gencodefnm)
        .filter(x => x.feaType == "gene")
        .filter(x => x.g.chrom != "chrM")
        .filter(x => x.g.chrom != "chrX" && x.g.chrom != "chrY")
        .filter(x => x.bioType == "protein_coding")

    val genes = gff3.map(x => {
      (
          g = x.g,
          name = x.gname,
          score = 0.0
      )
    })

    val geneBodys = genes.map(g => getGeneBody(g))

    val sc2K27me3bw: Map[String, BBFileReader] =
      ptscs
        .map(sc =>
          (sc,
              getBigWigReader(
                  ptDNAbwd / s"$sc.H3K27me3.e100.bs100.sm300.bw")))
        .toMap

    val sc2K9me3bw: Map[String, BBFileReader] =
      ptscs
        .map(sc =>
          (sc,
              getBigWigReader(
                  ptDNAbwd / s"$sc.H3K9me3.e100.bs100.sm1000.bw")))
        .toMap

    val sc2RNAbw: Map[String, BBFileReader] =
      ptscs.map(sc =>
          (
            sc,
            getBigWigReader(ptRNAbwd / s"$sc.RPKM.bw")
          ))
        .toMap

    val g2sc: Vector[Vector[Double]] =
      genes
        .map(x => {
          ptscs.par.map(sc => {
              getbwOneRegion(sc2RNAbw(sc), x.g, 0.0, true)
            })
            .toVector
        })
        .toVector

    val K9me3toSc: Vector[Vector[Double]] =
      geneBodys
        .map(x => {
          ptscs.par.map(sc => {
            getbwOneRegion(sc2K9me3bw(sc), x.g, 0.0, true)
          })
            .toVector
        })
        .toVector

    val K27me3toSc: Vector[Vector[Double]] =
      geneBodys
        .map(x => {
          ptscs.par.map(sc => {
            getbwOneRegion(sc2K27me3bw(sc), x.g, 0.0, true)
          })
            .toVector
        })
        .toVector

    // save results
    writeStrings2File(
        content = genes.map(_.name),
        to = outd / "genes.txt",
        overwrite = true,
        head = ""
    )

    writeStrings2File(
        content = ptscs,
        to = outd / "ptscs.txt",
        overwrite = true,
        head = ""
    )

    writeStrings2File(
        content = g2sc.map(x => x.mkString(sep = ",")),
        to = outd / "gene2ptscs_RNA.csv",
        overwrite = true,
        head = ""
    )

    writeStrings2File(
        content = K9me3toSc.map(x => x.mkString(sep = ",")),
        to = outd / "gene2ptscs_H3K9me3.csv",
        overwrite = true,
        head = ""
    )

    writeStrings2File(
        content = K27me3toSc.map(x => x.mkString(sep = ",")),
        to = outd / "gene2ptscs_H3K27me3.csv",
        overwrite = true,
        head = ""
    )

  } // end of main
}
