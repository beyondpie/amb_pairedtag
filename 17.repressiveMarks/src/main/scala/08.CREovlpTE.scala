import os._
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*
import SZUtils.{str2path, path2str, ifelse}
import SZUtils.writeStrings2File
import CEMBATAC.loadAllCRE
import bioscala.LightCoord.GenomeCoord.{GenomeCoord, GenomeCoords}
import bioscala.LightCoord.GenomeCoord.mkStringGenomeCoord
import bioscala.LightCoord.Bed.{LightBed, LightBedElement}
import bioscala.LightCoord.GenomeCoord.findSubjectOvlpGenomeCoordIgnoreStrand
import bioscala.LightCoord.TE.{TECategory, toTECategoryFromHomerString}
import bioscala.LightCoord.TE.mkStringTECat


object CREOvlpTE {

  // * meta data
  val projd = os.Path("/Users/szu/git-recipes/amb_pairedtag")
  val workd: Path = projd / "17.repressiveMarks"
  val outd: Path = workd / "out" / "ovlpTE"
  val wmbCREfnm = projd /
    "data" / "snATAC" / "wmb.snATAC.peak.srt.bed"
  val repeatfnm = os.Path("/Users/szu/softwares/homer") /
    "data" / "genomes" / "mm10" / "mm10.repeats"
  val autoChrs: Vector[String] = 1.to(19).map(i => s"chr$i").toVector
  val TEClass = Vector(
    "DNA",
    "LINE",
    "LTR",
    "SINE"
  )

  // * functions
  def loadRepeats(f: os.Path): LightBed = {
    os.read.lines.stream(f)
      .map(x => x.strip().split("\t"))
      .filter(x => x.nonEmpty)
      .filter(x => x.length >= 5)
      .map(x => (
        g = (
          chrom = x(1),
          coord = (startFrom = x(2).toInt, endTo = x(3).toInt),
          strand = ifelse(x(4).toInt > 0, "-", "+")
        ), name = x(0), score = x(4).toDouble))
      .filter(x => autoChrs.contains(x.g.chrom))
      .filter(x => TEClass.contains(toTECategoryFromHomerString(x.name).TEclass))
      .toVector
  }


  def main(args: Array[String]): Unit = {
    // 1. load all CREs
    val CREs = loadAllCRE(wmbCREfnm)
    val chr2CREs: Map[String, GenomeCoords] = CREs.groupBy(_.chrom)

    // 2. load all TEs
    val TEs = loadRepeats(repeatfnm)
    val TEg = TEs.map(_.g)
    val chr2TEg: Map[String, GenomeCoords] = TEg.groupBy(_.chrom)

    val TE2Cat: Map[GenomeCoord, TECategory] = TEs.map(
      x => (x.g, toTECategoryFromHomerString(x.name))
    ).toMap

    // 3. map CRE on TEs
    val CRE2TE = autoChrs.filter(i => chr2TEg.keySet.contains(i)).par.map(i => {
      val c = chr2CREs(i)
      val t = findSubjectOvlpGenomeCoordIgnoreStrand(
        query = c,
        subject = chr2TEg(i)
      )
      t.zipWithIndex.filter(
        (x, i) => x.coord.endTo != 0
      ).filter(
        (x, i) => TE2Cat.contains(x))
        .map((x, i) => (c(i), x, TE2Cat(x)))
    }).toVector.flatten

    // val TEsOvlpWithCRE = findSubjectOvlpGenomeCoordIgnoreStrand(
    //   query = CREs,
    //   subject = TEg
    // )

    // val CRE2TE = TEsOvlpWithCRE.zipWithIndex.filter(
    //   (x, i) => x.coord.endTo != 0
    // ).map((x, i) => (CREs(i), x, TE2Cat(x)))

    val CRE2TEstr = CRE2TE.map(
      (cre, te, tecat) => {
        val crestr = mkStringGenomeCoord(cre, sep = "\t", useStrand = false)
        val testr = mkStringGenomeCoord(te, sep = "\t", useStrand = false)
        val tecatstr = mkStringTECat(tecat, sep = "\t")
        s"$crestr\t$testr\t$tecatstr"
      }
    )

    writeStrings2File(
      content = CRE2TEstr, to = s"$outd/CRE-TE-TECat.tsv",
      overwrite = true,
      head = ""
    )


    // 4. add CREs' annotation
    //    - ChrA /ChrO
    //    - subclass-level
  }
}
