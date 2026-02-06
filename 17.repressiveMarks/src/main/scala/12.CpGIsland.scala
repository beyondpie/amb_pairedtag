import os._
import bioscala.LightCoord.GenomeCoord.GenomeCoords
import org.broad.igv.bbfile.BBFileHeader
import MetaData.readSubclassMeta
import SZUtils.{ifelse, str2path, path2str, writeStrings2File}
import bioscala.LightCoord.BedGraph.{getWigValue, getbwOneRegion}
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*
import bioscala.lightype.Strings
import bioscala.LightCoord.BedGraph.getBigWigReader
import org.broad.igv.bbfile.BBFileReader
import bioscala.LightCoord.GenomeCoord.mkStringGenomeCoord

object getCpGIslandRepressiveSignal {
  val projd = os.Path("/Users/szu/git-recipes/amb_pairedtag")
  val workd = projd / "17.repressiveMarks"
  val outd = workd / "out"
  val CpGIslandfnm = projd / "meta" / "mm10.cpgIslandExt.bed"
  val mCGbwd = projd / "data" / "mCGbw"
  val ptDNAbwd = projd / "data" / "ptDNAbw"

  val autoChrs: Vector[String] =
    1.to(19).map(i => s"chr$i").toVector

  val ptscMeta = readSubclassMeta(
    s"$projd/meta/PairedTagSubclassMetaFromCellMetaChromHMM.csv")
  val nMin: Int = 10
  val ptscs: Strings =
    ptscMeta.filter(x => x.h3K27me3.female >= nMin && x.h3K27me3.male >= nMin)
      .filter(x => x.meDNA)
      .map(x => x.name.PairedTagName)

  val sc2mCGbw: Map[String, BBFileReader] =
    ptscs
      .filter(sc => os.exists(mCGbwd / s"$sc.CGN.bw"))
      .map(sc =>
        (sc, getBigWigReader(mCGbwd / s"$sc.CGN.bw")))
      .toMap

  val sc2K27me3bw: Map[String, BBFileReader] =
    ptscs
      .filter(sc => sc2mCGbw.contains(sc))
      .map(sc => (sc, getBigWigReader(
        ptDNAbwd / s"$sc.H3K27me3.e100.bs100.sm300.bw"))
      )
      .toMap

  val CpGs: GenomeCoords =
    os.read
      .lines.stream(CpGIslandfnm)
      .map(x => x.strip().split("\t"))
      .filter(x => x.length >= 4)
      .map(x =>
        (chrom = x(0),
          coord = (startFrom = x(1).toInt, endTo = x(2).toInt),
          strand = "."
        ))
      .filter(x => autoChrs.contains(x.chrom))
      .toVector

  def saveFilterCpG(gs: GenomeCoords, outf: String): Unit = {
    writeStrings2File(
      content = gs.map(x => mkStringGenomeCoord(x, sep="\t", useStrand = false)),
      to = outf,
      overwrite = true,
      head = ""
    )
  }

  def getCpGmCGK27me3(scs: Strings, CpGs: GenomeCoords): Vector[(String, Double, Double)] = {
    scs.par.map(sc => {
      CpGs.map(r =>
        (sc,
          getbwOneRegion(sc2mCGbw(sc), r, 0.0, false),
          getbwOneRegion(sc2K27me3bw(sc), r, 0.0, true)
        ))
    })
      .toVector
      .flatten
  }

  def saveCpGmCGK27me3(r: Vector[(String, Double, Double)], out: String): Unit = {
    writeStrings2File(
      content = r.map(x => s"${x._1},${x._2.toString()},${x._3.toString()}"),
      to = out,
      overwrite = true,
      head = "subclass,mCG,H3K27me3"
    )
  }

  saveFilterCpG(CpGs, s"$outd/mm10autoChrCpG.bed")
  val scs = ptscs.filter(sc => sc2K27me3bw.contains(sc))
  val r = getCpGmCGK27me3(scs, CpGs)
  saveCpGmCGK27me3(r, s"$outd/ptscsCpG.mCG.K27me3.csv")

}
