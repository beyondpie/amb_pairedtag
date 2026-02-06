import os._
import bioscala.LightCoord.GenomeCoord.{GenomeCoord, GenomeCoords}
import bioscala.LightCoord.getCenter
import bioscala.LightCoord.coordOrd
import Genome.MouseGenome
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters._
import SZUtils.{str2path, path2str, ifelse, writeStrings2File}

type CREAnnot = (g: GenomeCoord, isDAR: String, isDistal: String, state: String, sc: String)

@main def calChrAChrOdist(): Unit = {
  // * meta
  val projd = "/projects/ps-renlab2/szu/projects/amb_pairedtag"
  val ptscMetafnm = os.Path(projd)/"meta"/
    "PairedTagSubclassMetaFromCellMetaChromHMM.csv"
  val workd = os.Path(projd) / "06.ChromHMM"
  val outd = workd / "out" / "distOfChrAChrO"
  if(!os.exists(outd)){
    os.makeDir(outd)
  }
  val ptscMeta = MetaData.readSubclassMeta(f = ptscMetafnm).filter(_.atac)
  val ptscs = ptscMeta.map(x => x.name.PairedTagName)
  val chromSize: Map[String, Int] = MouseGenome.chr2size
  val allCREfnm = os.Path(projd)/"data"/"chromHMM"/
    "allCRE.amb.PairedTag.annot.tsv"


  // * param
  val binSize: Int = 1E6.toInt
  // val state: String = "Chr-O"
  val useDiff: Boolean = true

  // * functions

  def readAllCRE: Vector[CREAnnot] = {
    os.read.lines.stream(allCREfnm)
    .drop(1)
    .map(x => x.strip().split("\t"))
    .filter(x => x.nonEmpty)
    .map(x => (g = (chrom = x.head, coord = (x(1).toInt, x(2).toInt), strand = "."),
      isDAR = x(3), isDistal = x(4), state = x(5), sc = x.last
    )).toVector
  }

  def distCRE(x: GenomeCoord, y: GenomeCoord): Double = {
    val a = getCenter(x.coord)
    val b = getCenter(y.coord)
    math.abs((a.startFrom - b.startFrom).toDouble)
  }

  def avgDist(x: GenomeCoords): Double = {
    if(x.length < 2) {
      Double.NaN
    } else {
      val y = x.combinations(2).toVector
      y.map(z => distCRE(z.head, z.last)).toVector.sum / y.length.toDouble
    }
  }

  /* Returns genomeBin Id for a given genome coord.
   * GenomeBin Id starts from zero. If a genomecoord overlapps with two bins,
   * left one is used.
   *     A          B
   *    ___       ----
   * |______|______|___| (a chromsome)
   * ^      ^          ^
   * |      |          | (end of chromsome, not used for binId.)
   * |      BinId of B
   * BinId of A
   *
   * @param x Genome coord
   * @param binSize genomeBin size
   */
  def getGenomeBinId(x: GenomeCoord, binSize: Int): Int = {
    x.coord.startFrom / binSize
  }

  def loadCRE(allCRE: Vector[CREAnnot], sc: String, state: String, useDiff: Boolean): GenomeCoords = {
    allCRE
      .filter(x => x.sc == sc)
      .filter(x => x.state == state)
      .filter(x => ifelse(useDiff, x.isDAR == "DAR", true))
      .map(x => x.g)
  }
  def getoutfnm(sc: String, state: String, useDiff: Boolean): String = {
    if (useDiff) {
      outd / s"$sc.$state.DAR.dist.csv"
    }
    else {
      outd / s"$sc.$state.dist.csv"
    }
  }
  
  def saveDist(x: Map[String, Vector[Double]], outf: String): Unit = {
    val r = 1.to(19)
      .map(i => s"chr$i")
      .filter(chr => x.contains(chr))
      .flatMap(chr => {
        x(chr).zipWithIndex.map((y, i) => s"$chr,$i,${y.toString}")
      })
    writeStrings2File(content = r, to = outf,
      head = "chrom,binId,avgDist", overwrite = true)
  }

  // * main
  val allCRE = readAllCRE

  Vector("Chr-A", "Chr-O").foreach(state => {
    ptscs.par.foreach(sc => {
      println(s"Start: $sc $state dist with DE $useDiff for binsize of $binSize.")
      val gs: GenomeCoords = loadCRE(allCRE, sc, state, useDiff)
      val r =
        gs
          .groupBy(_.chrom)
          .map((chrom, xs) => {
            val s = xs.sortBy(_.coord)(using coordOrd)
              .map(x => (getGenomeBinId(x, binSize), x))
              .groupBy(_._1)
              .filter((binId, ys) => ys.size > 1)
              .map((binId, ys) => avgDist(ys.map(x => x._2)))
              .toVector
              .filter(x => !x.isNaN)
            (chrom, s)
          })
      saveDist(r, getoutfnm(sc, state, useDiff))
      println(s"End: $sc $state dist with DE $useDiff for binsize of $binSize.")
    })
  })
} // end of main
