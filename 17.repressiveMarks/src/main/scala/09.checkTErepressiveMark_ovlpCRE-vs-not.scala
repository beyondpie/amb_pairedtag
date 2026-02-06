import os._
import bioscala.LightCoord.TE.loadHomerRepeats
import bioscala.LightCoord.GenomeCoord.{GenomeCoord, GenomeCoords}
import bioscala.LightCoord.Bed.{LightBed, LightBedElement}
import bioscala.LightCoord.TE.toTECategoryFromHomerString
import org.broad.igv.bbfile.BBFileReader
import bioscala.LightCoord.BedGraph.getBigWigReader
import MetaData.readSubclassMeta
import SZUtils.{ifelse, str2path, path2str, writeStrings2File}
import bioscala.LightCoord.BedGraph.{getWigValue, getbwOneRegion}
import bioscala.LightCoord.GenomeCoord.mkStringGenomeCoord
import bioscala.LightCoord.GenomeCoord.genomeCoordIgnoreStrandOrd
import MetaData.PairedTagCluster
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*
import scala.util.Random
import bioscala.lightype.{DMat, Strings}
import bioscala.LightCoord.Bed.mkStringLightBedElement
import scala.util.{Try, Success, Failure}

object RepressiveMarkOnTEOvlpCRE {
  // * meta data
  val projd    = os.Path("/Users/szu/git-recipes/amb_pairedtag")
  val workd    = projd / "17.repressiveMarks"
  val ptDNAbwd = projd / "data" / "ptDNAbw"
  val TEovlpCREfnm = workd / "out" / "ovlpTE" / "CRE-TE-TECat.tsv"
  val autoChrs: Vector[String] =
    1.to(19).map(i => s"chr$i").toVector
  val hs = Vector("H3K9me3", "H3K27me3", "H3K27ac", "H3K4me1")


  val TEClass = Vector("DNA", "LINE", "LTR", "SINE")

  val allTEfnm = os.Path("/Users/szu/softwares/homer") /
    "data" / "genomes" / "mm10" / "mm10.repeats"

  val ptscMeta = readSubclassMeta(
      s"$projd/meta/PairedTagSubclassMetaFromCellMetaChromHMM.csv")

  val nMin: Int = 10
  val ptscs: Strings =
    ptscMeta
      .filter(x => x.h3K9me3.female >= nMin && x.h3K9me3.male >= nMin)
      .filter(x => x.h3K27me3.female >= nMin && x.h3K27me3.male >= nMin)
      .filter(x => x.h3K27ac.female >= nMin && x.h3K27ac.male >= nMin)
      .filter(x => x.h3K4me1.female >= nMin && x.h3K4me1.male >= nMin)
      .map(x => x.name.PairedTagName)
  val bwd                = projd / "data" / "ptDNAbw"
  val emptyValue: Double = 0.0

  val CREAnnotfnm =
    projd / "data" / "CRE" / "allCRE.amb.PairedTag.annot.tsv"

  type AnnotCRE = (g: GenomeCoord, de: String, dist: String,
    fn: String, sc: String, cnsrv: String, a: String)

  val outd = workd / "out" / "ovlpTE"

  // * functions
  def loadSubclassesTrack(mark: String): Map[String, BBFileReader] = {
    ptscs
      .map(sc => {
        val suffix = ifelse(mark == "H3K9me3", "sm1000", "sm300")
        val f      = getBigWigReader(
            bwd / s"$sc.$mark.e100.bs100.$suffix.bw")
        (sc, f)
      })
      .toMap
  }

  def getbwMat(gs: GenomeCoords, scs: Strings, sc2bw: Map[String, BBFileReader]) = {
    gs.map(g => {
      scs.par
        .map(sc =>
          getbwOneRegion(sc2bw(sc), g, 0.0,
              weightByGivenRegion = true))
        .toVector
    })
  }

  def savebwMat(gs: GenomeCoords, scs: Strings, m: DMat, outd: os.Path): Unit = {
    if (!os.exists(outd)) {
      os.makeDir(outd)
    }
    writeStrings2File(
        content = m.map(x => x.mkString(sep = ",")),
        to = outd / "mat.csv",
        overwrite = true,
        head = ""
    )
    writeStrings2File(
        content = gs.map(g =>
          mkStringGenomeCoord(g, sep = "\t", useStrand = false)),
        to = outd / "te.bed",
        overwrite = false,
        head = ""
    )
    writeStrings2File(
        content = scs,
        to = outd / "scs.bed",
        overwrite = false,
        head = ""
    )
  }

  def loadDMat(d: os.Path): (GenomeCoords, Strings, DMat) = {
    val t: GenomeCoords = os.read.lines
      .stream(d / "te.bed")
      .map(x => x.strip().split("\t"))
      .filter(_.nonEmpty)
      .filter(_.length > 2)
      .map(x =>
        (chrom = x(0),
            coord = (startFrom = x(1).toInt, endTo = x(2).toInt),
            strand = "."))
      .toVector

    val s: Strings = os.read.lines
      .stream(d / "scs.bed")
      .map(_.strip())
      .filter(_.nonEmpty)
      .toVector

    val m: DMat = os.read.lines
      .stream(d / "mat.csv")
      .map(x => x.strip().split(","))
      .filter(_.nonEmpty)
      .filter(_.length == s.length)
      .map(x => x.map(i => i.toDouble).toVector)
      .toVector
    assert(t.size == m.size)
    (t, s, m)
  }

  def saveLightBed(x: LightBed, f: os.Path) = {
    val c = x.map(i => mkStringLightBedElement(i, sep = "\t"))
    writeStrings2File(content = c, to = f, overwrite = true,
        head = "")
  }

  // * main
  def main(args: Array[String]) = {

    // * load TE and seperate them into ovlpCRE and novlpCRE

    val CRE2TE: Vector[(GenomeCoord, GenomeCoord)] =
      os.read.lines
        .stream(TEovlpCREfnm)
        .map(x => x.strip().split("\t"))
        .filter(x => x.nonEmpty)
        .filter(x => x.length >= 2)
        .map(x =>
          (
              (chrom = x(0),
                  coord = (startFrom = x(1).toInt,
                      endTo = x(2).toInt), strand = "."),
              (chrom = x(3),
                  coord = (startFrom = x(4).toInt,
                      endTo = x(5).toInt), strand = ".")
          ))
        .toVector

    val allTEs: LightBed =
      loadHomerRepeats(allTEfnm)
        .filter(x => autoChrs.contains(x.g.chrom))
        .filter(x =>
          TEClass.contains(
              toTECategoryFromHomerString(x.name).TEclass))

    val TEovlpCRE: GenomeCoords  = CRE2TE.map(_._2).distinct
    val TEovlpCRESet             = TEovlpCRE.toSet
    val TEnovlpCRE: GenomeCoords =
      Random
        .shuffle(allTEs
              .filterNot(x => TEovlpCRESet.contains(x.g))
              .map(x => x.g))
        .take(TEovlpCRE.length)
        .sorted

    // * load subclass level repressive markers
    val sc2bwK9me3  = loadSubclassesTrack("H3K9me3")
    val sc2bwK27me3 = loadSubclassesTrack("H3K27me3")
    val sc2bwK27ac = loadSubclassesTrack("H3K27ac")
    val sc2bwK4me1 = loadSubclassesTrack("H3K4me1")

    // * get TE by subclass for each repressive markers
    hs.foreach(h => {
      val sc2bw = h match {
        case "H3K9me3" => sc2bwK9me3
        case "H3K27me3" => sc2bwK27me3
        case "H3K27ac" => sc2bwK27ac
        case _ => sc2bwK4me1
      }
      Vector("ovlpCRE", "nonOvlpCRE").foreach(o => {
        val outf = outd / s"TE${o}2sc${h}RPKM"
        if(os.exists(outf)) {
          println(s"$outf exist. skip.")
        } else {
          val tes = o match {
            case "ovlpCRE" => TEovlpCRE
            case _ => TEnovlpCRE
          }
          val m = getbwMat(tes, ptscs, sc2bw)
          savebwMat(tes, ptscs, m, outf)
          println(s"TE $o on $h signals are done.")
        }
      })
   })

    // * load subclass-level CREs, and CRE annotation
    val cres: Vector[AnnotCRE] =
      os.read.lines
      .stream(CREAnnotfnm)
      .slice(1, Int.MaxValue)
      .map(x => x.strip().split("\t"))
      .map(x =>
        (
            (
                g = (chrom = x(0),
                    coord = (x(1).toInt, x(2).toInt),
                    strand = "."),
                de = x(3),
                dist = x(4),
                fn = x(5),
                sc = x(6),
                cnsrv = x(7),
                a = x(8)
            )
        ))
      .toVector

    val te2cre: Map[GenomeCoord, GenomeCoord] =
      CRE2TE.map(x => (x._2, x._1)).toMap

    val sc2cre: Map[String, Set[GenomeCoord]] =
      cres.groupMap(_.sc)(_.g).map((k, v) => (k, v.distinct.toSet))

    val sc2creChrA: Map[String, Set[GenomeCoord]] =
      cres
        .filter(_.fn == "Chr-A")
        .groupMap(_.sc)(_.g)
        .map((k, v) => (k, v.distinct.toSet))

    val sc2creChrO: Map[String, Set[GenomeCoord]] =
      cres
        .filter(_.fn == "Chr-O")
        .groupMap(_.sc)(_.g)
        .map((k, v) => (k, v.distinct.toSet))

    // * load te x sc for each repressive markers
    val (teo, scs, omat9) =
      loadDMat(outd / "TEovlpCRE2scH3K9me3RPKM")

    val (_, _, omat27me3) = loadDMat(outd / "TEovlpCRE2scH3K27me3RPKM")

    val (_, _, omat27ac) = loadDMat(outd / "TEovlpCRE2scH3K27acRPKM")
    val (_, _, omat4) = loadDMat(outd / "TEovlpCRE2scH3K4me1RPKM")

    val (teno, _, nomat9) = loadDMat(outd / "TEnonOvlpCRE2scH3K9me3RPKM")

    val (_, _, nomat27me3) = loadDMat(outd / "TEnonOvlpCRE2scH3K27me3RPKM")

    val (_, _, nomat27ac) = loadDMat(outd / "TEnonOvlpCRE2scH3K27acRPKM")
    val (_, _, nomat4) = loadDMat(outd / "TEnonOvlpCRE2scH3K4me1RPKM")

    val teoRowIdMap: Map[GenomeCoord, Int] =
      teo.zipWithIndex.toMap
    val tenoRowIdMap: Map[GenomeCoord, Int] =
      teno.zipWithIndex.toMap
    val colIdMap: Map[String, Int] = scs.zipWithIndex.toMap

    // * create subclass-level TE signals
    val TE2name = allTEs.map(x => (
      x.g.toTuple.copy(_3 = ".").asInstanceOf[GenomeCoord], x.name)).toMap

    def organizeTE4sc(sc: String): Try[String] = Try {
      //hs.foreach(h => {
      Vector("H3K9me3").foreach(h => {
        val (omat, nomat) = h match {
          case "H3K9me3" => (omat9, nomat9)
          case "H3K27me3" => (omat27me3,nomat27me3)
          case "H3K4me1" => (omat4, nomat4)
          case _ => (omat27ac, nomat27ac)
        }
        val c     = cres.filter(x => x.sc == sc)

        // val teChrA: LightBed =
        //   teo
        //   .filter(t => {
        //     sc2creChrA(sc).contains(te2cre(t))
        //   })
        //   .map(t => {
        //     (
        //         g = t,
        //         name = TE2name(t),
        //         score = omat(teoRowIdMap(t))(colIdMap(sc))
        //     )
        //   })

        val teChrO: LightBed =
          teo
          .filter(t => {
            sc2creChrO(sc).contains(te2cre(t))
          })
          .map(t => {
            (
                g = t,
                name = TE2name(t),
                score = omat(teoRowIdMap(t))(colIdMap(sc))
            )
          })

        // val te: LightBed =
        //   teo
        //   .filter(t => {
        //     sc2cre(sc).contains(te2cre(t))
        //   })
        //   .map(t => {
        //     (
        //         g = t,
        //         name = TE2name(t),
        //         score = omat(teoRowIdMap(t))(colIdMap(sc))
        //     )
        //   })

        // // bg: backgrouand
        // val tebg: LightBed =
        //   Random
        //   .shuffle(teno)
        //   .take(te.length)
        //   .map(t => {
        //     (
        //         g = t,
        //         name = TE2name(t),
        //         score = nomat(tenoRowIdMap(t))(colIdMap(sc))
        //     )
        //   })

        // save all the te signals in a given path
        val o = outd / "sc2TE" / s"$sc.$h" 
        if(!os.exists(o)) {
          os.makeDir(o)
        }

        Map(
            //"teChrA" -> teChrA,
            "teChrO" -> teChrO,
            //"teCRE"  -> te,
            //"tebg"   -> tebg
        ).foreach((n, b) => {
          println(s"save LightBed for $sc - $h - $n.")
          if (b.nonEmpty) {
            saveLightBed(b, o / s"$n.tsv")
          } else {
            println("ignore since it's empty.")
          }
        }) // end of save data
      })   // end of histone foreach
      sc
    } // end of organize fn
    val r: Vector[Try[String]] = ptscs.map(sc => organizeTE4sc(sc))

    r.foreach(x => {
      x match {
        case Failure(e) => println(s"Failed: ${e.getMessage()}")
        case _ => 
      }
    })
  }        // end of main
}          // end of object
