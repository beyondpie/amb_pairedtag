import os._
import bioscala.LightCoord.Bed.LightBed
import bioscala.lightype.{DMat, Strings}
import bioscala.LightCoord.GenomeCoord.{GenomeCoord, GenomeCoords}
import bioscala.LightCoord.TE.toTECategoryFromHomerString
import bioscala.LightCoord.TE.TECategory
import RepressiveMarkOnTEOvlpCRE.{savebwMat, loadDMat}
import SZUtils.{str2path, path2str}
import SZUtils.writeStrings2File
import bioscala.LightCoord.TE.loadHomerRepeats


object Subclass2TEChrA {
  // * meta
  val projd = os.Path("/Users/szu/git-recipes/amb_pairedtag")
  val workd    = projd / "17.repressiveMarks"
  val outd = workd / "out" / "ovlpTE"
  val sc2TEd = workd / "out" / "ovlpTE" / "sc2TE"

  val ptscs: Strings =
    os.list(sc2TEd)
    .map(_.baseName).distinct
      .filter(x => x.length() > 1).toVector
  val hs = Vector("H3K9me3", "H3K27me3", "H3K27ac", "H3K4me1")

  val nan: Double = -1.0
  val g: String = "subfamily"

  val allTEfnm = os.Path("/Users/szu/softwares/homer") /
    "data" / "genomes" / "mm10" / "mm10.repeats"

  val TEClass = Vector("DNA", "LINE", "LTR", "SINE")
  val autoChrs: Vector[String] =
    1.to(19).map(i => s"chr$i").toVector

  val allTEs: LightBed =
    loadHomerRepeats(allTEfnm)
      .filter(x => autoChrs.contains(x.g.chrom))
      .filter(x =>
        TEClass.contains(
          toTECategoryFromHomerString(x.name).TEclass))

  val allTE2Subfam: Map[GenomeCoord, String] =
    allTEs
      .filter(x => !x.name.contains("Unknown"))
      .filter(x => !x.name.contains("Other"))
      .filter(x => !x.name.contains("?"))
      .map(x => (x, x.name.strip().split("\\|")))
      .filter(_._2.length == 3)
      .map((t, n) => {
        val r = "-HOMER\\d+$".r
        val fam = r.replaceAllIn(n.last, "")
        val sf = n.head
        (
          t.g.toTuple.copy(_3 = ".").asInstanceOf[GenomeCoord],
          s"$fam:$sf"
        )
      }).toMap

  // * functions

  def loadsc2TE(sc: String, h: String): LightBed = {
    //val f = sc2TEd / s"$sc.$h" / "teChrA.tsv"
    val f = sc2TEd / s"$sc.$h" / "teChrO.tsv"
    os.read.lines
      .stream(f)
      .map(x => x.strip().split("\t"))
      .filter(_.length > 1)
      .map(x => (
        g = (
          chrom = x(0),
          coord = (startFrom = x(1).toInt, endTo = x(2).toInt),
          strand = "."
        ),
        name = x(3),
        score = x(4).toDouble
      ))
      .toVector
  }

  def savesc2TEgroup(scs:Strings, g: Strings, r: DMat, outd: String): Unit = {
    if (!os.exists(outd)) {
      os.makeDir(outd)
    }
    writeStrings2File(
      content = scs,
      to = outd / "scs.txt",
      overwrite = true,
      head = ""
    )
    writeStrings2File(
      content = g,
      to = outd / "TEgroup.txt",
      overwrite = true,
      head = ""
    )
    writeStrings2File(
      content = r.map(x => x.mkString(sep = ",")),
      to = outd / "mat.csv",
      overwrite = true,
      head = ""
    )

  }

  // * main
  def main(args: Array[String]) = {
    hs.foreach(h => {
      val (tes, scs, teOvlpCRE2scMat) =
        loadDMat(outd / s"TEovlpCRE2sc${h}RPKM")
      val TE2RowId = tes.zipWithIndex.toMap
      val sc2ColId = scs.zipWithIndex.toMap
      // 249099
      val TECRE2Subfam: Map[GenomeCoord, String] =
        tes.filter(t => allTE2Subfam.contains(t))
          .map(t => (t, allTE2Subfam(t))).toMap

      val subFam2TECRE: Map[String, Iterable[GenomeCoord]] =
        TECRE2Subfam
          .toVector
          .groupMap(_._2)(_._1)

      val sc2TEChrASignal: Map[String, LightBed] = ptscs.map(sc => {
        (sc, loadsc2TE(sc, h))
      }).toMap

      // length of 123827
      val allTEChrAs: GenomeCoords =
        ptscs.map(sc => {
          val t = sc2TEChrASignal(sc)
          t.map(_.g)
        }).flatten
          .distinct
          .filter(g => TECRE2Subfam.contains(g))


      // val sc2TEChrAMat: DMat =
      //   ptscs.map(sc => {
      //     val colId = sc2ColId(sc)
      //     allTEChrAs.map(x => {
      //       teOvlpCRE2scMat(TE2RowId(x))(colId)
      //     })
      //   })

      // savebwMat(allTEChrAs, ptscs, sc2TEChrAMat, outd / s"ptsc2TEChrAMat$h")
      // println(s"finish sc2TEChrA matrix for $h.")

      val TEChrA2Cat: Map[GenomeCoord, String] =
        allTEChrAs.map(x => (x, TECRE2Subfam(x))).toMap

      // 990
      val subFam2TEChrA: Map[String, Iterable[GenomeCoord]] =
        TEChrA2Cat.map((t, c) => {
        (c, t)
      }).groupMap(_._1)(_._2)
        
      val subfams: Strings = subFam2TEChrA.keys.toVector

      // val sc2sfWithNaN: DMat = ptscs.map(sc => {
      //   subfams.map(sf => {
      //     val s = sc2TEChrASignal(sc).filter(t => {
      //       TEChrA2Cat(t.g).TEsubfamily == sf
      //     })
      //     if (s.length < 1) {
      //       nan
      //     } else {
      //       s.map(t => t.score).sum / s.length.toDouble
      //     }
      //   }) // end of subfams map
      // }) // end of ptscs map

      // val onan: String = outd / s"sc2TE${g}_${h}_withNaN"
      // savesc2TEgroup(ptscs, subfams, sc2sfWithNaN, onan)
      // println(s"finish sc2TE groupby $g for $h with NaN.")

      // val sc2sfWithoutNaN: DMat = ptscs.map(sc => {
      //   subfams.map(sf => {
      //     val s = sc2TEChrASignal(sc).filter(t => {
      //       TEChrA2Cat(t.g).TEsubfamily == sf
      //     })
      //     if (s.length < 1) {
      //       val tmp = subFam2TECRE(sf)
      //         .map(t => teOvlpCRE2scMat(TE2RowId(t))(sc2ColId(sc)))
      //         .toVector
      //       tmp.sum / tmp.length.toDouble
      //     } else {
      //       s.map(t => t.score).sum / s.length.toDouble
      //     }
      //   })
      // })
      // val oWithoutNaN: String = outd / s"sc2TE${g}_${h}_withoutNaN"
      // savesc2TEgroup(ptscs, subfams, sc2sfWithoutNaN, oWithoutNaN)
      // println(s"finish sc2TE groupby $g for $h without NaN.")

      // use all the TEs (ovlp with CRE) in a subfamily
      val sc2sfAll: DMat = ptscs.map(sc => {
        subfams.map(sf => {
          val t = subFam2TECRE(sf)
            .map(t => teOvlpCRE2scMat(TE2RowId(t))(sc2ColId(sc)))
            .toVector
          t.sum / t.length.toDouble
        })
      })
      // val oAll = outd / s"sc2TE${g}_${h}_All"
      val oAll = outd / s"sc2TE${g}_${h}_All_ChrO"
      savesc2TEgroup(ptscs, subfams, sc2sfAll, oAll)

    }) // end of hs foreach
  } // end of main
} // end of object
