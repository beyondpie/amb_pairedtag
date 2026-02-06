// Inspired by chromatin contact distances distribution from nn to neu.
import os._
import SZUtils.{str2path, path2str, ifelse, writeStrings2File}
import MetaData.{readSubclassMeta, TSCCMeta}
import Genome.MouseGenome
import bioscala.LightCoord.GenomeCoord.{GenomeCoord, GenomeCoords}
import project.Peak.loadPairedTagSubclassPeak

@main def CalRepressiveRegionRatio4Subclass(): Unit = {
  val projd = TSCCMeta.rawProjd
  val workd  = projd / "17.repressiveMarks"
  val outd   = workd / "out"
  val mods   = MetaData.TSCCMeta.modality
  val states =
    PairedTagChromHMM.chromHMMStates
    .filter(x => x != "ND")

  val mm10size: Double =
    MouseGenome.chr2size
    .filter((k, _) => k != "chrX" && k != "chrY")
    .values
    .map(x => x.toDouble)
    .sum

  val ptscMeta = readSubclassMeta().filter(x => x.atac)
  val ptscs    = ptscMeta.map(x => x.name.PairedTagName)

  val peakd = projd / "data" / "pairedtag_peak" / "subclass_peak"
  val chromHMMdensed = projd / "data" / "chromHMM" / "denseBed"

  def calGenomeRatio(x: GenomeCoords): Double = {
    x.map(x => x.coord.endTo - x.coord.startFrom)
      .sum
      .toDouble * 100 / mm10size
  }

  val sc2peakRatio: Vector[(String, Vector[Double])] =
    ptscs.map(sc => {
      val r =
        mods
          .map(mod => {
            val p = loadPairedTagSubclassPeak(sc, mod)
            calGenomeRatio(p)
          })
          .toVector
      (sc, r)
    })
  writeStrings2File(
    content = sc2peakRatio.map((sc, r) =>
      sc + "," + r.mkString(",")),
    to = outd / "pt.subclassPeakGenomeRaioCapture.prcnt.csv",
    overwrite = true, head = "subclass," + mods.mkString(","))

  // * check ChromHMM Chr-P and Chr-H
  def loadscChromHMM(sc: String): Map[String, GenomeCoords] = {
    val fnm: os.Path = chromHMMdensed / s"${sc}_18_dense.bed"
    os.read.lines
      .stream(fnm)
      .drop(1)
      .map(x => x.strip().split("\t"))
      .filter(x => x.nonEmpty)
      .filter(x => x.length >= 4)
      .map(x => {
        val g: GenomeCoord = (chrom = x.head,
            coord = (x(1).toInt, x(2).toInt), strand = ".")
        val s = x(3)
        (s, g)
      })
      .toVector
      .groupMap(_._1)(_._2)
  }

  // val sc2ChromHMMState =
  //   ptscs
  //     .map(sc => {
  //       (sc, loadscChromHMM(sc))
  //     })
  //     .toMap

  val sc2ChromHMMRatio: Vector[(String, Vector[Double])] =
    ptscs.map(sc => {
      val t = loadscChromHMM(sc)
      val r = states
        .map(state => {
          if (t.contains(state)) {
            calGenomeRatio(t(state))
          } else {
            0.0
          }
        })
        .toVector
      (sc, r)
    })

  writeStrings2File(content = sc2ChromHMMRatio.map((sc, r) =>
        sc + "," + r.mkString(",")),
      to = outd / "pt.subclassChromHMMGenomeRaioCapture.prcnt.csv",
      overwrite = true, head = "subclass," + states.mkString(","))

} // end of main
