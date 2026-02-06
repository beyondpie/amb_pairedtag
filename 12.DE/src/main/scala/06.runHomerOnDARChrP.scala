import os.*
import MetaData.readSubclassMeta
import GRange.{GenomicRange, filterByRanges}
import SZUtils.writeStrings2File
import Homer2.MotifFinderByHomer2
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*
import ChromHMM.loadDenseBed

@main def HomerOnDARChrP(): Unit = {
  val projd = MetaData.TSCCMeta.projd
  val workd = projd + "/12.DE"
  val flagd = workd + "/flag"
  val logd = workd + "/log"
  val densebedir = projd + "/06.ChromHMM/out/updateDenseBed"

  val DARH3K4me1d = projd + "/12.DE/out/DARH3K4me1"
  val DARChrPd = projd + "/12.DE/out/DARChrP"
  if (!os.exists(os.Path(DARChrPd))) {
    os.makeDir(os.Path(DARChrPd))
  }

  val outHomerd = projd + "/12.DE/out/DARChrPHomer"
  if (!os.exists(os.Path(outHomerd))) {
    os.makeDir(os.Path(outHomerd))
  }

  val ptscMeta = readSubclassMeta().filter(_.atac)
  val ptscs = ptscMeta.map(_.name.PairedTagName)
  // val ptscs = MetaData.representPaiedTagSubclass

  // 1. get ChrP overlapping with DARH3K4me1
  val sc2ChrP: Map[String, Vector[GenomicRange]] = ptscs.map(sc => {
    val densebedfnm = densebedir + s"/${sc}_18_dense.bed"
    val g = os.read.lines.stream(os.Path(densebedfnm))
      .drop(1)
      .map(x => x.strip().split("\t"))
      .filter(x => x(3) == "Chr-P")
      .map(x => new GenomicRange(x.head, x(1).toInt, x(2).toInt))
      .toVector
    (sc, g)
  }).toMap

  val sc2DARH3K4me1: Map[String, Vector[GenomicRange]] = ptscs.map(sc => {
    val bedfnm = DARH3K4me1d + s"/$sc.DAR.H3K4me1Peak.bed"
    val g = os.read.lines.stream(os.Path(bedfnm))
      .map(x => x.strip().split("\t"))
      .map(x => new GenomicRange(x.head, x(1).toInt, x.last.toInt))
      .toVector
    (sc, g)
  }).toMap

  val sc2DARChrP: Map[String, Vector[GenomicRange]] = sc2ChrP.map((sc, gs) => {
    val r = filterByRanges(gs, sc2DARH3K4me1(sc)).toVector
    (sc, r)
  })

  sc2DARChrP.foreach((sc, gs) => {
    val bedfnm = DARChrPd + s"/$sc.DARChrP.bed"
    writeStrings2File(
      content = gs.map(x => x.mkString("\t")),
      to = bedfnm,
      overwrite = true,
      head = ""
    )
  })

   // 2. perform DAR-ChrP Homer anlaysis
    val sc2homerTask = ptscs.map(sc => {
      val t = new MotifFinderByHomer2(
        inputfnm = DARChrPd + s"/$sc.DARChrP.bed",
        flagfnm = flagd + s"/$sc.DAR.ChrP.homer.done",
        logfnm = logd + s"/$sc.DAR.ChrP.homer.log",
        name = sc + "-DARChrP",
        outd = outHomerd + s"/$sc",
        skip = true,
        check = true,
        homer2 = "/projects/ps-renlab2/szu/softwares/homer/bin",
        seqsize = -1,
        bgfnm = None,
        denovo = "-nomotif",
        mask = "-mask",
        hyperGenomic = "",
        genome = "mm10",
        ncore = 1
      )
      (sc, t)
    }).toMap

    sc2homerTask.par.foreach((sc, t) => {t.run()})
  println("Done.")
}
