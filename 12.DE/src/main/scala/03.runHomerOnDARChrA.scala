// Perform homer based motif enrichment for
// DAR Enhs (ChrA-labeled) CREs per subclass.

// global background

import GRange.{GenomicRange,
  filterByRanges, mouseGenomicRangeOrd}

import SZUtils.writeStrings2File
import Homer2.MotifFinderByHomer2
import MetaData.readSubclassMeta
import Bed.BedElement4
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*

@main def HomerOnDARChrA(): Unit = {
  // 1. create bed files for DAR Enhs per subclass
  val projd = MetaData.TSCCMeta.projd
  val workd = projd + "/12.DE"
  val flagd =  workd + "/flag"
  val logd = workd + "/log"
  val outd = workd + "/out/DARChrAHomer"
  val DARChrAd = workd + "/out/DARChrA"
  val outHomerd = workd + "/out/DARChrAHomer"
  val DARfnm = projd + "/16.celloracle/out/" + "DAR.log2fd0.2.csv"
  val ChrAfnm = projd + "/data/enhancer/ChromHMMChrA.CRE.csv"
  val ptscMeta = readSubclassMeta().filter(x => x.atac)
  val scs = ptscMeta.map(x => x.name.PairedTagName)

  // a. get subclass-level DARChrA
  val allDAR: Vector[BedElement4] =
    os.read.lines
    .stream(os.Path(DARfnm))
    .slice(1, Int.MaxValue)
    .map(x => x.strip().split(","))
    .map(x => new BedElement4(
      GenomicRange.fromUCSCStr(x.head), x.last))
    .toVector

  val allChrA: Vector[BedElement4] = BedElement4.fromTwoCols(
    ChrAfnm, sep = ",", head = false)

  val sc2g: Map[String, Vector[GenomicRange]] = scs.par.map(sc => {
    val DARs = allDAR.filter(_.name == sc).map(_.x).sorted
    val ChrAs = allChrA.filter(_.name == sc).map(_.x).sorted
    val r = filterByRanges(ChrAs, DARs).toVector
    println(s"$sc is done.")
    (sc, r)
  }).toVector.toMap

  // b. save as bed files
  sc2g.foreach((sc, gs) => {
    writeStrings2File(
      content = gs.map(g => g.mkString(sep = "\t")),
      to = DARChrAd + s"/$sc.DARChrA.bed",
      overwrite = true,
      head = ""
    )
  })

  // 2. perform HOMER on these bed files independently
  val sc2homerTask = scs.map(sc => {
    val t = new MotifFinderByHomer2(
      inputfnm = DARChrAd + s"/$sc.DARChrA.bed",
      flagfnm = flagd + s"/$sc.DARChrA.homer.done",
      logfnm = logd + s"/$sc.DARChrA.homer.log",
      name = sc + "-DARChrA",
      outd = outHomerd + s"/$sc",
      skip = false,
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

  sc2homerTask.par.foreach((sc, t) => {
    t.run()
  })
  println("Done.")
}
