// Perform homer based motif enrichment for
// DAR Chr-O CREs per subclasses

import GRange.{GenomicRange,
  filterByRanges, mouseGenomicRangeOrd}

import SZUtils.writeStrings2File
import Homer2.MotifFinderByHomer2
import MetaData.readSubclassMeta
import Bed.BedElement4
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*

@main def HomerOnDARChrO(): Unit = {
  val projd = MetaData.TSCCMeta.projd
  val workd = projd + "/12.DE"
  val flagd = workd + "/flag"
  val logd = workd + "/log"
  //val outd = workd + "/out/DARChrOHomer"
  val DARChrOd = workd + "/out/DARChrO"
  val outHomerd = workd + "/out/DARChrOHomer"
  val DARfnm = projd + "/16.celloracle/out/" + "DAR.log2fd0.2.csv"
  val ChrOfnm = projd + "/data/enhancer/ChromHMMChrO.CRE.csv"

  // val scs = MetaData.representPaiedTagSubclass
  val ptscMeta = readSubclassMeta().filter(_.atac)
  val scs = ptscMeta.map(_.name.PairedTagName)

  // a. get subclass-level DARChrO
  val allDAR: Vector[BedElement4] =
    os.read.lines
      .stream(os.Path(DARfnm))
      .slice(1, Int.MaxValue)
      .map(x => x.strip().split(","))
      .map(x => new BedElement4(
        GenomicRange.fromUCSCStr(x.head), x.last))
      .toVector

  val allChrO: Vector[BedElement4] = BedElement4.fromTwoCols(
    ChrOfnm, sep = ",", head = false)

  val sc2g: Map[String, Vector[GenomicRange]] = scs.par.map(sc => {
    val DARs = allDAR.filter(_.name == sc).map(_.x).sorted
    val ChrOs = allChrO.filter(_.name == sc).map(_.x).sorted
    val r = filterByRanges(ChrOs, DARs).toVector
    println(s"$sc DARChrO collection is done.")
    (sc, r)
  }).toVector.toMap

  // b. save as bed files
  sc2g.foreach((sc, gs) => {
    writeStrings2File(
      content = gs.map(g => g.mkString(sep = "\t")),
      to = DARChrOd + s"/$sc.DARChrO.bed",
      overwrite = true,
      head = ""
    )
  })

  // 2. perform HOMER on these bed files independently
  val sc2homerTask = scs.map(sc => {
    val t = new MotifFinderByHomer2(
      inputfnm = DARChrOd + s"/$sc.DARChrO.bed",
      flagfnm = flagd + s"/$sc.DARChrO.homer.done",
      logfnm = logd + s"/$sc.DARChrO.homer.log",
      name = sc + "-DARChrO",
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
