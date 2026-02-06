// Perform homer based motif enrichment for
// All Chr-O CREs per subclass

import GRange.{GenomicRange,
  filterByRanges, mouseGenomicRangeOrd}

import SZUtils.writeStrings2File
import Homer2.MotifFinderByHomer2
import MetaData.readSubclassMeta
import Bed.BedElement4
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*

@main def HomerOnChrO(): Unit = {
  val projd = MetaData.TSCCMeta.projd
  val workd = projd + "/12.DE"
  val flagd = workd + "/flag"
  val logd = workd + "/log"
  val ChrOd = workd + "/out/ChrO"
  val outHomerd = workd + "/out/ChrOHomer"
  val ChrOfnm = projd + "/data/enhancer/ChromHMMChrO.CRE.csv"

  val ptscMeta = readSubclassMeta().filter(_.atac)
  val scs = ptscMeta.map(_.name.PairedTagName)

  val allChrO: Vector[BedElement4] = BedElement4.fromTwoCols(
    ChrOfnm, sep = ",", head = false)

  val sc2g: Map[String, Vector[GenomicRange]] = scs.par.map(sc => {
    val ChrOs = allChrO.filter(_.name == sc).map(_.x).sorted
    println(s"$sc ChrO collection is done.")
    (sc, ChrOs)
  }).toVector.toMap

  // b. save as bed files
  sc2g.foreach((sc, gs) => {
    writeStrings2File(
      content = gs.map(g => g.mkString(sep = "\t")),
      to = ChrOd + s"/$sc.ChrO.bed",
      overwrite = true,
      head = ""
    )
  })

  // 2. perform HOMER on these bed files independently
  val sc2homerTask = scs.map(sc => {
    val t = new MotifFinderByHomer2(
      inputfnm = ChrOd + s"/$sc.ChrO.bed",
      flagfnm = flagd + s"/$sc.ChrO.homer.done",
      logfnm = logd + s"/$sc.ChrO.homer.log",
      name = sc + "-ChrO",
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
