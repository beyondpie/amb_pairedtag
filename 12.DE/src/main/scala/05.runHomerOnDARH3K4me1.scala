// Perform Homer based motif enrichment for
// DAR Chr-P  H3K4me1-based peaks per subclass.

import GRange.GenomicRange
import os._
import SZUtils.writeStrings2File
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*
import Homer2.MotifFinderByHomer2
import MetaData.readSubclassMeta
import GRange.mouseGenomicRangeOrd

@main def HomerOnDARH3K4me1Peaks(): Unit = {
  val projd = MetaData.TSCCMeta.projd
  val workd = projd + "/12.DE"
  val flagd = workd + "/flag"
  val logd = workd + "/log"

  val DARH3K4me1d = workd + "/out/DARH3K4me1"
  if(!os.exists(os.Path(DARH3K4me1d))) {
    os.makeDir(os.Path(DARH3K4me1d))
  }
  val outHomerd = workd + "/out/DARH3K4me1Homer"
  if(!os.exists(os.Path(outHomerd))){
    os.makeDir(os.Path(outHomerd))
  }

  val H3K4me1DEd = workd + "/out/H3K4me1_sa2LRT_DE"

  //val scs = MetaData.representPaiedTagSubclass
  val ptscMeta = readSubclassMeta().filter(_.atac)
  val scs = ptscMeta.map(_.name.PairedTagName)


  // 1. get bed file for each subclass
  val sc2DAR: Map[String, Vector[GenomicRange]] = scs.map(sc => {
    val fnm = H3K4me1DEd + s"/${sc}_sc_diffpeaks_H3K4me1.txt"
    val r = os.read.lines
      .stream(os.Path(fnm))
      .drop(1)
      .map(x => x.strip().split(","))
      .filter(x => x(1).toDouble >= 0.2)
      .filter(x => x.last.toDouble <= 0.05)
      .map(x => GenomicRange.fromUCSCStr(x(0)))
      .filter(x => x.chrom != "chrX" && x.chrom != "chrY")
      .toVector
      .sorted(using mouseGenomicRangeOrd)
    (sc, r)
  }).toMap

  sc2DAR.foreach((sc, gs) => {
    println(s"Save $sc H3K4me1 DAR results as bed file.")
    writeStrings2File(
      content = gs.map(x => x.mkString("\t")),
      to = DARH3K4me1d + s"/$sc.DAR.H3K4me1Peak.bed",
      overwrite = true,
      head = ""
    )
  })

  // 2. perform HOMER on these bed files
//  val sc2homerTask = scs.map(sc => {
//    val t = new MotifFinderByHomer2(
//      inputfnm = DARH3K4me1d + s"/$sc.DAR.H3K4me1Peak.bed",
//      flagfnm = flagd + s"/$sc.DAR.H3K4me1Peak.homer.done",
//      logfnm = logd + s"/$sc.DAR.H3K4me1Peak.homer.log",
//      name = sc + "-DARH3K4me1Peak",
//      outd = outHomerd + s"/$sc",
//      skip = false,
//      check = true,
//      homer2 = "/projects/ps-renlab2/szu/softwares/homer/bin",
//      seqsize = -1,
//      bgfnm = None,
//      denovo = "-nomotif",
//      mask = "-mask",
//      hyperGenomic = "",
//      genome = "mm10",
//      ncore = 1
//    )
//    (sc, t)
//  }).toMap

//  sc2homerTask.par.foreach((sc, t) => {t.run()})
  println("Done.")
}
