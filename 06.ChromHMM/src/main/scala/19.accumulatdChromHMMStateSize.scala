import os._
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters._
import MetaData.TSCCMeta
import MetaData.readSubclassMeta
import GRange.{BedToolGenomicRangeMerge, GenomicRange}
import SZUtils.writeStrings2File
import ChromHMM.loadDenseBed

@main def countAccumulateStateSize(): Unit = {
  val projd      = "/projects/ps-renlab2/szu/projects/amb_pairedtag"
  val workd      = projd + "/06.ChromHMM"
  val denseBedir = workd + "/out/updateDenseBed"
  val outfnm     = workd + "/out/stat.ChromHMMState2size.csv"
  val ptscMetafnm =
    projd + "/meta/PairedTagSubclassMetaFromCellMetaChromHMM.csv"
  val ptscMeta = readSubclassMeta(ptscMetafnm).filter(_.atac)
  val ptscs    = ptscMeta.map(_.name.PairedTagName)
  val states =
    Vector("Chr-A", "Chr-O", "Chr-B", "Chr-P", "Chr-R", "Hc-H", "Hc-P")
  val mouseGenomefnm    = projd + "/meta/mm10.autosome.chromsize.tsv"
  val ordMouseAutosomes = 1.to(19).map(i => s"chr$i").toList

  // 1. load ChromHMM annotations on subclass level
  def loadsc2state(state: String): Map[String, List[GenomicRange]] = {
    ptscs
      .par
      .map(sc => {
        val densebedfnm = denseBedir + s"/${sc}_18_dense.bed"
        val r = loadDenseBed(densebedfnm)._2
          .filter(x => x.name == state)
          .map(_.g)
          .toList
        (sc, r)
      })
      .toList
      .toMap
  }

  val state2scs = states.map(s => (s, loadsc2state(s))).toMap
  println("Finish loading the genomic ranges for all the states.")

  // 2. at lease one logic to summarize the
  // state coverage
  val state2ranges: Map[String, List[GenomicRange]] = state2scs
    .map((s, sc2g) => {
      println(s"merge genomic ranges for state: $s.")
      val r: List[GenomicRange] = sc2g
        .values
        .toVector
        .par
        .reduce((x1, x2) => {
          val x = (x1 ::: x2)
          BedToolGenomicRangeMerge.merge(x, ordMouseAutosomes)
        })
      (s, r)
    })
    .toMap

  val state2size: Map[String, Int] = state2ranges.map((s, gs) => {
    (s, gs.map(x => x.endTo - x.startFrom).sum)
  })

  // 3. get mouse genomic coverage
  val mouseAutoChromsSize =
    os.read.lines
    .stream(os.Path(mouseGenomefnm))
    .map(x => x.strip().split("\t"))
    .map(x => x(1).toDouble)
    .toVector
    .sum

  val state2r: Vector[List[String]] = state2size
    .map((s, size) =>
      List(s, size.toString,
          (size.toDouble * 100.0 / mouseAutoChromsSize.toDouble).toString))
    .toVector

  // val state2r =
  //   os.read.lines.stream(os.Path(outfnm))
  //     .drop(1)
  //     .map(x => x.strip().split(","))
  //     .map(x => (x(0), x(1).toDouble))
  //     .map(x => List(x._1, x._2.toString,
  //       (x._2 * 100.0 / mouseAutoChromsSize).toString ))
  //     .toList

  // 4. output results
  writeStrings2File(
      content = state2r.map(x => x.mkString(",")),
      to = outfnm,
      overwrite = true,
      head = List("ChromHMMState", "size", "genomeCoverage").mkString(
          ",")
  )
  // takes about 1.5 hours
  println("Done")
}
