import os._
import MetaData.TSCCMeta
import GRange.{GenomicRange, mouseGenomicRangeOrd}
import SZUtils.writeStrings2File
import MetaData.readSubclassMeta
import ChromHMM.loadDenseBed

@main def intersectNotDeterminedRegion: Unit = {
  val projd      = "/projects/ps-renlab2/szu/projects/amb_pairedtag"
  val workd      = projd + "/06.ChromHMM"
  val denseBedir = workd + "/out/updateDenseBed"
  val outNDBedir = workd + "/out/NDBed"
  val outd       = workd + "/out"
  val tmpbedfnm  = s"${outd}/tmp0.bed"
  val intersectBedBin =
    "/home/szu/mambaforge/envs/sa2stable/bin/intersectBed"
  if (!os.exists(os.Path(outNDBedir))) {
    os.makeDir(os.Path(outNDBedir))
  }
  val ptscMeta = readSubclassMeta().filter(_.atac)
  val ptscs    = ptscMeta.map(_.name.PairedTagName)

  // 1. load subclass-level ND regions
  // val sc2ND: Map[String, Vector[GenomicRange]] = ptscs.map(sc => {
  //   val densebedfnm = denseBedir + s"/${sc}_18_dense.bed"
  //   val r = loadDenseBed(densebedfnm)
  //     ._2
  //     .filter(x => x.name == "ND")
  //     .map(x => x.g)
  //     .toVector
  //   (sc, r)
  // }).toMap

  // 2. save subclass-level ND into seperate bed files.
  // sc2ND.foreach((sc, gs) => {
  //   val bedfnm = outNDBedir + s"/$sc.ND.bed"
  //   val r = gs.sorted
  //   writeListOfString2File(
  //     content = r.map(x => x.mkString("\t")),
  //     to = bedfnm,
  //     overwrite = true,
  //     head = ""
  //   )
  // })

  // 3. use bedtools intersect to sequentially
  //    get the intersections
  def useBedTools4Intersect(x: String, y: String, out: String,
    bin: String): Unit = {
    os.proc(
        bin,
        "-a",
        x,
        "-b",
        y
    ).call(
        check = true,
        stdout = os.Path(out)
    )
    println(s"$out is generated")
  }

  // copy first ptscs ND file as initial tmpfnm
  // val sc0 = ptscs(0)
  // val fnm0 = outd + "/tmp0.bed"
  // if (!os.exists(os.Path(fnm0))) {
  //   os.proc("cp", outNDBedir + s"/${sc0}.ND.bed", fnm0).call(check = true)
  // }

  // ptscs.foldLeft(0)((curId, sc) => {
  //   val tmpbedfnm = s"$outd/tmp$curId.bed"
  //   val nextId = curId + 1
  //   useBedTools4Intersect(
  //     x = tmpbedfnm,
  //     y = outNDBedir + s"/$sc.ND.bed",
  //     out = outd + s"/tmp$nextId.bed",
  //     bin = intersectBedBin
  //   )
  //   // remove previous file.
  //   os.remove(os.Path(tmpbedfnm))
  //   nextId
  // })

  val fnm = outd + s"/tmp${ptscs.length}.bed"
  val r = os.read.lines
    .stream(os.Path(fnm))
    .map(x => x.strip().split("\t"))
    .map(x => new GenomicRange(x(0), x(1).toInt, x(2).toInt))
    .toVector
    .sorted(using mouseGenomicRangeOrd)
  val s: Double = r.groupMapReduce(x => x.chrom)(x => (x.endTo - x.startFrom))(
      (x, y) => { x + y }).foldLeft(0)((s, x) => {
        s + x._2
      }).toDouble
  println(s"Total genomic size: $s.")

  val gs = os.read.lines.stream(os.Path(projd + "/meta/mm10.chrom.sizes.lite"))
    .map(x => x.strip().split("\t"))
    .filter(x => x(0) != "ChrX" || (x(0) != "ChrY"))
    .map(x => x(1).toDouble)
    .toVector.sum

  val ratio: Double = s * 100.0 / gs 
  println(s"Capture $ratio percnt of mouse genome.")
}
