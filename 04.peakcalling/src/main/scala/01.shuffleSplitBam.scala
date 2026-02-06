import os._
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*
import MetaData.{PairedTagBarcodeMeta, toStringOfAllenIdLabel}
import MetaData.{TSCCMeta, TSCCTools}
import SZUtils.writeStrings2File
import MetaData.TSCCMeta.{ptclusterDNAbamd, ptscDNAbamd}
import scala.util.Random
import MetaData.isNeu

def readPairedTagCellMetaFile(path:String = TSCCMeta.ptCellMetaFile): List[PairedTagBarcodeMeta] = {
  os.read.lines
    .stream(os.Path(path))
    .slice(1, Int.MaxValue)
    .toList
    .map(x => PairedTagBarcodeMeta(x))
}

object ShuffleSplitBam{
  // paths
  val prefix = "/tscc/projects/ps-renlab2/szu"
  val projd = s"${prefix}/projects/amb_pairedtag"
  val cellMetafnm = s"${projd}/meta/pairedtag.cell.meta.all.240626.csv"
  val workd = s"${projd}/04.peakcalling"
  val flagd = s"${workd}/flag/shuffleSplitBam"
  val splitBarcoded = s"${workd}/out/shuffleSplitBam/splitBarcoded"
  val rawBamd = s"${projd}/data/ptDNAbam/bam"
  val outBamd = s"${projd}/data/ptDNAbam/bam"
  val subsetBamScript = s"${projd}/extratools/singlecell/subset_bam.sh"
  val tmpdir = s"${prefix}/tmpdir"
  val rscd = s"${workd}/src/main/resource"

  List(flagd, splitBarcoded, tmpdir).foreach(
    x => {
      if(!os.exists(os.Path(x))) {
        os.makeDir.all(os.Path(x))
    }
   })

  // meta
  lazy val cellMeta = readPairedTagCellMetaFile(cellMetafnm).
    filter(x => x.annotQuality == "Good")

  // 1. get group to shuffled barcodes
  lazy val neuG2Shfb =
    cellMeta.filter(isNeu)
    .groupMap(
      x => (
        toStringOfAllenIdLabel(x.annot.sc.get),
        List(x.exp.modularity, x.sample.sex).mkString("-")
      ))(x => x.barcode)
    .map((k, v) => {
      val b = Random.shuffle(v).splitAt(1.max(v.length / 2))
      List(
        ((k._1, k._2 + "-shufA") -> b._1),
        ((k._1, k._2 + "-shufB") -> b._2)
      )
    })
  .flatten.toMap

  // 2. write the shuffled barcodes to splitBarcded
  lazy val neuG2Shfnm = neuG2Shfb.map( (k, v) =>
    k -> s"${splitBarcoded}/${k._1}-${k._2}.barcode.txt"
  )
  neuG2Shfb.foreach( (k, v) => {
    writeStrings2File(v, neuG2Shfnm(k), overwrite = true)
  })

  // 3. get group to flag and target bam file
  lazy val neuG2Shflag = neuG2Shfb.map( (k, v) =>
    k -> s"${flagd}/${k._1}-${k._2}.done"
  )
  lazy val neuG2Shfbam = neuG2Shfb.map((k,v) =>
    k -> s"${outBamd}/${k._1}/${k._2}.shuffle.srt.bam"
  )

  // 4. get group to source bam file
  lazy val neuG2Rawbam = neuG2Shfb.map( (k, v) => {
    val prefix = k._2.split("-").take(2).mkString("-")
    k -> s"${rawBamd}/${k._1}/${prefix}.srt.bam"
  })

  // 5. run split
  neuG2Shfnm.par.foreach { (k, v) => {
    os.proc(
      "bash", subsetBamScript,
      neuG2Rawbam(k),
      neuG2Shfnm(k),
      neuG2Shfbam(k),
      tmpdir
    ).call(check = true)
    os.proc("touch", neuG2Shflag(k)).call(check = false)
    println(s"shuffle done: ${k._1}-${k._2}")
  }}

  // 6. check result
  lazy val neuJobStatus = neuG2Shflag.map(
    (k, v) => (k, os.exists(os.Path(v)))
  )

  // 7. save group4shuffle to resource for peak calling
  lazy val neuGStrs = neuG2Shfnm.map((k, v) => {
    val k2 = k._2.split("-")
    List(k._1, k2.head, k2.tail.mkString("-")).mkString(",")
  }).toList
  writeStrings2File(neuGStrs, s"${rscd}/neu.subclass.SexShufRep.record.csv")


  // 1. get group to shuffled barcodes
  lazy val NNG2Shfb =
    cellMeta.filter(x => !isNeu(x))
    .groupMap(
      x => (
        toStringOfAllenIdLabel(x.annot.sp.get),
        List(x.exp.modularity, x.sample.sex).mkString("-")
      ))(x => x.barcode)
    .map((k, v) => {
      val b = Random.shuffle(v).splitAt(1.max(v.length / 2))
      List(
        ((k._1, k._2 + "-shufA") -> b._1),
        ((k._1, k._2 + "-shufB") -> b._2)
      )
    })
  .flatten.toMap

  // 2. write the shuffled barcodes to splitBarcded
  lazy val NNG2Shfnm = NNG2Shfb.map( (k, v) =>
    k -> s"${splitBarcoded}/${k._1}-${k._2}.barcode.txt"
  )
  NNG2Shfb.foreach( (k, v) => {
    writeStrings2File(v, NNG2Shfnm(k), overwrite = true)
  })

  // 3. get group to flag and target bam file
  lazy val NNG2Shflag = NNG2Shfb.map( (k, v) =>
    k -> s"${flagd}/${k._1}-${k._2}.done"
  )
  lazy val NNG2Shfbam = NNG2Shfb.map((k,v) =>
    k -> s"${outBamd}/${k._1}/${k._2}.shuffle.srt.bam"
  )

  // 4. get group to source bam file
  lazy val NNG2Rawbam = NNG2Shfb.map( (k, v) => {
    val prefix = k._2.split("-").take(2).mkString("-")
    k -> s"${rawBamd}/${k._1}/${prefix}.srt.bam"
  })

  // 5. run split
  NNG2Shfnm.par.foreach { (k, v) => {
    os.proc(
      "bash", subsetBamScript,
      NNG2Rawbam(k),
      NNG2Shfnm(k),
      NNG2Shfbam(k),
      tmpdir
    ).call(check = true)
    os.proc("touch", NNG2Shflag(k)).call(check = false)
    println(s"shuffle done: ${k._1}-${k._2}")
  }}

  // 6. check result
  lazy val NNJobStatus = NNG2Shflag.map(
    (k, v) => (k, os.exists(os.Path(v)))
  )

  // 7. save group4shuffle to resource for peak calling
  lazy val NNGStrs = NNG2Shfnm.map((k, v) => {
    val k2 = k._2.split("-")
    List(k._1, k2.head, k2.tail.mkString("-")).mkString(",")
  }).toList
  writeStrings2File(NNGStrs, s"${rscd}/NN.supertype.SexShufRep.record.csv")

  def main(args: Array[String]) = {
    println("Hello")
  }
}
