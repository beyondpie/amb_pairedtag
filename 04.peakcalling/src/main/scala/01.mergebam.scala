import os._
import scala.collection.parallel.CollectionConverters.*
import MetaData.TSCCMeta.{ptscDNAbamd, ptclusterDNAbamd, ptCellMetaFile}
import MetaData.PairedTagBarcodeMeta
import MetaData.toStringOfAllenIdLabel
import MetaData.TSCCTools.samtools
import MetaData.isNeu
import MetaData.PairedTagBarcodeMeta.readPairedTagCellMetaFile
import MetaData.TSCCMeta.getBarcode2Bam
import scala.collection.immutable.HashMap
import scala.util.Random

/**
 Organize single-cell DNA bam files into different clusters.
 Under the outdir, the results are organized as followed
 - [cell-type] (subclass / supertype)
 - [DNA modalirty] (H3K27ac, H3K27me3, H3K4me1, K3K9me3)
 - [sex] (male / female)
 Under [cell type] / [DNA modality] / [sex] directory,
 we typically have 5 files: all, repA, repB, pseudorepA, pseudorepB
 for pseudorep, we equally partition the dataset.

 All the bam files are sorted based on coordinate.
 Mapping score under 10 will be ignored, except H3K9me3.

 Task before merging bams [3/3]:
 1. check cellmeta is right [Done]
 2. check IMN including cells from non-neuronal cells [Done]
 3. filter low-quality mapping reads[Done]
 */

val cellMeta = readPairedTagCellMetaFile()
val barcode2Bam: Map[String, String] = getBarcode2Bam(cellMeta)

def listBamFilesToFile(listOfBams: List[String], toFile: String): Unit = {
  val out = os.Path(toFile)
  if (os.exists(out)) {
    // println(s"${toFile} exists, and will be removed.")
    os.remove(out)
  }
  os.write(out, listOfBams.mkString("\n"))
}

def mergeBamCore(
  from: String,
  to: String,
  filterMAPQ: Boolean = true,
  MAPQ: Int = 10
): os.CommandResult = {
  if (filterMAPQ) {
    println(s"merge bam with filter MAPQ under ${MAPQ}.")
    os.proc(samtools, "merge", "-", "-b", from)
      .pipeTo(os.proc(samtools, "view", "-b", "-q", MAPQ, "-o", to))
      .call(check = false)
  } else {
    println("merge bam without filterMAPQ.")
    os.proc(samtools, "merge", "-f", "-b", from, "-o", to)
      .call(check = false)
  }
}

// Sorted bams will be sorted bam after merge.
// http://www.htslib.org/doc/samtools-merge.html
def mergeBamFile(
  tagFile: String,
  listOfBamsFile: String,
  toBamFile: String,
  filterMAPQ: Boolean = true,
  MAPQ: Int = 10,
  force: Boolean = false
): Unit = {
  if ((!force) && os.exists(os.Path(tagFile))) {
    println("Detect the tag file: ")
    println(tagFile)
    println("Skip the merging job.")
  } else {
    println("Prepare to generate: ")
    println(toBamFile)
    val res: os.CommandResult =
      mergeBamCore(listOfBamsFile, toBamFile, filterMAPQ, MAPQ)
    if (res.exitCode != 0) {
      println("Error: job failed.")
      println(toBamFile)
    } else {
      println("Merge is done. Now generate tag.")
      os.proc("touch", tagFile).call(check = false)
      println(s"All done for ${toBamFile}.")
    }
  }
}

def key2String(k: List[String]): String = {
  List(toStringOfAllenIdLabel(k(0)), k(1), k(2)).mkString("-")
}
def key2UniPath(k: List[String]): String = {
  List(
    ptclusterDNAbamd,
    "bam",
    toStringOfAllenIdLabel(k(0)),
    s"${k(1)}-${k(2)}.srt.bam"
  ).mkString("/")
}
def key2ListOfBamPath(k: List[String]): String = {
  List(
    ptclusterDNAbamd,
    "meta",
    s"${key2String(k)}.barcode.bams.txt"
  ).mkString("/")
}

def key2TagPath(k: List[String]): String = {
  List(
    ptclusterDNAbamd,
    "tag",
    s"${key2String(k)}.merge.done"
  ).mkString("/")
}

def getGroup(
  x: PairedTagBarcodeMeta,
  cellClass: String,
  groupBy: String
): List[String] = {
  List(
    cellClass match {
      case "neu" => x.annot.sc.get
      case _     => x.annot.sp.get
    },
    x.exp.modularity,
    groupBy match {
      case "sex"       => x.sample.sex
      case "rep"       => x.sample.rep
      case "pseusorep" => x.sample.rep
    }
  )
}

abstract class GroupMergeBam {
  def group2Barcodes: Map[List[String], List[String]]
  // barcodes hare are single-cell bam files
  val g2b = group2Barcodes

  def group2UniBam: Map[List[String], String] = {
    g2b.map((k, v) => (k, key2UniPath(k)))
  }
  val g2u = group2UniBam

  def group2ListOfBams: Map[List[String], String] = {
    g2b.map((k, v) => (k, key2ListOfBamPath(k)))
  }
  val g2l = group2ListOfBams

  def group2Tag: Map[List[String], String] = {
    g2b.map((k, v) => (k, key2TagPath(k)))
  }
  val g2t = group2Tag

  def apply(n: Int = Int.MaxValue, MAPQ: Int = 10, force: Boolean = false): Unit = {
    println("prepare directories.")
    g2l.values
      .map(v => os.Path(v) / os.up)
      .toList
      .distinct
      .foreach(dir => {
        if (!os.exists(dir)) {
          os.makeDir(dir)
        }
      })
    println("start generating merged bam files.")
    g2l.take(n).par.foreach { (k, v) =>
      {
        if (k(1) != "H3K9me3") {
          mergeBamFile(g2t(k), v, g2u(k), true, MAPQ, force)
        } else {
          println("H3K9me3 do not need to filter.")
          mergeBamFile(g2t(k), v, g2u(k), false, MAPQ, force)
        }
      }
    }
  }
  def updateModality(
    n: Int = Int.MaxValue,
    MAPQ: Int = 10,
    modality: String = "H3K9me3",
    force: Boolean = true
  ): Unit = {
    println(s"Update Modality on: ${modality}.")
    g2l.filter((k, v) => k(1) == modality).take(n).par.foreach { (k, v) =>
      {
        if (k(1) != "H3K9me3") {
          mergeBamFile(g2t(k), v, g2u(k), true, MAPQ, force)
        } else {
          println("H3K9me3 do not need to filter.")
          mergeBamFile(g2t(k), v, g2u(k), false, MAPQ, force)
        }
      }
    }
  }
}

def mapBarcodeToBam(
  x: PairedTagBarcodeMeta,
  cellClass: String,
  groupBy: String
): String = {
  barcode2Bam(x.barcode)
}

def mapSexRepToBam(
  x: PairedTagBarcodeMeta,
  cellClass: String,
  groupBy: String
): String = {
  val sex = getGroup(x, cellClass, "sex")
  val rep = List(sex(0), sex(1), x.sample.rep)
  key2UniPath(rep)
}

class Group4MergeBam(
  cellClass: String,
  updateListOfBam: Boolean = true,
  groupBy: String = "rep",
  pairedTag2BamMap: (PairedTagBarcodeMeta, String, String) => String
) extends GroupMergeBam {
  def group2Barcodes = {
    cellMeta
      .filter(x => x.annotQuality == "Good")
      .filter(x => {
        cellClass match {
          case "neu" => isNeu(x)
          case _     => !isNeu(x)
        }
      })
      .groupMap(x => getGroup(x, cellClass, groupBy))(x =>
        pairedTag2BamMap(x, cellClass, groupBy)
      )
      .map((k, v) => (k, v.distinct))
  }

  if (updateListOfBam) {
    println("Update List of Bams for each group.")
    g2l.keys.par.foreach(k => listBamFilesToFile(g2b(k), g2l(k)))
  } else {
    println("NOT update List of Bams for each group.")
  }
}

object mergeBamClusterLevel {
  def main(args: Array[String]) = {
    // val neuRepGroup = new Group4MergeBam("neu", false, "rep", mapBarcodeToBam)
    // neuRepGroup()
    // neuRepGroup.updateModality(modality = "H3K9me3", force = true)
    // val nnRepGroup = new GroupSexRep4MergeBam("nn", false, "rep", mapBarcodeToBam)
    // nnRepGroup(n = Int.MaxValue, MAPQ = 10)
    // nnRepGroup.updateModality(modality = "H3K9me3", force = true)
    // val neuSexGroup =
    //  new Group4MergeBam("neu", updateListOfBam = true, groupBy = "sex", mapSexRepToBam)
    // neuSexGroup.apply(n = Int.MaxValue, MAPQ = 10, force = false)
    val nnSexGroup =
      new Group4MergeBam("nn", updateListOfBam = true, groupBy = "sex", mapSexRepToBam)
    nnSexGroup.apply(n = Int.MaxValue, MAPQ = 10, force = false)
  }
}
