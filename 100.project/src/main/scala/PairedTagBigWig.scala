package PairedTagBigWig

import Bam.removeDuplicatedPileUp
import GRange.GenomicRange
import MetaData.TSCCMeta.projd
import SZUtils.{TaskElement, writeStrings2File}
import SimpleBigWig.{
  BedGraphElement, loadBigWigOneChr, loadChromosomeData,
  mapBigWigOnRegionOneChr
}

val ptRNAbwd = s"$projd/data/ptRNAbam/bigwig"
val ptDNAbwd = s"$projd/data/ptDNAbam/bigwig"

def getptmodSuffix(mod: String) = mod match {
  case "H3K9me3" => "e100.bs100.sm1000"
  case _         => "e100.bs100.sm300"
}

def loadptRNAbw(sc: String,
  chrs: List[String]): Map[String, List[BedGraphElement]] = {
  loadChromosomeData(s"$ptRNAbwd/${sc}.RPKM.bw", chrs)
}

def loadptRNAbwOneChr(sc: String,
  chr: String): Vector[BedGraphElement] = {
  loadBigWigOneChr(s"$ptRNAbwd/${sc}.RPKM.bw", chr)
}

def loadptDNAbw(sc: String, mod: String,
  chrs: List[String]): Map[String, List[BedGraphElement]] = {
  loadChromosomeData(
      s"$ptDNAbwd/${sc}.${mod}.${getptmodSuffix(mod)}.bw", chrs)
}

def loadptDNAbwOneChr(sc: String, mod: String,
  chr: String): Vector[BedGraphElement] = {
  loadBigWigOneChr(
      s"$ptDNAbwd/${sc}.${mod}.${getptmodSuffix(mod)}.bw", chr)
}

@deprecated("No need now: check latest lightBedGraph Module in bioscala.")
class BigWigScoreTask(val bwf: String,
  val chr2CREs: Map[String, Vector[GenomicRange]],
  val outd: String, val flagd: String, val prefix: String)
    extends TaskElement {
  val skip: Boolean   = true
  val flagfnm         = s"$flagd/${prefix}.done"
  val outf            = s"$outd/${prefix}.csv"
  def runCore(): Unit = {
    val r = 1
      .to(19)
      .map(i => s"chr$i")
      .filter(chr => chr2CREs.contains(chr)).flatMap(chr => {
        val bw = loadBigWigOneChr(bwf, chr)
        if (bw.length < 1) {
          chr2CREs(chr)
            .map(g => new BedGraphElement(g, 0.0))
            .toVector
        } else {
          mapBigWigOnRegionOneChr(bw, chr2CREs(chr),
            statistic = "mean", emptyValue = 0.0,
            keepOrder = true)
        }
      })
      .toList
    writeStrings2File(
        content = r.map(x => x.mkString(",")),
        to = outf,
        overwrite = true,
        head = BedGraphElement.mkHead(",")
    )
  } // end of runcore
}   // end of class BigWigScoreTask

class CleanDuplicatedPileups(
  val inbam: String,
  val outbam: String,
  val logfnm: String,
  val flagfnm: String,
  val cutoff: Int = 10,
  val skip: Boolean = true,
  val check: Boolean = false
) extends TaskElement {
  def runCore(): Unit = {
    println(s"Index ${inbam}.")
    os.proc("samtools", "index", inbam)
      .call(check = check, stdout = os.Path(logfnm),
          stderr = os.Path(logfnm))
    println(s"Remove duplicated pileups for ${inbam}.")
    removeDuplicatedPileUp(inbam, outbam, cutoff = cutoff,
        skip = true)
    println(s"Generated ${outbam}.")
  }
} // end of class CleanDuplicatedPileups


class BamToBigWig(
  val inbam: String,
  val bw: String,
  val logfnm: String,
  val flagfnm: String,
  val skip: Boolean = true,
  val check: Boolean = false,
  val skipbw: Boolean = true,
  val extSize: Int = 100,
  val binSize: Int = 100,
  val samtools: String = "samtools",
  val bamCoverage: String = "bamCoverage",
  val smoothLen: Int = -1
) extends TaskElement {
  def runCore(): Unit = {
    val m       = inbam.split("/").last.split("\\.").head

    // setting smoothLen if not setting
    val sm: Int = if (smoothLen < 0) {
      if (m == "H3K9me3") {
        1000
      } else {
        300
      }
    } else {
      smoothLen
    }

    if (os.exists(os.Path(bw)) && skipbw) {
      println(s"${bw} already exists, and skip it.")
    } else {
      println(s"Index ${inbam}.")
      os.proc(samtools, "index", inbam)
        .call(
            check = false,
            stdout = os.Path(logfnm),
            stderr = os.Path(logfnm)
        )
      println(s"Generate bigwig from ${inbam} with smooth ${sm}.")
      os.proc(
          bamCoverage,
          "-b",
          inbam,
          "-o",
          bw,
          "-p",
          1,
          "-e",
          extSize.toString,
          "-bs",
          binSize.toString,
          "--smoothLength",
          sm.toString,
          "--normalizeUsing",
          "RPKM"
      ).call(
          check = check,
          stdout = os.Path(logfnm),
          stderr = os.Path(logfnm)
      )
      println(s"Generated ${bw}.")
    }
  }
} // end of class BamToBigWig
