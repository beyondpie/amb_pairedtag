import os._
import SZUtils.writeStrings2File
import SZUtils.TaskElement
import SZUtils.SimpleCommandTask
import SZUtils.readTable
import Peak.Peak.readBed
import SZUtils.str2path
import SZUtils.Extensions.*

/*
   Add binarization after SPM quantile 0.25 filtering (2025-01-21)
 */
object ReRunPeakPartitionAfterSPMFiltering {
  // about chromHMM
  val chromHMMJar =
    "/projects/ps-renlab2/szu/softwares/ChromHMM/ChromHMM.jar"

  val chromSizefnm =
    "/projects/ps-renlab2/szu/softwares/ChromHMM/CHROMSIZES/mm10.txt"

  val binSize: Int = 200
  val nstate: Int  = 18
  val ncpu: Int    = 40

  // meta
  val allMods =
    List("ATAC", "H3K27me3", "H3K27ac", "H3K4me1", "H3K9me3")
  val projd     = "/projects/ps-renlab2/szu/projects/amb_pairedtag"
  val workd     = s"${projd}/06.ChromHMM"
  val chromHMMd = s"${projd}/data/chromHMM"
  val peakd     = s"${chromHMMd}/subclass_peak_SPMq25"
  val bBedr     = s"${workd}/out/bBed_bypeak_b${binSize}_SPMq25"
  // os.makeDir(bBedr)
  val cellmarktable = s"${peakd}/cellmarkfiletable.tsv"

  // os.makeDir(outd)
  val logd  = s"${workd}/log"
  val flagd = s"${workd}/flag"

  val binarizeLog  = s"${logd}/binarizeSPMq25.log"
  val binarizeFlag = s"${flagd}/binarizeSPMq25.done"

  def main(args: Array[String]) = {
    val binarizeGenome =
      new SimpleCommandTask(commands = List("java", "-mx500000M",
              "-jar", chromHMMJar, "BinarizeBed", "-b",
              binSize.toString, "-peaks", chromSizefnm, peakd,
              cellmarktable, bBedr), binarizeLog, binarizeFlag,
          skip = true, check = true, verbose = true,
          name = "Binarize signals")
    binarizeGenome.run()

  } // end of main
}   // end of object runChromHMM
