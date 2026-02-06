import os._
import SZUtils.writeStrings2File
import SZUtils.TaskElement
import SZUtils.SimpleCommandTask
import SZUtils.readTable
import Peak.Peak.readBed
import SZUtils.str2path
import SZUtils.Extensions.*
import AllenMetaData.Allen.AllenClMeta
import AllenMetaData.Allen.subclassNameCEMBATAC
import AllenMetaData.Allen.subclassNamePairedTag
import AllenMetaData.Allen.subclassNameDNAMeth

/* Run ChromHMM using the partition defined by each peaks.
   - Only autosome included
   - Only subclasses with all the modalities are used
     - 151 subclasses included finally (14 removed)
 */
object RunChromHMM_PeakPartition {
  // about chromHMM
  val chromHMMJar =
    "/projects/ps-renlab2/szu/softwares/ChromHMM/ChromHMM.jar"

  val chromSizefnm =
    "/projects/ps-renlab2/szu/softwares/ChromHMM/CHROMSIZES/mm10.txt"

  val binSize: Int = 200
  val nstate: Int  = 18
  val ncpu: Int    = 40

  // meta
  val allMods = List("ATAC", "H3K27me3", "H3K27ac", "H3K4me1", "H3K9me3")
  val projd   = "/projects/ps-renlab2/szu/projects/amb_pairedtag"
  val workd   = s"${projd}/06.ChromHMM"
  val chromHMMd = s"${projd}/data/chromHMM"
  val peakd     = s"${chromHMMd}/subclass_peak"
  val hstpeakd  = s"${projd}/data/pairedtag_peak/subclass_peak"
  val atacpeakd = s"${projd}/data/snATAC/sa2.subclassv3.final.peak.srt"
  val bBedr     = s"${workd}/out/bBed_bypeak_b${binSize}"
  // os.makeDir(bBedr)
  val oldTableFile  = s"${chromHMMd}/all_bed/cellmarkfiletable.txt"
  val cellmarktable = s"${peakd}/cellmarkfiletable.tsv"

  val outd  = s"${workd}/out/model_bypeak_b${binSize}_s${nstate}"
  // os.makeDir(outd)
  val logd  = s"${workd}/log"
  val flagd = s"${workd}/flag"

  val binarizeLog    = s"${logd}/binarize.log"
  val binarizeFlag   = s"${flagd}/binarize.done"
  val learnModelLog  = s"${logd}/learnModelLog.log"
  val learnModelFlag = s"${flagd}/learnModelLog.done"

  // Allen cluster meta
  val scPairedTag2ATAC =
    AllenClMeta
      .map(x => (x.subclass.id_label, x.subclass.label))
      .map(x =>
        (subclassNamePairedTag(x._1), subclassNameCEMBATAC(x._2)))
      .toMap

  val scPairedTag2DNAMethy =
    AllenClMeta
      .map(x => (x.subclass.id_label, x.subclass.label))
      .map(x => (subclassNamePairedTag(x._1), subclassNameDNAMeth(x._2)))
      .toMap

  // functions
  def getSubclassesWithoutFullModalities: List[String] = {
    val scs = readTable(oldTableFile, sep = "\t", head = false)
      .map(x => x(0))
      .distinct
    scs.filter(sc => {
      var state = false
      for (i <- allMods) {
        i match {
          case "ATAC" => {
            val tsc = scPairedTag2ATAC(sc)
            val fnm = s"${atacpeakd}/${tsc}.bed"
            if (!os.exists(fnm)) {
              println(s"no ATAC peaks in ${sc}.")
              state = true
            }
          }
          case "H3K27ac" => {
            val fnm = s"${hstpeakd}/${sc}-${i}.BestSPM.peak"
            if (!os.exists(fnm)) {
              println(s"no ${i} peaks in ${sc}.")
              state = true
            }
          }
          case _ => {
            val fnm = s"${hstpeakd}/${sc}-${i}.bedtoolmerge.peak"
            if (!os.exists(fnm)) {
              println(s"no ${i} peaks in ${sc}.")
              state = true
            }
          }
        }
      }
      state
    })
  }

  def createCellMarkTableFile(): Unit = {
    val oldContents =
      readTable(oldTableFile, sep = "\t", head = false)
    val rmscs = getSubclassesWithoutFullModalities
    val contents = oldContents
      .map(x => List(x(0), x(1), x(2).replace("chromHMM", "peak")))
      .filter(x => !rmscs.contains(x(0)))
      .map(x => x.mkString(sep = "\t"))
    writeStrings2File(contents, to = cellmarktable,
        overwrite = true, head = "")
  }

  def isSexChromosome(chr: String): Boolean = {
    chr.contains("chrX") || chr.contains("chrY")
  }

  def createHistonePeakBed(sc: String, h: String): Unit = {
    val fromfnm = h match {
      case "H3K27ac" => s"${sc}-${h}.BestSPM.peak"
      case _         => s"${sc}-${h}.bedtoolmerge.peak"
    }
    val peakInBed =
      readTable(s"${hstpeakd}/${fromfnm}", "\t", head = true)
        .map(x => List(x(0), x(1), x(2), x(3)))
        .filter(x => !isSexChromosome(x(0)))
        .map(x => x.mkString("\t"))

    val tofnm = s"${peakd}/${sc}.${h}.peak.bed"
    writeStrings2File(content = peakInBed, to = tofnm,
        overwrite = true, head = "")
  }

  def createATACPeakBed(sc: String): Unit = {
    val tsc     = scPairedTag2ATAC(sc)
    val fromfnm = s"${atacpeakd}/${tsc}.bed"
    val tofnm   = s"${peakd}/${sc}.ATAC.peak.bed"
    val rawPeaks = readTable(fromfnm, "\t", true, head = false)
      .filter(x => !isSexChromosome(x(0)))
      .map(x => x.mkString("\t"))
    writeStrings2File(rawPeaks, to = tofnm, overwrite = true,
        head = "")
  }

  def runCreatePeakBed(scs: List[String]): Unit = {
    for (
      sc <- scs;
      h  <- allMods
    ) {
      if (h == "ATAC") {
        createATACPeakBed(sc)
      } else {
        createHistonePeakBed(sc, h)
      }
    }
  } // end of runCreatePeakBed

  def main(args: Array[String]) = {
    // 1. prepare binarzie by peak
    // createCellMarkTableFile()
    // // prepare peaks for genome partition
    // val allCombs = readTable(cellmarktable, sep = "\t", head = false)
    // allCombs.foreach { x =>
    //   {
    //     println(s"${x(0)} ${x(1)}")
    //     x(1) match {
    //       case "ATAC" => createATACPeakBed(x(0))
    //       case _      => createHistonePeakBed(x(0), x(1))
    //     }
    //   }
    // }

    // 2. start to binarize by peak
    // val binarizeGenome =
    //   new SimpleCommandTask(commands = List("java", "-mx500000M", "-jar",
    //           chromHMMJar, "BinarizeBed", "-b", binSize.toString, "-peaks",
    //           chromSizefnm, peakd, cellmarktable, bBedr), binarizeLog,
    //       binarizeFlag, skip = true, check = true, verbose = true,
    //       name = "Binarize signals")
    // binarizeGenome.run()

    // 3. run chromHMM
    // this takes time and let's run this step in screen
    val learnModel =
      new SimpleCommandTask(commands = List("java", "-mx800000M", "-jar",
              chromHMMJar, "LearnModel", "-p", ncpu.toString, bBedr, outd,
              nstate.toString, "mm10"), learnModelLog, learnModelFlag,
          skip = true, check = true, verbose = true,
          name = "Learn ChromHMM Model")
    learnModel.run()
  } // end of main
}   // end of object runChromHMM
