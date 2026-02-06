import os._
import MetaData.TSCCTools.MergeBamUsingSamtools
import SZUtils.writeStrings2File
import MetaData.TSCCTools.Bam2Bed
import AllenMetaData.sc2spNN
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*

object getPairedTagSubclassBed4ChromHMM {

  val projd = "/projects/ps-renlab2/szu/projects/amb_pairedtag"
  val bamd = s"${projd}/data/ptDNAbam/bam"
  val chromHMMd = s"${projd}/data/chromHMM"

  val logd = s"${chromHMMd}/tmp/log"
  val flagd = s"${chromHMMd}/tmp/flag"

  val subclasses = os.list(os.Path(bamd)).map(x => x.baseName)

  def setsc2bamlistfnm(hm: String): Map[String, String] = {
    subclasses.map(sc => (sc, s"${logd}/${sc}_${hm}_bams.list")).toMap
  }

  def setsc2bam(hm: String): Map[String, String] = {
    subclasses.map(sc => (sc, s"${bamd}/${sc}/${hm}.srt.bam")).toMap
  }

  def setsc2bedfnm(hm: String): Map[String, String] = {
    subclasses
      .map(sc => (sc, s"${chromHMMd}/pt_subclass/${sc}.${hm}.chromHMM.bed"))
      .toMap
  }

  def writeBamList(sc2bamlist: Map[String, String], hm: String): Unit = {
    subclasses.foreach { sc =>
      {
        val f1 = s"${bamd}/${sc}/${hm}-Male.srt.bam"
        val f2 = s"${bamd}/${sc}/${hm}-Female.srt.bam"
        val fs = List(f1, f2).filter(f => os.exists(os.Path(f)))
        if (fs.length > 0) {
          writeStrings2File(content = fs, to = sc2bamlist(sc))
        }
      }
    }
  } // end of function of writeBamList

  def setMergeBamTask(
    sc2bamlist: Map[String, String],
    sc2bam: Map[String, String],
    hm: String
  ): Map[String, MergeBamUsingSamtools] = {
    subclasses
      .filter(sc => os.exists(os.Path(sc2bamlist(sc))))
      .map { sc =>
        {
          val r = new MergeBamUsingSamtools(
            samtools = "/home/szu/mambaforge/bin/samtools",
            flagfnm = s"${flagd}/${sc}.${hm}.mergesexbam.done",
            logfnm = s"${logd}/${sc}.${hm}.mergesexbam.log",
            fnmofListOfBams = sc2bamlist(sc),
            toBamfile = sc2bam(sc),
            filterMAPQ = false,
            skip = true
          )
          (sc, r)
        }
      }
      .toMap
  }

  def setBam2BedTask(
    sc2bam: Map[String, String],
    sc2bed: Map[String, String],
    hm: String
  ): Map[String, Bam2Bed] = {
    sc2bam
      .map((sc, bam) => {
        val r = new Bam2Bed(
          bedtools = "/home/szu/mambaforge/bin/bedtools",
          flagfnm = s"${flagd}/${sc}.${hm}.bam2bed.done",
          logfnm = s"${logd}/${sc}.${hm}.bam2bed.log",
          bamfnm = sc2bam(sc),
          bedfnm = sc2bed(sc),
          skip = true
        )
        (sc, r)
      })
      .toMap
  }

  def main(args: Array[String]) = {

    // for (hm <- List("H3K27ac", "H3K4me1", "H3K27me3", "H3K9me3")) {
    //   println(hm)
    //   val sc2bam = setsc2bam(hm)
    //   val sc2bed = setsc2bedfnm(hm)

    //   // merge male and female bams
    //   val sc2bamlist = setsc2bamlistfnm(hm)
    //   writeBamList(sc2bamlist, hm)
    //   val sc2mergebam = setMergeBamTask(sc2bamlist, sc2bam, hm)
    //   sc2mergebam.par.foreach((sc, t) => {
    //     println(s"Merge bam for ${sc}.")
    //     t.run()
    //   })
    //   // transform bam to bed
    //   val sc2bamtobed = setBam2BedTask(sc2bam, sc2bed, hm)
    //   sc2bamtobed.par.foreach((sc, t) => {
    //     println(s"Bam2bed for ${sc}.")
    //     t.run()
    //   })
    //   println(s"${hm} is done.")
    // } // end of for loop of 4 histone modifications

    // 2024-10-02
    // Fix bam of "045_OB_STR_CTX_Inh_IMN"
    // especially on H3K27ac
    val sc = "045_OB_STR_CTX_Inh_IMN"
    val hm = "H3K27ac"
    val bamlistfnm = s"${logd}/${sc}_${hm}_bams.list"
    val bamfnm = s"${bamd}/${sc}/${hm}.srt.bam"
    val r = new MergeBamUsingSamtools(
      samtools = "/home/szu/mambaforge/bin/samtools",
      flagfnm = s"${flagd}/${sc}.${hm}.mergesexbam.done",
      logfnm = s"${logd}/${sc}.${hm}.mergesexbam.log",
      fnmofListOfBams = bamlistfnm,
      toBamfile = bamfnm,
      filterMAPQ = false,
      skip = false
    )
    r.run()

    // NN
    val nnsps: List[String] = sc2spNN.values.flatten.toList
    val nnscs: List[String] = sc2spNN.keys.toList
    // merge NN supertype bed to subclass level
    // for (hm <- List("H3K27ac", "H3K4me1", "H3K27me3", "H3K9me3")) {
    //   val sp2bed: Map[String, String] = nnsps.map(
    //     sp => (sp, s"${chromHMMd}/pt_subclass/${sp}.${hm}.chromHMM.bed")
    //   ).toMap
    //   val sc2bed: Map[String, String] = nnscs.map(
    //     sc => (sc, s"${chromHMMd}/pt_subclass/${sc}.${hm}.chromHMM.bed")
    //   ).toMap

    //   nnscs.foreach{
    //     sc => {
    //       println(s"Merge beds for ${sc}_${hm}.")
    //       val spbeds = sc2spNN(sc).map(sp => sp2bed(sp)).filter(x => os.exists(os.Path(x)))
    //       if (spbeds.length < 1) {
    //         println("No sp beds for merging. Skip.")
    //       } else {
    //         val c1 = (List("cat") ::: spbeds).map(x => os.Shellable(List(x)))
    //         os.proc(c1*)
    //           .pipeTo(os.proc("sort", "-k1,1V", "-k2,2n", "-k3,3n"))
    //           .call(check = false, stdout = os.Path(sc2bed(sc)))
    //       } // end of if-else
    //     } // sc foreach function
    //   } // end of foreach
    // } // end of hm list for NN
    for (hm <- List("H3K4me1", "H3K27me3", "H3K9me3")) {
      val sp2bam: Map[String, String] =
        nnsps.map(sp => (sp, s"${bamd}/${sp}/${hm}.srt.bam")).toMap
      val sc2bam: Map[String, String] =
        nnscs.map(sc => (sc, s"${bamd}/${sc}/${hm}.srt.bam")).toMap
      val sc2bamlist: Map[String, String] =
        nnscs.map(sc => (sc, s"${logd}/${sc}_${hm}_bams.list")).toMap
      // update sc2bamlist
      nnscs.foreach { sc =>
        {
          val fs =
            sc2spNN(sc).map(sp => sp2bam(sp)).filter(x => os.exists(os.Path(x)))
          if (fs.length < 1) {
            println(s"${sc}_${hm} has no supertype bams. Skip it")
          } else {
            writeStrings2File(content = fs, to = sc2bamlist(sc))
          }
        }
      }
      // merge bam
      nnscs.par.foreach { sc =>
        {
          if (os.exists(os.Path(sc2bamlist(sc)))) {
            println(s"Merge bams for ${sc} ${hm}.")
            var r = new MergeBamUsingSamtools(
              samtools = "/home/szu/mambaforge/bin/samtools",
              flagfnm = s"${flagd}/${sc}.${hm}.mergespbam.done",
              logfnm = s"${logd}/${sc}.${hm}.mergespbam.log",
              fnmofListOfBams = sc2bamlist(sc),
              toBamfile = sc2bam(sc),
              filterMAPQ = false,
              skip = true
            )
            val d = os.Path(sc2bam(sc)) / ".."
            if (!os.exists(d)) {
              println(s"Create dir: ${d}.")
              os.makeDir(d)
            }
            r.run()
          } // end of if
        } // end of fn inside foreach element.
      } // end of foreach of merging bam
    } // end of hm for nn merging bams

  } // end of Main
}
