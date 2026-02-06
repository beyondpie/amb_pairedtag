import os._
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*
import MetaData.TSCCTools.BamCoverage2BigWig
import MetaData.TSCCTools.Bed2Bam
import MetaData.TSCCMeta.projd

object GetATACBamAndBigwigFromBed {
  val bedir = s"${projd}/data/snATAC/subclass_bed_zst"
  val bamd = s"${projd}/data/snATAC/subclass_bam_frombed"
  val bwd = s"${projd}/data/snATAC/subclass_bigwig_bamCoverage"
  val genome = s"${projd}/meta/mm10.chrom.sizes.lite"
  val logd = s"${projd}/data/snATAC/tmp/log"
  val flagd = s"${projd}/data/snATAC/tmp/flag"

  def main(args: Array[String]) = {
    val extSize: Int = 100
    val binSize: Int = 100
    val smoothLength: Int = 300
    val normal = "RPKM"
    val tag = s"ATAC.e${extSize}.bs${binSize}.sm${smoothLength}"
    // set up sc2bedzst map
    val scs: List[String] =
      os.list(os.Path(bedir))
      .map(x => x.baseName.split("\\.").head)
      .toList

    val sc2bedzst: Map[String, String] =
      scs.map(s => (s, s"${bedir}/${s}.bed.zst")).toMap

    val sc2bam: Map[String, String] =
      scs.map(s => (s, s"${bamd}/${s}.srt.bam")).toMap

    val sc2bw: Map[String, String] =
      scs.map(s => (s, s"${bwd}/${s}.${tag}.bw")).toMap

    // set up sc bed2bam task
    // val sc2bamTask = scs.map(s =>
    //   new Bed2Bam(
    //     inbed = sc2bedzst(s),
    //     outbam = sc2bam(s),
    //     genomefnm = genome,
    //     logfnm = s"${logd}/${s}.${tag}.bam.log",
    //     flagfnm = s"${flagd}/${s}.${tag}.bam.done",
    //     skip = true,
    //     check = true
    //   )
    // )

    // sc2bamTask.par.foreach{
    //   t => {
    //     t.run()
    //   }
    // }

    // set up sc bam2bigwig task
    val sc2bwTask = scs.map(s =>
      new BamCoverage2BigWig(
        inbam = sc2bam(s),
        bw = sc2bw(s),
        ext = extSize,
        bs = binSize,
        sm = smoothLength,
        normal = normal,
        ncpu = 1,
        logfnm = s"${logd}/${s}.${tag}.bw.log",
        flagfnm = s"${flagd}/${s}.${tag}.bw.done",
        skip = true,
        check = true,
        skipbw = false
      )
    )
    sc2bwTask.par.foreach{
      t => t.run()
    }

  } // end of main
} // end of object with main
