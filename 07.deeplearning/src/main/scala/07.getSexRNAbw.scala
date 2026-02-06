import MetaData.readSubclassMeta
import PairedTagBigWig.BamToBigWig
import SZUtils.path2str
import os.*
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*

@main def GetSexRNABigWig(): Unit = {
  val projd =
    os.Path("/projects/ps-renlab2/szu/projects/amb_pairedtag")
  val workd = projd / "07.deeplearning"
  val logd  = workd / "log"
  val flagd = workd / "flag"
  val RNAbamd   = projd / "data" / "ptRNAbam"
  val ptscMeta  = readSubclassMeta().filter(_.atac)
  val ptscs     = ptscMeta.map(x => x.name.PairedTagName)
  val sexes     = Vector("Male", "Female")
  val outRNAbwd = RNAbamd / "bigwig"

  val scsexbam2bwTask: Vector[(String, String, BamToBigWig)] = {
    ptscs.flatMap(sc => {
      sexes.map(sex => {
        val name           = s"$sc.$sex.RNA.bam2bw"
        val bamfnm: String = RNAbamd / "bam" / s"${sc}-${sex}.bam"
        val bwfnm: String  =
          outRNAbwd / s"${sc}.${sex}.RPKM.s300.bw"
        val logfnm: String  = logd / s"$name.log"
        val flagfnm: String = flagd / s"$name.done"
        val t               = new BamToBigWig(
            inbam = bamfnm,
            bw = bwfnm,
            logfnm = logfnm,
            flagfnm = flagfnm,
            skip = false,
            check = true,
            skipbw = true,
            smoothLen = 300
        )
        (sc, sex, t)
      })
    })}

    scsexbam2bwTask.par.foreach(x => {
      x._3.run()
    })

}
