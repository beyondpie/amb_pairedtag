import os._
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*
import SZUtils.TaskElement
import MetaData.TSCCMeta
import SZUtils.readTable

object UpSnapATAC2OfsnATAC {
  val workd = s"${TSCCMeta.projd}/05.CRE"
  val logd = s"${workd}/log/snATACsa2sample"
  val flagd = s"${workd}/flag/snATACsa2sample"
  val snATACsamplefnm = s"${workd}/src/main/resource/snATACsample.txt"

  val script = s"${workd}/src/main/python/sup.get.sa2v26.snATAC.py"
  val python = "/home/szu/miniforge3/envs/sa2/bin/python"

  val samples: List[String] =
    os.read.lines
    .stream(os.Path(snATACsamplefnm))
    .slice(0, Int.MaxValue)
    .toList

  class GenSnapATAC2AnnData(val sample: String) extends TaskElement {
    val flagfnm = s"${flagd}/${sample}.done"
    val logfnm = s"${logd}/${sample}.log"
    val skip = true
    def runCore(): Unit = {
      println(s"start GenSnapATAC2AnnData: ${sample}.")
      os.proc(
        python,
        script,
        sample).call(
        check = false, stdout = os.Path(logfnm), stderr = os.Path(logfnm))
      println(s"finish GenSnapATAC2AnnData: ${sample}.")
    }
  }

  def main(args: Array[String]) = {
    val sampleGenAnn = samples.map(s => new GenSnapATAC2AnnData(s))
    val forkJoinPool = new java.util.concurrent.ForkJoinPool(5)
    val parSampleGenAnn = sampleGenAnn.par
    parSampleGenAnn.tasksupport = new ForkJoinTaskSupport(forkJoinPool)
    parSampleGenAnn.foreach{
      t => t.run()
    }
  }
}
