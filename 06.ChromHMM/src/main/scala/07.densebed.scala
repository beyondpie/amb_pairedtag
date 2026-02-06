import os._
import SZUtils.{readHead, readTable}
import SZUtils.writeStrings2File
import SZUtils.colorHex2RGB
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*

object GenChromHMMDenseBed {
  val projd         = "/projects/ps-renlab2/szu/projects/amb_pairedtag"
  val workd         = s"${projd}/06.ChromHMM"
  val rawDenseBedir = s"${workd}/out/model_bypeak_b200_s18"
  val outd          = s"${workd}/out/updateDenseBed"
  val scs =
    os.list(os.Path(rawDenseBedir))
      .map(x => x.toString)
      .filter(x => x.contains("_dense.bed"))
      .map(x => x.split("/").toList.last.replace("_dense.bed", ""))
      .toList

  val stateMapf = s"${projd}/meta/chromHMM18statesAnnotation.csv"
  val stateMap  = readTable(stateMapf, sep = ",", head = true)
  val state2group: Map[String, String] =
    stateMap
    .map(x => (x(0), x(2)))
    .toMap

  val group2RGB: Map[String, String] =
    stateMap
    .map(x => (x(2), colorHex2RGB(x(4))))
    .toMap

  def genDenseBed(sc: String): Unit = {
    val rawf = s"${rawDenseBedir}/${sc}_dense.bed"
    val outf = s"${outd}/${sc}_dense.bed"
    val head = readHead(rawf)
    val content = readTable(rawf, sep = "\t", head = true)
      .map(x => {
        val chr       = x(0)
        val startFrom = x(1)
        val endTo     = x(2)
        val group     = state2group(x(3))
        val score     = "0"
        val strand    = "."
        val RGB       = group2RGB(group)
        List(chr, startFrom, endTo, group, score, strand, startFrom,
            endTo, RGB)
      })
      .map(x => x.mkString("\t"))
    writeStrings2File(
        content = content,
        to = outf,
        overwrite = true,
        head = head
    )
  }

  def main(args: Array[String]) = {
    scs.par.foreach { sc =>
      {
        println(s"Gen dense bed of ${sc}...")
        genDenseBed(sc)
      }
    }
  } // end of main
  
}
