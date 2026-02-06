import os._
import SZUtils.{readHead, readTable}
import SZUtils.writeStrings2File
import SZUtils.colorHex2RGB
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*

object GenTransitionMap {
  val projd = "/projects/ps-renlab2/szu/projects/amb-pairedtag"
  val workd = s"${projd}/06.ChromHMM"
  val stated = s"${workd}/out"
  val outd = s"${workd}/out/updateTransitionMap"

  val chrs = 1.to(19).map(i => s"chr${i}")

  def getTransitionCount(
    states: List[String]): Map[(String, String), Int] = ???

  def main(args: Array[String]) = {
    val scs: List[String] = ???
    val states: Map[String, Map[String, List[String]]] = ???
    // val transitions: Map[String, Map[String, List[String]]] =
    //   scs.foreach{ sc => {
    //     ???
    //   }}
  }
}
