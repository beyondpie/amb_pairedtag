import os._
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*
import Bam.getNRead
import MetaData.TSCCMeta
import SZUtils.writeStrings2File

object StatRead {
  def main(args: Array[String]) = {
  val projd = TSCCMeta.projd
  val bamd = s"${projd}/data/ptDNAbam/bam"
  val groups =
    os.list(os.Path(bamd)).map(x => x.toString.split("/").last)
  val r =
    groups
      .map(
        g => (g,
          os.list(os.Path(bamd) / g).map(x => x.toString.split("/").last)))
      .par
      .map(
        (g, names) => {
          println(g)
          val rr = names.par.map(y => (g, y,
            getNRead(List(bamd, g, y).mkString("/")))).toList
          println(s"done: ${g}.")
          rr
        }
      )
      .toList
  val rr =
    r.map(x => x.map(y => y.toList.mkString(",")))
      .flatten
  val outfnm = s"${projd}/04.peakcalling/src/main/resource/bam2nread.full.csv"
  writeStrings2File(rr, outfnm, overwrite = true)
    
  }
}

