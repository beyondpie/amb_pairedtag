import os._
import scala.collection.mutable.ListBuffer
import GRange.{GenomicRange, mouseChrOrd}
import bioscala.LightCoord.findOvlpOneChrSorted
import CEMBATAC._
import PairedTagBigWig._

object TestGenomicRange {

  def main(args: Array[String]) = {
    val sc   = "334_Microglia_NN"
    val chrs = List("chr1", "chr2")
    val bw   = loadptRNAbw(sc, chrs)
    val CREs = loadscATAC(sc).map(_.x)

    // test ovlp in one region
    val subject =
      bw("chr1")
      .sortBy(_.g)
      .map(x => (x.g.startFrom, x.g.endTo))
      .toVector
    val query =
      CREs
      .groupBy(_.chrom)("chr1")
      .sorted
      .map(x => (x.startFrom, x.endTo))
      .toVector
    val ovlp = findOvlpOneChrSorted(query = query, subject = subject)

    // val query   = Vector((0, 100), (200, 300), (500, 600))
    // val subject = Vector((0, 10), (20, 100), (550, 600))
    // val si = 0
    // val q = (query(0), 0)
    // val r = ListBuffer.empty[(Int, Vector[Int])]
    // query.zipWithIndex.foldLeft[Int](0)((si, q) => {
    //   val sovlp =
    //     subject.zipWithIndex
    //       .dropWhile((x, id) => {
    //         (id < si) || (x._2 < q._1._1) || (x._1 >= q._1._2)
    //       })
    //       .takeWhile((x, id) => {
    //         (x._2 >= q._1._1) && (x._1 < q._1._2)
    //       })
    //   if (sovlp.length >= 1) {
    //     r.addOne((q._2, sovlp.map(_._2)))
    //     sovlp.head._2
    //   } else {
    //     si
    //   }
    // })

  } // end of main
}
