import os._
import org.broad.igv.bbfile.BBFileReader
import org.broad.igv.bbfile.WigItem
import java.io.IOException
import java.util.Iterator
import scala.collection.mutable.{ListBuffer, Map}
import SimpleBigWig._
import GRange.{GenomicRange, mouseChrOrd}
import bioscala.LightCoord.findOvlpOneChrSorted
import CEMBATAC._
import PairedTagBigWig._

object TestBigWig {
  def testRead():Unit = {
    // * test reading bigwig
    val projd      = "/projects/ps-renlab2/szu/projects/amb_pairedtag"
    val bwd        = s"${projd}/data/ptDNAbam/bigwig"
    val sc         = "001_CLA_EPd_CTX_Car3_Glut"
    val h          = "H3K4me1"
    val bgf        = s"${bwd}/${sc}.${h}.e100.bs100.sm300.bw"
    val read       = new BBFileReader(bgf)
    val chromosome = "chr2"
    val start      = 100000
    val endTo      = 4782800
    val wigIterator =
      read.getBigWigIterator(chromosome, 0, chromosome, endTo, true)
    while (wigIterator.hasNext()) {
      val item = wigIterator.next();
      println(List(item.getStartBase().toString,
        item.getEndBase().toString, item.getWigValue().toString))
    }
  } // end of fn testRead

  def testScore(): Unit = {
    val sc   = "334_Microglia_NN"
    val chrs = List("chr1", "chr2")
    val bw   = loadptRNAbw(sc, chrs)
    val CREs = loadscATAC(sc).map(_.x)
    val subject =
      bw("chr1")
      .sortBy(_.g)
      .toVector
    val query =
      CREs
      .groupBy(_.chrom)("chr1")
      .sorted
      .toVector
    val scores = mapBigWigOnRegionOneChr(bw = subject, region = query)


  }
    
}


