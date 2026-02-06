import os._
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*
import SZUtils.SimpleCommandTask
import PairedTagChromHMM.{stateAnnotfnm, loadChromHMMStateAnnots}
import PairedTagChromHMM.chromHMMStates as ordAnnots
import MetaData.TSCCMeta.projd
import ChromHMM.{SegmentElement, loadSegmentBed}
import SZUtils.readTable
import SZUtils.writeStrings2File
import SimpleBigWig.BedGraph.{toBigWig, toBedGraph}
import SimpleBigWig.BedGraphElement
import GRange.GenomicRange
import SZUtils.TaskElement

object ChromHMMarginalProb2BigWig {
  val root        = "/projects/ps-renlab2/szu"
  val chromHMMjar = s"${root}/softwares/ChromHMM/ChromHMM.jar"
  val bg2bw = "/home/szu/mambaforge/envs/deeptools/bin/bedGraphToBigWig"
  val chromSizef = s"${projd}/meta/mm10.chrom.sizes.lite"
  val densebedir = s"${projd}/data/chromHMM/denseBed"
  val ordChrs    = 1.to(19).map(i => "chr" + i.toString)

  val modelfnm =
    s"${projd}/06.ChromHMM/out/model_bypeak_b200_s18/model_18.txt"

  val inputd = s"${projd}/06.ChromHMM/out/bBed_bypeak_b200"

  val outd        = s"${projd}/06.ChromHMM/out/segmentProb"
  val probd       = s"${outd}/POSTERIOR"
  val chrs        = 1.to(19).map(i => "chr" + i.toString)
  val stateAnnots = loadChromHMMStateAnnots(stateAnnotfnm)

  val stat2annot =
    stateAnnots.map(x => ("E" + x.stateid, x.groupName)).toMap

  // Index from 1 in ChromHMM
  val ordStates = 1.to(18).map(i => s"E${i}")
  // index start from 0
  val annot2Index = ordAnnots.zipWithIndex.toMap

  // for task
  val flagd = s"${projd}/06.ChromHMM/flag"
  val logd  = s"${projd}/06.ChromHMM/log"

  def calcMarginalProbTask(modelfnm: String, inputd: String,
    outd: String): Unit = {
    os.proc("java", "-mx800000M", "-jar", chromHMMjar,
        "MakeSegmentation", "-b", "200", "-printposterior",
        "-printstatebyline", modelfnm, inputd, outd)
      .call(check = true)
  }

  def checkStateOrdered(states: List[String]): Boolean = {
    // match index in ChromHMM
    states.zipWithIndex.forall((s, i) => s == s"E${i + 1}")
  }

  def loadStateProb(fnm: String): Vector[Vector[Double]] = {
    val r = readTable(fnm, sep = "\t", head = true).toVector
    if (checkStateOrdered(r(0))) {
      r.slice(from = 1, until = r.length)
        .map(x => x.map(y => y.toDouble).toVector)
        .toVector
    } else {
      println("state order not match, reorder.")
      val states   = r(0).zipWithIndex.toMap
      val ordIndex = ordStates.map(s => states(s))
      r.slice(from = 1, until = r.length)
        .map(x => {
          val raw = x.map(y => y.toDouble).toVector
          ordIndex.map(i => raw(i)).toVector
        })
        .toVector
    }
  }

  def mapStateProb2AnnotProb(x: Vector[Double]): Vector[Double] = {
    val r = x.zipWithIndex
      .map((p, i) => (stat2annot(s"E${i + 1}"), p))
      .groupBy(_._1)
      .map((a, s) => (a, s.map(_._2).sum.min(1.0)))
      .toMap
    ordAnnots.map(a => r(a)).toVector
  }

  def getStateAnnot2Prob(
    raws2p: Map[String, Vector[Vector[Double]]]) = {
    raws2p.map((chr, probs) =>
      (chr, probs.map(p => mapStateProb2AnnotProb(p))))
  }

  def toGenomeBin(i: Int): (Int, Int) = {
    (i * 200, (i + 1) * 200)
  }

  def write2bedGraph(annotProb: Map[String, Vector[Vector[Double]]],
    prefix: String): Unit = {
    ordAnnots.foreach { a =>
      {
        val chrs = ordChrs.filter(chr => annotProb.contains(chr))
        val r = chrs
          .map(chr =>
            annotProb(chr)
              .map(x => x(annot2Index(a)))
              .zipWithIndex
              .filter((x, i) => x > 0.0)
              .map((x, i) =>
                new BedGraphElement(
                    g = new GenomicRange(chr,
                        startFrom = toGenomeBin(i)._1,
                        endTo = toGenomeBin(i)._2),
                    s = x
                )))
          .flatten
          .toList
        toBedGraph(r, outf = s"${prefix}_${a}.bedgraph")
      }
    }
  }

  class getTask(val sc: String) extends TaskElement {
    val flagfnm = s"${flagd}/${sc}_prob2bigwig.done"
    val logfnm  = s"${flagd}/${sc}_prob2bigwig.log"
    val skip    = true
    def runCore(): Unit = {
      println(s"loading probs from ${sc}.")
      val rawStateProb =
        ordChrs.par
          .map(chr =>
            (chr,
                loadStateProb(
                    s"${probd}/${sc}_18_${chr}_posterior.txt")))
          .toList
          .toMap
      val annotProb = getStateAnnot2Prob(rawStateProb)
      write2bedGraph(annotProb, prefix = s"${outd}/bedgraph/${sc}")
      // save result per annots.
      println(s"save bigwig for ${sc}.")
      ordAnnots.par.foreach { a =>
        {
          toBigWig(
              bg2bw,
              chromSizef,
              bgf = s"${outd}/bedgraph/${sc}_${a}.bedgraph",
              bwf = s"${outd}/bigwig/${sc}_${a}.bw"
          )
        }
      }
    }
  }

  def main(args: Array[String]) = {
    // each one may need 5-10min.
    // should be run at very first.
    // calcMarginalProbTask(modelfnm, inputd, outd)

    val scs =
      os.list(os.Path(densebedir))
        .map(x => x.baseName.replace("_18_dense", ""))

    val tasks = scs.map(sc => new getTask(sc))
    tasks.foreach{t => t.run()}

  } // end of main

}
