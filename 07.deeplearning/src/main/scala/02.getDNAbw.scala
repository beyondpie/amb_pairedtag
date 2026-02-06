import os._
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*
import scala.collection.mutable.LinkedHashMap

import java.io.File
import htsjdk.samtools.SamReader
import htsjdk.samtools.SamReaderFactory
import htsjdk.samtools.SAMRecord

import SZUtils.TaskElement
import MetaData.TSCCMeta.projd
import MetaData.TSCCMeta.modality

import Bam.getPileUpMap
import Bam.removeDuplicatedPileUp
import PairedTagBigWig.{CleanDuplicatedPileups, BamToBigWig}

object GetPairedTagDNABigWig4DeepLearningAndChromHMM {
  val cutoff: Int = 10
  val bamd = s"${projd}/data/ptDNAbam/bam"
  val bwd = s"${projd}/data/ptDNAbam/bigwig"
  val logd = s"${projd}/data/ptDNAbam/log"
  val flagd = s"${projd}/data/ptDNAbam/tag"
  val groups =
    os.list(os.Path(bamd))
      .map(x => x.baseName)
      .toList
  // sc2bamfnm
  val g2bamfnms: Map[String, Array[String]] =
    groups
      .map(g => {
        val fs = modality
          .map(m => s"${bamd}/${g}/${m}.srt.bam")
          .filter(f => os.exists(os.Path(f)))
        (g, fs.toArray)
      })
      .toMap

  val g2cleanBams =
    g2bamfnms
      .map((g, fs) =>
        (
          g,
          fs.map(f => {
            val m = f.split("\\/").last.split("\\.").head
            s"${bamd}/${g}/${m}.cleanpileup.srt.bam"
          }).toArray
        )
      )
      .toMap

  val g2cleanBamLogs =
    g2bamfnms
      .map((g, fs) =>
        (
          g,
          fs.map(f => {
            val m = f.split("\\/").last.split("\\.").head
            s"${logd}/${g}.${m}.cleanpileup.log"
          }).toArray
        )
      )
      .toMap

  val g2cleanBamFlags =
    g2bamfnms
      .map((g, fs) =>
        (
          g,
          fs.map(f => {
            val m = f.split("\\/").last.split("\\.").head
            s"${flagd}/${g}.${m}.cleanpileup.done"
          }).toArray
        )
      )
      .toMap

  val g2bw =
    g2bamfnms
      .map((g, fs) =>
        (
          g,
          fs.map(f => {
            val m = f.split("\\/").last.split("\\.").head
            val sm: Int = if (m == "H3K9me3") {
              1000
            } else {
              300
            }
            s"${bwd}/${g}.${m}.e100.bs100.sm${sm}.bw"
          }).toArray
        )
      )
      .toMap

  val g2bwlogs =
    g2bamfnms
      .map((g, fs) =>
        (
          g,
          fs.map(f => {
            val m = f.split("\\/").last.split("\\.").head
            val sm: Int = if (m == "H3K9me3") {
              1000
            } else {
              300
            }
            s"${logd}/${g}.${m}.e100.bs100.sm${sm}.log"
          }).toArray
        )
      )
      .toMap

  val g2bwflags =
    g2bamfnms
      .map((g, fs) =>
        (
          g,
          fs.map(f => {
            val m = f.split("\\/").last.split("\\.").head
            val sm: Int = if (m == "H3K9me3") {
              1000
            } else {
              300
            }
            s"${flagd}/${g}.${m}.e100.bs100.sm${sm}.done"
          }).toArray
        )
      )
      .toMap

  def main(args: Array[String]) = {
    // val g2cleanPileupTasks = g2bamfnms.map{
    //   (g, fs) => {
    //     (g, fs.zipWithIndex.map(
    //       (f, i) => new CleanDuplicatedPileups(
    //         inbam = f,
    //         outbam = g2cleanBams(g)(i),
    //         logfnm = g2cleanBamLogs(g)(i),
    //         flagfnm = g2cleanBamFlags(g)(i),
    //         skip = true,
    //         check = false
    //       )
    //     ))
    //   }
    // }

    // g2cleanPileupTasks.par.foreach(
    //   (g, ts) => {
    //     if(ts.length > 0) {
    //       ts.par.foreach{
    //         t => t.run()
    //       }
    //     }
    //   }
    // )

    // bigwig
    val g2bwTasks = g2bw.map((g, fs) => {
      (
        g,
        fs.zipWithIndex.map((f, i) =>
          new BamToBigWig(
            inbam = g2cleanBams(g)(i),
            bw = f,
            logfnm = g2bwlogs(g)(i),
            flagfnm = g2bwflags(g)(i),
            skip = false,
            check = false,
            skipbw = true
          )
        )
      )
    }) // end of map to define g2bwTasks

    g2bwTasks.par.foreach { (g, ts) =>
      {
        if (ts.length > 0) {
          ts.par.foreach { t =>
            t.run()
          }
        }
      } // end of defination of fn for foreach
    } // end of foreach of g2bwTask

    // val ts = g2bwTasks("1193_Endo_NN_1")(0)

  } // end of main
} // end of object
