import os.*

import scala.util.Random
import scala.collection.parallel.CollectionConverters.*
import MetaData.{PairedTagBarcodeMeta, toStringOfAllenIdLabel}
import MetaData.PairedTagBarcodeMeta.readPairedTagCellMetaFile
import MetaData.{TSCCMeta, TSCCTools}
import SZUtils.{writeStrings2File}
import MetaData.TSCCMeta.{ptclusterDNAbamd, ptscDNAbamd}

/*
 using DG Glut neurons, and subset 100,000 (344 subclass/mod/sexrep filtered out, out of 2,376 total)
 ; 200,000 (551 filtered out); 300,000 (706 filtered out);
 500,000 (908 filtered out); 750,000 (1082 filtered out);
 and 1,000,000 (1,234 filtered out) reads, to call peaks for each modality.

 1. Generate the merged bam files with downsampling for DG Glut
 - different modalities independently
 2. call peaks using macs3 with fixed parameters
    - generate the tracks
 3. saving peaks then statistic the peaks callback using overlapping
 (let's use R for the overlapping analysis)
 4. visualize the peaks with bigwig tracks

 */

object SaturationAnalysis {
// meta data
  val workd = s"${TSCCMeta.projd}/04.peakcalling"
  lazy val cellMeta = readPairedTagCellMetaFile()
  lazy val barcode2Bam = TSCCMeta.getBarcode2Bam(cellMeta)
  lazy val NeuscModBamStat = TSCCMeta.getRawBamStat(
    TSCCMeta.NeuscSexRepBamStatFile
  )
  lazy val NNspModBamStat = TSCCMeta.getRawBamStat(
    TSCCMeta.NNspSexRepBamStatFile
  )

  val nReads: List[Int] = List(
    100_000, 200_000, 500_000, 750_000, 1_000_000, 1_250_000, 1_500_000,
    1_750_000, 2_000_000, 2_500_000, 4_000_000, 5_000_000, 10_000_000,
    15_000_000, 20_000_000
  )

  def downsample[T](c: List[T], n: Int): List[T] = {
    Random.shuffle(c).take(n)
  }

  def sampleBarcode[T](ns: List[Int])(from: List[T]) = {
    ns.map(n => downsample(from, n))
  }

  class SaturationAnalysis(pattern: String, isNeu: Boolean = true) {
    val outd = s"${workd}/out/saturation/${pattern}"
    if (!os.exists(os.Path(outd))) {
      os.makeDir(os.Path(outd))
    }
    println(s"prepare results to ${outd}.")

    val mod2BarcodeMeta =
      cellMeta
        .filter(_.annotQuality == "Good")
        .filter(x =>
          if (isNeu) {
            toStringOfAllenIdLabel(x.annot.sc.get).contains(pattern)
          } else {
            toStringOfAllenIdLabel(x.annot.sp.get).contains(pattern)
          }
        )
        .groupBy(_.exp.modularity)
    val bamStat =
      if (isNeu) {
        NeuscModBamStat.filter((k, v) =>
          (toStringOfAllenIdLabel(k._1).contains(pattern))
        )
      } else {
        NNspModBamStat.filter((k, v) =>
          (toStringOfAllenIdLabel(k._1).contains(pattern))
        )
      }

    val mod2NdsCells: Map[String, List[Int]] =
      bamStat
        .map((k, v) =>
          (
            k,
            nReads
              .map(x =>
                math
                  .ceil(x.toFloat / v._3)
                  .toInt
                  .min(mod2BarcodeMeta(k._2).length)
              )
              .distinct
          )
        )
        .map((k, v) => (k._2, v))
    val mod2DsCells =
      mod2BarcodeMeta.map((k, v) =>
        (k, sampleBarcode(mod2NdsCells(k))(v.map(_.barcode)))
      )
    val mod2DsFiles = mod2DsCells.map { (k, v) =>
      (
        k,
        v.map(x =>
          x.length ->
            s"${outd}/${pattern}_${k}_${x.length}.listofbam.txt"
        ).toMap
      )
    }
    val mod2DsMergedBam: Map[String, Map[Int, String]] =
      mod2NdsCells.map((k, v) =>
        (
          k,
          v.map(x =>
            (x ->
              s"${outd}/${pattern}_${k}_${x}.srt.bam")
          ).toMap
        )
      )

    def updateListOfBams: Unit = {
      mod2DsCells
        .map((k, v) => (k, v.map(b => b.map(bb => barcode2Bam(bb)))))
        .foreach { (modality, listOfListOfFiles) =>
          {
            listOfListOfFiles.foreach { x =>
              writeStrings2File(x, mod2DsFiles(modality)(x.length))
            }
          }
        }
    }

    def mergeBamsInPar: Unit = {
      mod2DsFiles.par.foreach { (k, v) =>
        {
          v.foreach((i, j) => {
            TSCCTools.mergeBamCore(
              from = mod2DsFiles(k)(i),
              to = mod2DsMergedBam(k)(i),
              filterMAPQ = true,
              MAPQ = 10
            )
          })
        }
      }
    }

    def runMACS3InPar: Unit = {
      mod2DsMergedBam.par.foreach { (k, v) =>
        {
          v.foreach { (size, bamfnm) =>
            {
              val name = s"${pattern}_${k}_ds${size}"
              val logfnm = s"${outd}/${name}_macs3.log"
              val broad = if (k == "H3K27ac") false else true
              TSCCTools.callPeak(
                macs3 = TSCCTools.macs3,
                bamfnm,
                name = name,
                shift = "-100",
                extsize = "200",
                genome = "mm",
                qvalue = "0.05",
                outdir = outd,
                logfnm = logfnm,
                broad = broad,
                broad_cutoff = "0.1",
                cutoff_analysis = false
              )
            }
          }
        }
      }
    }
  } // end of class SaturationAnalysis

  def main(args: Array[String]) = {
    val OPCNN1 = new SaturationAnalysis(pattern = "OPC_NN_1", isNeu = false)
    OPCNN1.updateListOfBams
    OPCNN1.mergeBamsInPar
    OPCNN1.runMACS3InPar

    OPCNN1.mod2DsMergedBam.par.foreach { (k, v) =>
      {
        v.foreach { (size, bamfnm) =>
          {
            val outd = OPCNN1.outd
            val pattern = "OPC_NN_1"
            val name = s"${pattern}_${k}_ds${size}"
            val logfnm = s"${outd}/${name}_macs3.log"
            val broad = if (k == "H3K27ac") false else true
            TSCCTools.callPeak(
              TSCCTools.macs3,
              bamfnm,
              name = name,
              shift = "-100",
              extsize = "200",
              genome = "mm",
              qvalue = "0.05",
              outdir = outd,
              logfnm = logfnm,
              broad = broad,
              broad_cutoff = "0.1",
              cutoff_analysis = false
            )
          }
        }
      }
    }

    val SstChodl =
      new SaturationAnalysis(pattern = "Sst_Chodl_Gaba", isNeu = true)
    SstChodl.updateListOfBams
    SstChodl.mergeBamsInPar
    SstChodl.runMACS3InPar

    SstChodl.mod2DsMergedBam.par.foreach { (k, v) =>
      {
        v.foreach { (size, bamfnm) =>
          {
            val outd = SstChodl.outd
            val pattern = "Sst_Chodl_Gaba"
            val name = s"${pattern}_${k}_ds${size}"
            val logfnm = s"${outd}/${name}_macs3.log"
            val broad = if (k == "H3K27ac") false else true
            TSCCTools.callPeak(
              macs3 = TSCCTools.macs3,
              bamfnm,
              name = name,
              shift = "-100",
              extsize = "200",
              genome = "mm",
              qvalue = "0.05",
              outdir = outd,
              logfnm = logfnm,
              broad = broad,
              broad_cutoff = "0.1",
              cutoff_analysis = false
            )
          }
        }
      }
    }

  }
}
