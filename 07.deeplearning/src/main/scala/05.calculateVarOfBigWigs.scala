import os._
import MetaData.TSCCMeta
import SZUtils.readTable
import MetaData.readSubclassMeta
import PairedTagBigWig.{CleanDuplicatedPileups, BamToBigWig}
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*
import Bed.BedElement4
import Bed.BedElement4.readBed4
import SimpleBigWig.BedGraphElement
import SimpleBigWig.BedGraph.{toBigWig, toBedGraph}
import SimpleBigWig.mapBigWigOnRegionOneChr
import SimpleBigWig.loadChromosomeData
import SZUtils.writeStrings2File
import SimpleBigWig.mapBigWigOnRegion
import Genome.MouseGenome.ordChrs
import GRange.{GenomicRange, mouseGenomicRangeOrd}

object VarOfBigWig4DL {

  def decomposefnm(s: String): (String, String, String) = {
    val sc    = s.split("/").takeRight(2).head
    val h_rep = s.split("/").last.split("\\.")(0)
    val h     = h_rep.split("-")(0)
    val rep   = h_rep.split("-")(1)
    (sc, h, rep)
  }

  def main(args: Array[String]) = {
    val projd = TSCCMeta.projd
    val genomicBin4Varfnm =
      s"$projd/data/dl/intervals_524kb.bed"
    // 1. select subclasses based on #cells.
    // - mainly based on top number of cells
    //   will cover some major types
    // val scMeta =
    //   readSubclassMeta.sortBy(x => -(x.all.male + x.all.female))

    // NN subcleass have no replicated files.
    // val scs =
    //   Vector("319_Astro_TE_NN", "318_Astro_NT_NN", "327_Oligo_NN",
    //       "326_OPC_NN", "334_Microglia_NN", "053_Sst_Gaba",
    //       "061_STR_D1_Gaba", "062_STR_D2_Gaba", "037_DG_Glut",
    //       "016_CA1_ProS_Glut", "017_CA3_Glut", "008_L2_3_IT_ENT_Glut",
    //       "007_L2_3_IT_CTX_Glut", "006_L4_5_IT_CTX_Glut",
    //       "142_HY_Gnrh1_Glut", "135_STN_PSTN_Pitx2_Glut")

    val scs = Vector("053_Sst_Gaba",
          "061_STR_D1_Gaba", "062_STR_D2_Gaba", "037_DG_Glut",
          "016_CA1_ProS_Glut", "017_CA3_Glut", "008_L2_3_IT_ENT_Glut",
          "007_L2_3_IT_CTX_Glut", "006_L4_5_IT_CTX_Glut",
          "142_HY_Gnrh1_Glut", "135_STN_PSTN_Pitx2_Glut")
    // 2. then generate the bigwigs for each replicates
    //    - [female|male][A|B]
    val workd     = s"$projd/07.deeplearning"
    val logd      = s"$workd/log"
    val flagd     = s"$workd/flag"
    val rawbamd   = s"$projd/data/ptDNAbam/bam"
    val outbamd   = s"$workd/out/bam"
    val bedgraphd = s"$workd/out/bedgraph"
    val bigwigd   = s"$workd/out/bigwig"

    val allRawbams: Vector[String] =
      scs
      .map(sc => {
        TSCCMeta.rep
          .map(r => {
            TSCCMeta.modality.map(h => s"$rawbamd/$sc/$h-$r.srt.bam")
          })
          .flatten
      })
      .flatten
      .toVector

    val cleanDupPileupTasks = allRawbams.map(b => {
        val (sc, h, rep) = decomposefnm(b)
        new CleanDuplicatedPileups(
            inbam = b,
            outbam = s"$rawbamd/$sc/$h-$rep.cleanpileup.srt.bam",
            logfnm = s"$logd/${sc}_${h}_${rep}.cleanpileup.log",
            flagfnm = s"${flagd}/${sc}_${h}_${rep}.cleanpipleup.done",
            cutoff = 10,
            skip = true,
            check = true
        )
      })
    cleanDupPileupTasks.par.foreach(_.run())

    val bam2bwTasks: Vector[BamToBigWig] =
      allRawbams.map(b => {
        val (sc, h, rep) = decomposefnm(b)
        val upbam        = s"$rawbamd/$sc/$h-$rep.cleanpileup.srt.bam"
        new BamToBigWig(
          inbam = upbam,
          bw = s"$bigwigd/${sc}_${h}_${rep}.bw",
          logfnm = s"$logd/${sc}_${h}_${rep}.bigwig.log",
          flagfnm = s"$flagd/${sc}_${h}_${rep}.bigwig.done",
          skip = true,
          check = true,
          skipbw = true
        )
      })
    bam2bwTasks.par.foreach(_.run())

    // 3. load the intervals for variance calculations.
    val intervals: Map[String, List[BedElement4]] =
      readBed4(genomicBin4Varfnm)
        .groupBy(_.name)
        .map((k, v) => (k, v.sortBy(_.x)(using mouseGenomicRangeOrd)))

    // 4. calculate the variances
    // For each train / test / valid groups of intervals,
    // we calculate, each subclass and histone modification, the bigwig scores
    // along the intervals.
    val chrom2size = readTable(
      fnm = s"$projd/meta/mm10.chrom.sizes.lite",
      sep = "\t", head = false
    )

    val values: Map[String, Map[(String, String), Vector[Vector[Double]]]] =
      intervals.map((mltag, bedElmns) => {
        val sch2score = scs.par.map(sc => {
          TSCCMeta.modality.map(h => {
            val scores = TSCCMeta.rep.toVector.par.map(rep => {
              val bw =
                loadChromosomeData(
                  bigWigFilePath = s"$bigwigd/${sc}_${h}_${rep}.bw",
                  chromosomes = chrom2size.map(x => x(0)))
                  .map((k, v) => (k, v.toVector))
              val r = mapBigWigOnRegion(
                bw = bw,
                region = bedElmns.map(_.x),
                statistic = "mean",
                emptyValue = 0.0,
                keepOrder = true)
                .map(_.s)
              println(s"finish mapping for $sc $h $rep for $mltag")
              r
            })
            val r = scores.toVector.transpose
            ((sc, h), r)
          }) // end of modality
        }).toList
        val r = sch2score.flatten.toMap
        (mltag, r)
      })

    // subclass, histone, train / test / valid
    val sds: Map[(String, String, String), Vector[Double]] =
      values.flatMap((k, v) => {
        v.map((kk, vv) => {
          val intervalSds: Vector[Double] = vv.map(v3 => {
            val m: Double = v3.sum / v3.length
            val v: Double = v3.map(
              x => math.pow(x - m, 2)).sum / (v3.length.toDouble - 1.0)
            math.sqrt(v)
          })
          ((kk._1, kk._2, k), intervalSds)
        })
      })

    // write sds as bed graphs
    sds.foreach((k, v) => {
      val (sc, h, mltag) = k
      val bgElmns = intervals(mltag)
        .zip(v)
        .map((be, s) => {
          BedGraphElement(g = be.x, s = s)
        })
      writeStrings2File(
          content = bgElmns.map(x => x.mkString("\t")),
          to = s"$bedgraphd/${sc}_${h}_${mltag}.std.acrossSampleRep.bdg",
          overwrite = true,
          head = ""
      )
    })

    // 5. output as bigwig files.
    val chromSizef = s"$projd/meta/mm10.chrom.sizes.lite"
    val bg2bw = "/home/szu/mambaforge/" +
      "envs/deeptools/bin/" + "bedGraphToBigWig"

    scs.foreach(sc => {
      TSCCMeta.modality.foreach(h => {
        intervals.keys.foreach(mltag => {
          val mykey = (sc, h, mltag)
          val bgf = s"$bedgraphd/${sc}_${h}_${mltag}.std.acrossSampleRep.bdg"
          val bwf = s"$bigwigd/${sc}_${h}_${mltag}.std.acrossSampleRep.bw"
          toBigWig(bg2bw = bg2bw, chromSizef = chromSizef, bgf = bgf,
              bwf = bwf)
        }) // end of intervals' keys foreach
      })   // end of modality foreach
    })     // end of scs foreach
  }        // end of main
}          // end of object VarOfBigWig
