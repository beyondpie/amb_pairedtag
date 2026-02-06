import os._
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*
import MetaData.TSCCMeta
import SZUtils.TaskElement
import MetaData.TSCCTools
import MergePeak.readReproPeak
import Peak.Peaks
import Bed.GFF
import SZUtils.writeStrings2File
import scala.compiletime.ops.boolean

object RunROSE {
  val projd: String = TSCCMeta.projd
  val bamd: String = s"${projd}/data/ptDNAbam/bam"
  val peakd: String = s"${projd}/data/pairedtag_peak/subclass_peak"
  val workd: String = s"${projd}/10.superEnhancer"
  val outd: String = s"${workd}/out/ROSE"
  val logd: String = s"${workd}/log/ROSE"
  val python: String = "/home/szu/mambaforge/bin/python"
  val rscript: String = "/home/szu/mambaforge/envs/seurat/bin/Rscript"
  val ROSE_stitchCRE: String =
    s"${workd}/src/main/python/ROSE_stitchCRE.py"
  val ROSE_bam2GFF: String = s"${workd}/src/main/python/ROSE_bamToGFF.py"
  val ROSE_calDensity: String =
    s"${workd}/src/main/python/ROSE_calDensity.py"
  val ROSE_callSuper: String = s"${workd}/src/main/R/ROSE_callSuper.R"
  val ROSE_geneMapper: String =
    s"${workd}/src/main/python/ROSE_geneMapper.py"
  val annotfnm: String = s"${workd}/src/main/resource/mm10_refseq.ucsc"

  class SubclassROSE(
    val sc: String,
    val stitch: Int = 12500,
    val tssWindow: Int = 2000,
    val bam: String,
    val peakfnm: String
  ) {
    val logfnm = s"${logd}/${sc}_runROSE.log"
    val peakGFF = s"${outd}/${sc}/${sc}.H3K27ac.peak.gff"
    val stitchedName = s"${sc}_${stitch.toFloat / 1000.0}KB_STITCHED"
    val stitchedGFF = s"${outd}/${sc}/${stitchedName}.gff"
    val stitchedMapName = s"${stitchedName}_TSS_DISTAL"
    val bamName = os.Path(bam).last
    val bamStitchedGFF =
      s"${outd}/${sc}/${stitchedMapName}_${bamName}_MAPPED.gff"
    val bamPeakGFF = s"${outd}/${sc}/${stitchedName}_${bamName}_MAPPED.gff"
    val densityfnm = s"${outd}/${sc}/${stitchedName}_REGION_MAP.txt"
    val superTableFile = s"${outd}/${sc}/${sc}_SuperStitched.table.txt"
    val allTableFile = s"${outd}/${sc}/${sc}_AllStitched.table.txt"

    if (!os.exists(os.Path(s"${outd}/${sc}"))) {
      os.makeDir.all(os.Path(s"${outd}/${sc}"))
    }

    // 1. prepare the GFF of peaks
    def transformPeaks2GFF(skip: Boolean = true): Unit = {
      if (!os.exists(os.Path(peakGFF)) || !skip) {
        val ps: Peaks = readReproPeak(peakfnm, head = true)
        val r = ps
          .map(x => x.toGFF(fea = "CRE", frame = "."))
          .map(x => x.toSingleStringLine(sep = "\t"))
        writeStrings2File(r, to = peakGFF, overwrite = true, head = "")
      }
    } // end of transformPeak2GFF

    // 2. put the major part in Python, and run through command line.
    def stitchRegions(): Unit = {
      os.proc(
        python,
        ROSE_stitchCRE,
        peakGFF,
        stitch.toString,
        tssWindow.toString,
        stitchedGFF,
        annotfnm
      ).call(check = true, stderr = os.Path(logfnm), stdout = os.Path(logfnm))
    }

    // 3. bam2GFF
    def mapBam2GFF(stitched: Boolean = true): Unit = {
      val fromGFF = if (stitched) {
        stitchedGFF
      } else { peakGFF }
      val outGFF = if (stitched) {
        bamStitchedGFF
      } else { bamPeakGFF }
      os.proc(
        python,
        ROSE_bam2GFF,
        "-f",
        "1",
        "-e",
        "200",
        "-r",
        "-m",
        "1",
        "-b",
        bam,
        "-i",
        fromGFF,
        "-o",
        outGFF
      ).call(check = true, stderr = os.Path(logfnm), stdout = os.Path(logfnm))
    }

    // 4. calculate density
    def calDensity(): Unit = {
      os.proc(
        python,
        ROSE_calDensity,
        stitchedGFF,
        peakGFF,
        bam,
        bamPeakGFF,
        bamStitchedGFF,
        densityfnm,
        stitchedMapName
      ).call(check = true, stderr = os.Path(logfnm), stdout = os.Path(logfnm))
    }

    // 5. rank Super
    def rankSuper(): Unit = {
      os.proc(
        rscript,
        ROSE_callSuper,
        s"${outd}/${sc}/",
        densityfnm,
        sc,
        "NONE"
      ).call(check = true, stderr = os.Path(logfnm), stdout = os.Path(logfnm))
    }

    // 6. map Gene
    def mapGene(): Unit = {
      List(superTableFile, allTableFile).foreach { x =>
        {
          os.proc(python, ROSE_geneMapper, "-g", "mm10", "-i", x, "-r", "TRUE")
            .call(
              check = true,
              stderr = os.Path(logfnm),
              stdout = os.Path(logfnm)
            )
        }
      }
    } // end of function mapGene

  } // end of class

  def main(args: Array[String]) = {
    List(outd, logd).map(d => os.Path(d)).foreach { d =>
      if (!os.exists(d)) {
        os.makeDir.all(d)
      }
    }

    val scs: List[String] =
      os.list(os.Path(bamd))
        .map(d => d.baseName)
        .filter(d => os.exists(os.Path(s"${bamd}/${d}/H3K27ac.srt.bam")))
        .filter(
          d => os.exists(os.Path(s"${peakd}/${d}-H3K27ac.BestSPM.peak")))
        .filter(sc => sc != "327_Oligo_NN")
        .toList

    // val scs = List("326_OPC_NN", "327_Oligo_NN")

    val sc2bam =
      scs.map(sc => (sc, os.Path(s"${bamd}/${sc}/H3K27ac.srt.bam"))).toMap

    val sc2peak =
      scs.map(sc => (
        sc, os.Path(s"${peakd}/${sc}-H3K27ac.BestSPM.peak"))).toMap

    val sc2rose: Map[String, SubclassROSE] =
      scs
        .map(sc =>
          (
            sc,
            new SubclassROSE(
              sc = sc,
              bam = sc2bam(sc).toString,
              peakfnm = sc2peak(sc).toString
            )
          )
        )
        .toMap

    // DEBUG
    // val sc = "001_CLA_EPd_CTX_Car3_Glut"
    // val stitch = 12500
    // val tssWindow = 2000
    // val bam = s"${bamd}/${sc}/H3K27ac.srt.bam"
    // val peakfnm = s"${peakd}/${sc}-H3K27ac.BestSPM.peak"

    // transformPeaks2GFF(skip = true)
    // stitchRegions()
    // mapBam2GFF(stitched = true)
    // mapBam2GFF(stitched = false)
    // calDensity()
    // rankSuper()

    // val t = sc2rose("137_PH_an_Pitx2_Glut")

    //  runROSE
    scs.par.foreach(sc => {
      val t: SubclassROSE = sc2rose(sc)
      if (!os.exists(os.Path(t.allTableFile))) {
        println(s"1. Transform Peak2GFF for ${sc}.")
        t.transformPeaks2GFF(skip = true)
        println(s"2. Stitch peak regions for ${sc}.")
        t.stitchRegions()
        println(s"3. Map bam on stitichedGFF for ${sc}.")
        t.mapBam2GFF(stitched = true)
        println(s"3. Map bam on peakGFF for ${sc}.")
        t.mapBam2GFF(stitched = false)
        println(s"4. Calculate density for ${sc}.")
        t.calDensity()
        println(s"5. Rank super enhancers for ${sc}.")
        t.rankSuper()
        println(s"6. Map genes of super enhancers for ${sc}.")
        t.mapGene()
        println(s"ROSE for ${sc} is done. Good luck!")
      } else {
        println(s"${t.allTableFile} exists, skip it.")
      }
    }) // end of foreach
  } // end of main
} // end of RunROSE
