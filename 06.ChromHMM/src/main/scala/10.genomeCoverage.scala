import os._
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*
import scala.util.Random
import GRange.GenomicRange
import SZUtils.readTable
import SZUtils.writeStrings2File

object GetGenomeCoverageHeatmap4ChromHMM {
  val root  = "/projects/ps-renlab2/szu"
  val projd = s"${root}/projects/amb_pairedtag"
  val workd = s"${projd}/06.ChromHMM"
  val ATACbwd = s"${projd}/data/snATAC/subclass_bigwig_bamCoverage"
  val Histonesbwd = s"${projd}/data/ptDNAbam/bigwig"

  val computeMatBin =
    "/home/szu/mambaforge/envs/deeptools/bin/computeMatrix"

  val plotHeatmapBin =
    "/home/szu/mambaforge/envs/deeptools/bin/plotHeatmap"

  val scs: List[String] = List("005_L5_IT_CTX_Glut", "319_Astro_TE_NN")
  val sc2bwf: Map[String, List[String]] = scs.map(
    sc => 
      (sc, List(
        s"${ATACbwd}/${sc}.ATAC.e100.bs100.sm300.bw",
        s"${Histonesbwd}/${sc}.H3K27ac.e100.bs100.sm300.bw",
        s"${Histonesbwd}/${sc}.H3K27me3.e100.bs100.sm300.bw",
        s"${Histonesbwd}/${sc}.H3K4me1.e100.bs100.sm300.bw",
        s"${Histonesbwd}/${sc}.H3K9me3.e100.bs100.sm1000.bw",
      ))
  ).toMap

  // val states: List[Int] = List(5, 10, 15, 18, 20, 23, 25)
  val states: List[Int] = List(18)

  def splitSegments(f: String, n: Int): Map[String, List[List[String]]] = {
    readTable(f, sep = "\t", head = false)
      .groupBy(x => x(3))
      .map((k, v) => {
        val r = if (v.length > n) {
          Random.shuffle(v).take(n)
        } else {
          v
        }
        (k, r)
      })
      .toMap
  }

  def computeSignalsCenterOnStates(bwfs: List[String], segmentfs: List[String], outf: String, ncore: Int = 5): Unit = {
    val params: List[String] =
      (s"-a 2500 -b 2500 -p ${ncore}" +
        " --referencePoint center --missingDataAsZero")
        .split(" ")
        .toList
    val commands: List[String] =
      List(computeMatBin, "reference-point", "-S")
        ::: bwfs
        ::: List("-R")
        ::: segmentfs
        ::: params ::: List("-out", outf)
    os.proc(commands.map(x => os.Shellable(List(x)))*)
      .call(check = true)
  }

  def main(args: Array[String]) = {
    // parallel here
    states.par.foreach(n => {
      // val modeld = s"${workd}/out/model_bypeakSPMq25_b200_s${n}"
      val modeld = s"${workd}/out/model_bypeak_b200_s${n}"
      scs.foreach(sc => {
        val outd = s"${modeld}/${sc}"
        if (!os.exists(os.Path(outd))) {
          os.makeDir(os.Path(outd))
        }
        val segmentf      = s"${modeld}/${sc}_${n}_segments.bed"

        val state2segment = splitSegments(segmentf, n = 50000)

        val state2segmentfs: Map[String, String] =
          1.to(n)
          .map(i => (s"E${i}", s"${outd}/E${i}.bed"))
          .toMap

        val outmatf = s"${outd}/${sc}.s${n}.computeMatrix.matrix"

        // save state2segment to segmentfs
        state2segment.foreach(
          (s, ss) => {
            writeStrings2File(
              content = ss.map(x => x.take(3).mkString("\t")),
              to = state2segmentfs(s),
              overwrite = true,
              head = ""
            )
          }
        )


        // start deepTools
        println(s"computeMatrix for ${sc} on state ${n}.")
        // under 10 cores for 20-30 minutes
        computeSignalsCenterOnStates(
            bwfs = sc2bwf(sc),
            segmentfs =
              1.to(n)
              .map(i => state2segmentfs(s"E${i}"))
              .toList,
            outf = outmatf
        )

        println(s"computeMatrix for ${sc} on state ${n}: done.")
        os.proc(List(plotHeatmapBin, "--matrixFile", outmatf,
                "--outFileName", s"${outd}/${sc}.s${n}.plotHeatmap.pdf",
                "--colorMap", "Greys", "Reds", "Oranges", "Greens",
                "Purples").map(x => os.Shellable(List(x)))).call(check = true)
        println(s"PlotHeatmap for ${sc} on state ${n}: done.")
      }) // end of scs
    })   // end of states
  }      // end of main
}        // end of object
