import os._
import SZUtils.{str2path, path2str, ifelse}
import SZUtils.writeStrings2File
import bioscala.LightCoord.Bed.{
  LightBedElement, mkStringLightBedElement
}
import scala.collection.parallel.*
import scala.collection.parallel.CollectionConverters.*

object PlotHeatmapTE {
  val projd = "/Users/szu/git-recipes/amb_pairedtag"
  val workd = os.Path(projd) / "17.repressiveMarks"
  val rscd  = workd / "src" / "main" / "resources"
  val TEd   = rscd / "repeats"
  val sTEClasses: Vector[String] = Vector("LINE", "LTR")
  val chrs                       = 1.to(19).map(i => s"chr$i")
  val TEBedPath                  = rscd / "repeatBed"
  val npart: Int                 = 4
  val minTE: Int                 = 20

  // * folders for data and output
  val logd = workd / "log"
  val flagd = workd / "flagd"
  val outd = workd / "out" / "TEHeatmap"

  val bwd: String =
    os.Path(projd) / "data" / "bw" / "bw_PairedTag_DNA"

  val ascs: Vector[String] = Vector(
      "016_CA1_ProS_Glut",
      "017_CA3_Glut",
      // "001_CLA_EPd_CTX_Car3_Glut",
      "037_DG_Glut",
      //"010_IT_AON_TT_DP_Glut",
      //"002_IT_EP_CLA_Glut",
      "007_L2_3_IT_CTX_Glut",
      "008_L2_3_IT_ENT_Glut",
      //"009_L2_3_IT_PIR_ENTl_Glut",
      // "019_L2_3_IT_PPP_Glut",
      "011_L2_IT_ENT_po_Glut",
      "006_L4_5_IT_CTX_Glut",
      //"003_L5_6_IT_TPE_ENT_Glut",
      "022_L5_ET_CTX_Glut",
      "005_L5_IT_CTX_Glut",
      "030_L6_CT_CTX_Glut",
      "004_L6_IT_CTX_Glut",
      "029_L6b_CTX_Glut"
      //"014_LA_BLA_BMA_PA_Glut",
      // "151_TH_Prkcd_Grin2c_Glut"
  )

  val bgscs: Vector[String] = Vector(
      "318_Astro_NT_NN",
      "319_Astro_TE_NN",
      "334_Microglia_NN",
      "326_OPC_NN",
      "327_Oligo_NN",
      "053_Sst_Gaba",
      "052_Pvalb_Gaba",
      "049_Lamp5_Gaba"
  )
  // ATAC, H3K27ac, H3K27me3, H3K4me1, H3K9me3
  // val colors = List(
  //   "Greys", "Reds", "Oranges", "Greens", "Purples"
  // )

  // H3K27ac, H3K27me3, H3K4me1, H3K9me3
  val colors = List(
      "Reds",
      "Oranges",
      "Greens",
      "Purples"
  )

  val h2sm: Vector[(String, Int)] = Vector(
      ("H3K27ac", 300),
      ("H3K27me3", 300),
      ("H3K4me1", 300),
      ("H3K9me3", 1000)
  )

  val histones = h2sm.map((x, y) => x).toList

  // * Applications
  val deepTools: String =
    "/Users/szu/miniforge/envs/deeptools/bin/deeptools"

  val computeMatrix: String =
    "/Users/szu/miniforge/envs/deeptools/bin/computeMatrix"

  val plotHeatmap: String =
    "/Users/szu/miniforge/envs/deeptools/bin/plotHeatmap"

  // * functions
  def getTEClass(x: String): String = {
    x.split("\\|")(1)
  }

  def getTEFamily(x: String): String = {
    x.split("\\|")(2)
  }

  def getRawTESubfamily: Vector[String] = {
    os.list(TEd)
      .filter(x => os.isFile(x))
      .map(x => x.baseName)
      .filter(x => x.contains(".ann"))
      .map(x => x.replace(".ann", ""))
      .filter(x => x.contains("|") && (!x.contains("?")))
      .filter(x => x.split("\\|").length > 2)
      .filter(x => sTEClasses.contains(getTEClass(x)))
      .toVector
  }

  def getTEClassGroup(g: Iterable[String]) = {
    g.map(x => (getTEClass(x), x))
      .groupMap(_._1)(_._2)
  }

  def getTEFamilyGroup(g: Iterable[String]) = {
    g.map(x => (getTEFamily(x), x))
      .groupMap(_._1)(_._2)
  }

  def getRegionBodyLength(f: String): Int = ???

  def getRawTEfnm(s: String): String = {
    s"$TEd/$s.ann.txt"
  }

  def genTEBedfnmFromRawTEfnm(f: String) = {
    f.replace("txt", "bed")
  }

  def getRawTEfnms(g: Iterable[String]) = {
    g.map(x => getRawTEfnm(x)).toList
  }

  def transformRawTEtxt2Bed(f: String, overwrite: Boolean = false): Int = {
    val elements = os.read.lines
      .stream(f)
      .filter(x => x.nonEmpty)
      .map(x => x.split("\t"))
      .filter(x => chrs.contains(x(1)))
      .map(x =>
        (g = (
                chrom = x(1),
                coord = (startFrom = x(2).toInt,
                    endTo = x(3).toInt),
                strand = ifelse(x(4).toInt == 0, "+", "-")
            ), name = x(0), score = 1.0))
      .toVector
    val content = elements.map(x => mkStringLightBedElement(x))
    val outf = genTEBedfnmFromRawTEfnm(f)
    if ( !os.exists(outf) || overwrite) {
      writeStrings2File(content, to = outf,
        overwrite = overwrite, head = "")
    }
    val s =
      elements.map(x =>
        (x.g.coord.endTo - x.g.coord.startFrom).toDouble).sum
    (s /  elements.length.toDouble).toInt
  }

  def getTEbedfnms(g: Iterable[String]) = {
    g.map(x => s"$TEd/$x.ann.bed").toList
  }

  def getbwfnms(sc: String) = {
    h2sm
      .map(x => {
        s"$bwd/$sc.${x._1}.e100.bs100.sm${x._2}.bw"
      })
      .toList
  }

  def runComputeMatrix(bedfnms: List[String], regionSize: Int, bwfnms: List[String], outfnm: String, logfnm: String) = {
    val t = List(
        computeMatrix,
        "scale-regions",
        "-S"
    ) ::: bwfnms ::: List("-R") ::: bedfnms :::
      List("-a", "1000", "-b", "1000", "-p", "1",
          "--regionBodyLength", (regionSize / 10 * 10).toString,
          "--missingDataAsZero",
          "-out", outfnm)
    val cmd = t.map(x => os.Shellable(List(x)))
    os.proc(cmd*)
      .call(check = true, stdout = os.Path(logfnm),
          stderr = os.Path(logfnm))
  }

  def runPlotHeatmap(title: String, k:Int = 4,matfnm: String, outfnm: String, logfnm: String) = {
    val t = List(plotHeatmap, "--matrixFile", matfnm,
      "--outFileName", outfnm, "--colorMap") ::: colors :::
    List("--samplesLabel") ::: histones :::
    List("--startLabel", "u", "--endLabel", "d", "-T", title,
      "--legendLocation", "lower-center", "--kmeans", k.toString)
    val cmd = t.map(x => os.Shellable(List(x)))
    os.proc(cmd*)
      .call(check = true, stdout = os.Path(logfnm),
        stderr = os.Path(logfnm))
  }

  // * data classes
  case class TESubFamily(name: String, rawfnm: String,
    bedfnm: String, size: Int)

  // * main
  def main(args: Array[String]) = {

    val TESubfmls =
      getRawTESubfamily
        .map(x => (x, getRawTEfnm(x)))
        .filter(x => os.exists(x._2))
        .map(x =>
          TESubFamily(
              name = x._1,
              rawfnm = x._2,
              bedfnm = genTEBedfnmFromRawTEfnm(x._2),
              size = transformRawTEtxt2Bed(x._2, overwrite = false)
          ))
        .filter(x => x.size >= minTE)

    // test
    // NOTE: during test, HAL1M8|LINE|L1.ann.txt was replaced by
    // bed file, so I move it to repeatBed, and no
    // HAL1M8|LINE|L1.ann.txt in the original directory.

    // val k  = "L1"
    // val sc = ascs(0)
    // val te = TESubfmls(0)

    // // * run plotHeatmap
    val scs = ascs ++ bgscs
    TESubfmls.foreach(te => {
      scs.par.foreach { sc =>
        {
          val tenm           = te.name
          val bedfnm         = te.bedfnm
          val bwfnms         = getbwfnms(sc)
          val prefix         = s"$sc.$tenm"
          val logfnm: String =
            logd / s"$prefix.TE.plotHeatmap.log"
          val flagfnm: String =
            flagd / s"$prefix.TE.plotHeatmap.done"
          if (!os.exists(flagfnm)) {
            val matfnm: String = outd / s"$prefix.TE.mat"
            val figfnm: String = outd / s"$prefix.TE.heatmap.pdf"
            runComputeMatrix( List(bedfnm), te.size, bwfnms, matfnm,
              logfnm)
            if (os.exists(matfnm)) {
              val r = runPlotHeatmap(
                s"$sc: ${te.name}:${te.size}bp", 2, matfnm, figfnm, logfnm)
              println(s"Done: $sc: ${te.name}:${te.size}bp.")
            } else {
              println(s"Error: cannot find $matfnm.")
            }
          }} // end if flagfnm check
      } // end of scs foreach
    })  // end of TESubfmls foreach
  }     // end of main
}       // end of object
