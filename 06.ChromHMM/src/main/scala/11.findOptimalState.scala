import os._
import scala.math.*
import smile.read
import smile.data.DataFrame
import smile.math.matrix.Matrix
import smile.math.MathEx.median
import smile.clustering.kmeans
import smile.clustering.KMeans
import SZUtils.writeStrings2File

object FindOptimialChromHMMState {
  val chromHMMjar =
    "/projects/ps-renlab2/szu/softwares/ChromHMM/ChromHMM.jar"
  val projd = "/projects/ps-renlab2/szu/projects/amb_pairedtag"
  val workd = s"${projd}/06.ChromHMM"

  def getEmf(s: Int): String = {
    s"${workd}/out/model_bypeak_b200_s${s}/emissions_${s}.txt"
  }

  def sumSquare(x: Array[Double], y: Array[Double]): Double = {
    x.zip(y).map((x, y) => pow(x - y, 2)).sum
  }

  def copyEms2Emsd(fromd: String, states: List[Int], outd: String): Unit = {
    if (!os.exists(os.Path(outd))) {
      os.makeDir(os.Path(outd))
    }
    states.foreach { s =>
      {
        val emf = states.map(s => getEmf(s))
        os.proc("cp", "-t", outd, emf).call(check = true)
      }
    }
  }

  def runChromHMMCompareModels(refEmf: String, otherEmsd: String, outprefix: String): Unit = {
    os.proc(
        List("java", "-jar", chromHMMjar, "CompareModels", "-color",
            "255,0,0", refEmf, otherEmsd, outprefix).map(x =>
          os.Shellable(List(x)))*
    ).call(check = true)
  }

  def getMedianCors(compModelf: String): Vector[Double] = {
    val r  = read.csv(file = compModelf, "\t", header = true)
    val r2 = r.toMatrix().submatrix(1, 1, r.nrow() - 1, r.ncol() - 1)
    0.to(r2.ncol() - 1).map(i => median(r2.col(i))).toVector
  }

  def readMultEms(emsf: List[String]): DataFrame = {
    val rs = emsf.map(f => read.csv(file = f, "\t", header = true))
    val left = rs.tail.map(s => s.select(rs.head.names*))
    // work under smile 4.1.0
    // in smile 4.2.0, no union in 
    //rs.head.union(left*)
    rs.head.concat(left*)
  }

  def measureKMeans(km: KMeans, data: Array[Array[Double]]): Double = {
    val m = Matrix.of(data).colMeans()
    val cs = km.centroids
    val cs2n = km.y.groupBy(x => x).map((cs, v) => (cs, v.length)).toMap
    val totalss: Double = data.map(d => sumSquare(d, m)).sum
    val btwss: Double = cs.zipWithIndex.map(
      (d, cs) => cs2n(cs).toDouble * sumSquare(d, m)).sum
    btwss / (totalss + 1e-10)
  }

  def main(args: Array[String]) = {
    // run compareModels
    val states    = 5.to(24).toList
    val refEmf    = getEmf(24)
    val outd      = s"${workd}/out/compareModel_bypeak_b200"
    val outprefix = s"${outd}/compareModel"

    // Strategy 1: correlation
    copyEms2Emsd(s"${workd}/out", states, outd)
    runChromHMMCompareModels(refEmf, outd, outprefix)
    val compModelf = s"${outd}/compareModel.txt"
    val cors       = getMedianCors(compModelf)
    val s2cor = states.zip(cors)

    // Strategy 2: K-Means
    val em = readMultEms(states.map(s => getEmf(s)))
    val dataPoints = em.select(em.names().tail*).toMatrix()
    val kms = states.map(n => kmeans(dataPoints.toArray(),
      k = n, maxIter = 200))
    val bwSSOverTSS = kms.map(k => measureKMeans(k, data = dataPoints.toArray()))
    val ratio = bwSSOverTSS.zip(states).map(
       (x,s) => (s, x, x * 100 / bwSSOverTSS.last))
    val r = ratio.zip(s2cor).map( (x1, x2) => (x2._1, x2._2, x1._2, x1._3) )
    writeStrings2File(content = r.map(x => x.toList.mkString(",")),
      to = s"${outd}/evalChromHMM.cor.kmeans.csv",
      head = "state,medianCor,bwtSS,ssratio")

  }
}
