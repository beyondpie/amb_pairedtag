import os._
import SZUtils.{readTable, ifelse}
import SZUtils.writeStrings2File

object eRNANaiveDEAnlaysis {
  val projd = "/projects/ps-renlab2/szu/projects/amb_pairedtag"
  val workd = s"${projd}/13.cicero"
  val outd = s"${workd}/out/eRNA"
  val eps = 1e-10

  def getFoldChange(s: Vector[Double]): Vector[Double] = {
    val sum = s.sum
    val n = s.length
    s.map(x => {
      val bg = (sum - x) / (n - 1)
      ifelse(
        test = (x <= eps || (bg <= eps)),
        yesValue = 0.0,
        noValue = (x / bg).min(10.0)
      )
    })
  }

  def main(args: Array[String]) = {
    val fd  =
      os.read.lines.stream(os.Path(s"$outd/mat/mat.csv"))
        .map(x => x.strip().split(","))
        .map(x => x.map(y => y.toDouble).toVector)
        .zipWithIndex.map((x, i) => {
          println(s"run naive fold change for ${i}.")
          getFoldChange(x)
        })
        .toVector
    writeStrings2File(
      content = fd.map(x => x.map(y => y.toString).mkString(",")),
      to = s"${outd}/mat/fd.csv",
      overwrite = true,
      head = ""
    )
  }
}
