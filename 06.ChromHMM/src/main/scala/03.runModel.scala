import os._
import SZUtils.writeStrings2File

object RunChromHMModel {
  val chromHMMJar = "/projects/ps-renlab2/szu/softwares/ChromHMM/ChromHMM.jar"
  val projd = "/projects/ps-renlab2/szu/projects/amb_pairedtag"
  val chromHMMd = s"${projd}/data/chromHMM"

  val binSize: Int = 200
  val nstate: Int = 10
  val ncpu: Int = 30
  def main(args: Array[String]) = {
    val binarizeBedir = s"${projd}/06.ChromHMM/out/binarizeBed_${binSize}"
    val flagfnm =
      s"${projd}/06.ChromHMM/out/flag/runModel_b${binSize}_s${nstate}.done"
    val logfnm =
      s"${projd}/06.ChromHMM/out/log/runModel_b${binSize}_s${nstate}.log"
    val outd = s"${projd}/06.ChromHMM/out/ChromHMModel_b${binSize}_s${nstate}"
    if (!os.exists(os.Path(flagfnm))) {
      os.proc(
        "java",
        "-mx800000M",
        "-jar",
        chromHMMJar,
        "LearnModel",
        "-p",
        ncpu,
        binarizeBedir,
        outd,
        nstate,
        "mm10"
      ).call(check = true, stderr = os.Path(logfnm), stdout = os.Path(logfnm))
      os.proc("touch", flagfnm).call(check = true)
      println("Finish runChromHMModel. Good luck!")
    } else {
      println("ChromHMM has finished.")
    }

  }

}
