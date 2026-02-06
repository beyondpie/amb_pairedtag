import os._
import SZUtils.writeStrings2File

object PrepareBinarzieBedFiles {
  // * meta
  val chromHMMJar =
    "/projects/ps-renlab2/szu/softwares/ChromHMM/ChromHMM.jar"
  val chromSizefnm =
    "/projects/ps-renlab2/szu/softwares/ChromHMM/CHROMSIZES/mm10.txt"

  val projd     = "/projects/ps-renlab2/szu/projects/amb_pairedtag"
  val chromHMMd = s"${projd}/data/chromHMM"
  val allbedir  = s"${chromHMMd}/all_bed"

  val ptbedir = s"${chromHMMd}/pt_subclass"
  val atacd   = s"${chromHMMd}/snATAC_subclass"

  val datafnm = s"${allbedir}/cellmarkfiletable.txt"

  def get_histone_scs(hm: String): List[String] = {
    os.list(os.Path(ptbedir))
      .map(x => x.baseName)
      .filter(_.contains(hm))
      .map(x => x.split("\\.")(0))
      .toList
  }
  def soft_link(
    fnm: String, fromd: String, tod: String, outfnm: String = ""
  ): Unit = {
    val tofnm = if (outfnm.length < 1) {
      s"${tod}/${fnm}"
    } else {
      s"${tod}/${outfnm}"
    }

    if (!os.exists(os.Path(tofnm))) {
      os.proc("ln", "-s", s"${fromd}/${fnm}", tofnm).call(check = true)
    } else {
      println(s"${tofnm} exists, and skip it.")
    }
  }

  def main(args: Array[String]) = {
    // val H3K27ac_scs = get_histone_scs("H3K27ac")
    // val H3K27me3_scs = get_histone_scs("H3K27me3")
    // val H3K4me1_scs = get_histone_scs("H3K4me1")
    // val H3K9me3_scs = get_histone_scs("H3K9me3")

    // val pt_scs =
    //   (H3K27ac_scs ::: H3K27me3_scs :::
    //     H3K4me1_scs ::: H3K9me3_scs).sorted.distinct
    // val pt_files: List[String] =
    //   os.list(os.Path(ptbedir))
    //     .map(x => x.baseName)
    //     .toList
    // val atac_files: List[String] =
    //   os.list(os.Path(atacd))
    //     .map(x => x.baseName)
    //     .toList

    // // * soft link bed files into one direcotry
    // pt_files.foreach { f =>
    //   soft_link(s"${f}.bed", ptbedir, allbedir)
    // }
    // atac_files.foreach { f =>
    //   {
    //     val sc = f.split("\\.")(0)
    //     soft_link(s"${f}.bed", atacd, allbedir, s"${sc}.ATAC.chromHMM.bed")
    //   }
    // }

    // // * prepare tab-sep files to describe the data
    // val bedInfo: List[List[String]] =
    //   os.list(os.Path(allbedir))
    //     .map(f => {
    //       val x = f.baseName
    //       val s = x.split("\\.")
    //       List(s(0), s(1), f"${x}.bed")
    //     })
    //     .toList

    // writeListOfString2File(
    //   content = bedInfo.map(x => x.mkString("\t")),
    //   to = datafnm,
    //   overwrite = true
    // )

    // * run ChromHMM BinarizeBed
    val binSize: Int = 200
    val binarizeBedir =
      s"${projd}/06.ChromHMM/out/binarizeBed_${binSize}"
    val flagfnm =
      s"${projd}/06.ChromHMM/out/flag/binarzieBed_${binSize}.done"
    val logfnm =
      s"${projd}/06.ChromHMM/out/log/binarizeBed_${binSize}.log"
    if (!os.exists(os.Path(flagfnm))) {
      os.proc(
          "java",
          "-mx500000M",
          "-jar",
          chromHMMJar,
          "BinarizeBed",
          "-b",
          binSize,
          chromSizefnm,
          allbedir,
          datafnm,
          binarizeBedir
      ).call(check = true, stderr = os.Path(logfnm),
          stdout = os.Path(logfnm))
      os.proc("touch", flagfnm).call(check = true)
      println("BinarizeBed is done.")
    } else {
      println("BinarizeBed has been performed. Skip.")
    }

  }
}
