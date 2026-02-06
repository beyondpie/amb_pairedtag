import scala.collection.immutable.HashMap
import os._
import java.io.File
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*

import Peak.{Peak, Peaks}
import GRange.{GenomicRange, mouseGenomicRangeOrd}
import SZUtils.{hereByGit, getFileSize, TaskElement}
import SZUtils.{readTable, writeStrings2File}
import MetaData.TSCCMeta
import project.Peak._


object RunGetRobustPeak {
  def main(args: Array[String]) = {
    val projd = TSCCMeta.projd
    val workd = s"${projd}/04.peakcalling"
    val from = s"${workd}/out/macs3_nolambda"
    val toReproPeak = s"${workd}/out/reproPeak_nolambda"
    val reproPeakFlagd = s"${workd}/flag/reproPeak_nolambda"

    val groups =
      os.list(os.Path(from)).map(x => x.toString.split("/").last)
    groups
      .filter(g => !ignoredGroup.values.toList.contains(g))
      .par
      .foreach(g =>
        TSCCMeta.modality.foreach(m => {
          List("Male", "Female").foreach { sex =>
            {
              val r = new GroupOfPeaks(
                from = from,
                name = g,
                sex = sex,
                modality = m,
                neglog10qval = 2,
                outd = toReproPeak,
                flagd = reproPeakFlagd,
                skipTask = true
              )
              r.run()
            }
          }
        })
      ) // end of groups foreach
    // test
    // val name = "001_CLA_EPd_CTX_Car3_Glut"
    // val sex = "Male"
    // val modality = "H3K27ac"
    // val r = new GroupOfPeaks(
    //   from = from,
    //   name = name,
    //   sex = sex,
    //   modality = modality,
    //   neglog10qval = 2,
    //   outd = toReproPeak,
    //   flagd = reproPeakFlagd,
    //   skipTask = true
    // )
    // val tp = r.getTreatPeaks
    // val rp = r.getReproPeak
    // val fnms = r.fnms.tail
    // val name =
    //   fnms(0).toString.split("/").last.split("_").toList.init.mkString("_")
    // group2nread.isDefinedAt(name)
    // group2nread(name)

  } // end of main
}
