import os._
import GRange.{GenomicRange, mouseGenomicRangeOrd}
import SZUtils.{readTable, writeMap2File}
import ChromHMM.ChromHMMState
import ChromHMM.PeakChromHMMStateAnnot as PeakStateAnnot

object MapChromHMM2ATACPeak {

  /**
   * Get the 8 states information on each subclass's genomic bins.
   *
   * NOTE: ignore ND state.
   *
   * @param scfnm
   * @return
   */
  def getChromHMMStateFromDenseBed(
    scfnm: String): Map[String, Vector[ChromHMMState]] = {
    readTable(fnm = scfnm, sep = "\t", head = true)
      .filter(x => x(3) != "ND")
      .map(x =>
        new ChromHMMState(
            name = x(3),
            g = new GenomicRange(chrom = x(0), startFrom = x(1).toInt,
                endTo = x(2).toInt)
        ))
      .groupBy(_.g.chrom)
      .map((k, v) => (k, v.sortBy(x => x.g).toVector))
  }

  def loadATACPeak(scfnm: String,
    ext: Int = 250): Map[String, List[GenomicRange]] = {
    readTable(fnm = scfnm, sep = "\t", head = false)
      .map(x => {
        val l = x(1).toInt
        val r = x(2).toInt
        val s = l + math.floor((r - l) / 2).toInt
        // FIXME: startFrom is one base left compared with the raw one.
        // we manually fix the summit annotatin later
        // but in the future, for extenstion, we should correct it
        // from code, like put it in GenomicRange.
        new GenomicRange(chrom = x(0), startFrom = s - ext,
            endTo = s + ext)
      })
      .groupBy(_.chrom)
      .map((k, v) => (k, v.sorted.toList))
  }

  def mapChromHMMStateOnATACPeak(p: Map[String, List[GenomicRange]],
    s: Map[String,
        Vector[ChromHMMState]]): Map[String, List[PeakStateAnnot]] = {
    p.map((chr, gs) => {
      if (!s.contains(chr)) {
        val r = gs.map(ele => {
          new PeakStateAnnot(p = ele,
              a = Vector(new ChromHMMState(name = "ND", g = ele)))
        })
        (chr, r)
      } else {
        val ss = s(chr)
        val r = gs.foldLeft((List.empty[PeakStateAnnot], 0))(
            (lastState, ele) => {
              val (acc, preIndex) = lastState
              val s3 = ss.zipWithIndex
                .slice(preIndex, ss.length)
                .takeWhile((x, id) => x.g.startFrom < ele.endTo)
              val r = if (s3.isEmpty) {
                (new PeakStateAnnot(p = ele,
                        a = Vector(new ChromHMMState(name = "ND",
                                g = ele))), ss.length)
              } else {
                val s4 = s3.filter(_._1.g.endTo >= ele.startFrom)
                if (s4.isEmpty) {
                  (new PeakStateAnnot(p = ele,
                          a = Vector(new ChromHMMState(name = "ND",
                                  g = ele))), s3.last._2 + 1)
                } else {
                  (new PeakStateAnnot(p = ele,
                          a = s4
                            .map(_._1)
                            .map(x =>
                              new ChromHMMState(
                                  name = x.name,
                                  g = new GenomicRange(
                                      chrom = x.g.chrom,
                                      startFrom = x.g.startFrom.max(
                                          ele.startFrom),
                                      endTo = x.g.endTo.min(ele.endTo))
                              ))), s4.head._2)
                }
              } // end of creation of PeakStateAnnot
              // return the foldLeft op result (update the states)
              (acc :+ r._1, r._2)
            }) // end of foldLeft
        (chr, r._1)
      } // end of ifelse
    })
  }

  def main(args: Array[String]) = {
    val projd = "/projects/ps-renlab2/szu/projects/amb_pairedtag"
    val workd = s"${projd}/13.cicero"
    val chromHMMDenseBedir = s"${projd}/data/chromHMM/denseBed"
    val ATACPeakd          = s"${projd}/data/chromHMM/subclass_peak"
    val outd               = s"${workd}/out/ATACPeakChromHMMAnnot"
    val chrs               = 1.to(19).map(i => s"chr${i}")

    val scs: List[String] =
      os.list(os.Path(ATACPeakd))
        .map(x => x.toString.split("/").last)
        .filter(x => x.contains("ATAC.peak.bed"))
        .map(x => x.split("\\.").head)
        .toList

    val exts: List[Int] = List(250, 500, 750, 1000)
    scs.slice(1, scs.length).foreach { sc =>
      {
        val states = getChromHMMStateFromDenseBed(
            scfnm = s"${chromHMMDenseBedir}/${sc}_18_dense.bed")
        exts.foreach { e =>
          {
            // test
            // val sc = scs(0)
            println(
                s"map ChromHMM 8 states to ATAC peaks: ${sc} for ext ${e}.")
            val peaks =
              loadATACPeak(scfnm = s"${ATACPeakd}/${sc}.ATAC.peak.bed",
                  ext = e)

            val annPeaks = mapChromHMMStateOnATACPeak(peaks, states)

            val content =
              annPeaks.map((k, v) => (k, v.map(x => x.toString())))

            writeMap2File(content = content,
                out = os.Path(
                    s"${outd}/${sc}.chromHMM8state.annot.ATAC.peakSummit.ext${e}.csv"),
                overwrite = true, ordedKeys = chrs.toList,
                head = "CRE,ChromHMMAnnot")
          }
        }
      }
    } // end of foreach
  }   // end of main
}     // end of object
