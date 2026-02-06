import os._
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*
import GRange.{GenomicRange, mouseGenomicRangeOrd}
import ChromHMM.ChromHMMState
import SZUtils.readTable
import SZUtils.writeStrings2File

object CalStateChangesArosccSubclass {
  type ChromHMMStates = Map[String, Vector[ChromHMMState]]

  def loadChromHMMState(densefnm: String): Vector[ChromHMMState] = {
    readTable(densefnm, sep = "\t", head = true)
      .map(x =>
        new ChromHMMState(name = x(3),
            new GenomicRange(chrom = x(0), startFrom = x(1).toInt,
                endTo = x(2).toInt)))
      .toVector
  }

  def groupChromHMMState(
    x: Vector[ChromHMMState]): Map[String, ChromHMMStates] = {
    x.groupBy(_.name)
      .toMap
      .map((k, v) =>
        (k,
            v.groupBy(_.g.chrom)
              .toMap
              .map((k, v) => (k, v.sortBy(_.g)))))
  }

  def statChromHMMStateChangesPerChrom(base: Vector[ChromHMMState],
    query: Vector[ChromHMMState]): Vector[Int] = {
    base.map(s => {
      // BUG: last element in the chrom would be ignored in this case
      // for the same subclass
      val indexFrom  = query.indexWhere(x => x.g.endTo > s.g.startFrom)
      val indexUntil = query.indexWhere(x => x.g.startFrom > s.g.endTo)
      if (indexFrom < 0 || indexUntil < 0) {
        0
      } else {
        val qs = query.slice(indexFrom, indexUntil)
        qs match {
          case Vector() => 0
          case _ =>
            qs.map(q =>
              s.g.overlap(q.g) match {
                case Some(v) => v.endTo - v.startFrom
                case None    => 0
              })
              .max
        }
      }
    })
  }

  def statChromHMMStateChanges(base: ChromHMMStates,
    query: ChromHMMStates): Vector[Int] = {
    val chrs = base.keySet.intersect(query.keySet).toList
    chrs.par
      .map(chr =>
        statChromHMMStateChangesPerChrom(
            base = base(chr),
            query = query(chr)
        ))
      .toVector
      .flatten
  }

  def sumChromHMMStateCoverage(x: ChromHMMStates): Double = {
    x.map((k, v) => (k, v.map(x => x.g.endTo - x.g.startFrom).sum))
      .values
      .map(x => x.toDouble)
      .sum
  }

  def getOvlpRatio(base: ChromHMMStates,
    query: ChromHMMStates): Double = {
    statChromHMMStateChanges(base, query).map(x => x.toDouble).sum /
      sumChromHMMStateCoverage(base)
  }

  def getnstate(x: ChromHMMStates): Double = {
    x.map((chr, ss) => (chr, ss.length))
      .map((chr, size) => size.toDouble)
      .sum
  }

  def getnovlpPerChrom(x: Vector[ChromHMMState],
    y: Vector[ChromHMMState]): Double = {
    x.map(s => {
      if (s.g.isInSortedRanges(y.map(t => t.g).toList)) {
        1.0.toDouble
      } else {
        0.0.toDouble
      }
    }).sum
  }

  def pgetnovlp(x: ChromHMMStates, y: ChromHMMStates): Double = {
    val chrs = x.keySet.intersect(y.keySet).toList.distinct
    chrs.par
      .map(chr => getnovlpPerChrom(x(chr), y(chr)))
      .toList
      .map(x => x.toDouble)
      .sum
  }

  def getOvlpJaccard(x: ChromHMMStates, y: ChromHMMStates): Double = {
    val nx = getnstate(x)
    // val ny = getnstate(y)
    val n = pgetnovlp(x, y)
    // when n > ny, this will be great than 1
    // n / (nx + ny - n)
    n / nx
  }

  def main(args: Array[String]) = {
    val projd  = "/tscc/projects/ps-renlab2/szu/projects/amb_pairedtag"
    val workd  = s"${projd}/06.ChromHMM"
    val densed = s"${workd}/out/updateDenseBed"
    val scs =
      os.list(os.Path(densed))
        .map(x => x.baseName)
        .map(x => x.replace("_18_dense", ""))
    val ustate =
      List("Chr-R", "ND", "Hc-H", "Chr-A", "Hc-P", "Chr-O", "Chr-P",
          "Chr-B")

    // val allSim = scs.par.map( sc1 => {
    //   val c1 = groupChromHMMState(
    //     loadChromHMMState(s"${densed}/${sc1}_18_dense.bed"))
    //   scs.par.map(sc2 => {
    //     val c2 = groupChromHMMState(
    //       loadChromHMMState(s"${densed}/${sc2}_18_dense.bed"))
    //     val r = ustate.map(s => {
    //       println(s"${sc1} ${sc2} ${s}...")
    //       val r = getOvlpRatio(base = c1(s),query = c2(s))
    //       s"${s}:${r.toString}"
    //     }).toList.mkString(";")
    //     s"${sc1},${sc2},${r}"
    //   }).toList
    // }).toList.flatten

    // writeListOfString2File(
    //   content = allSim,
    //   to = s"${workd}/out/chromHMMStateConservation.csv",
    //   head = "sc1,sc2,conservation"
    // )

    // // * check all the dense bed files have chr1-chr19s.
    // val sc2nchr = scs.par.map(sc => {
    //   val r = readTable(fnm = s"${densed}/${sc}_18_dense.bed", sep = "\t", head = true).map(
    //     x => x(0)
    //   ).distinct.length
    //   (sc, r)
    // }).toList
    // // all have 19 chrs

    // * use Jaccard Index
    val sc1 = args(0)
    val c1 = groupChromHMMState(
        loadChromHMMState(s"${densed}/${sc1}_18_dense.bed"))
    val allSim = scs.par
      .map(sc2 => {
        val c2 = groupChromHMMState(
            loadChromHMMState(s"${densed}/${sc2}_18_dense.bed"))
        val r = ustate
          .map(s => {
            println(s"${sc1} ${sc2} ${s}...")
            val r = getOvlpJaccard(x = c1(s), y = c2(s))
            s"${s}:${r.toString}"
          })
          .toList
          .mkString(";")
        s"${sc1},${sc2},${r}"
      })
      .toList

    writeStrings2File(
        content = allSim,
        to = s"${workd}/out/varOfState/${sc1}.chromHMMStateConservation.Jaccard.csv",
        head = "sc1,sc2,conservation"
    )

    // test
    // val sc1 = "001_CLA_EPd_CTX_Car3_Glut"

    // run in local
    // val allSim = scs.par.map(sc1 => {
    //   val c1 = groupChromHMMState(loadChromHMMState(s"${densed}/${sc1}_18_dense.bed"))
    //   scs.par.map(sc2 => {
    //     val c2 = groupChromHMMState(loadChromHMMState(s"${densed}/${sc2}_18_dense.bed"))
    //     val r = ustate.map(s => {
    //       println(s"${sc1} ${sc2} ${s}...")
    //       val r = getOvlpJaccard(x = c1(s), y = c2(s))
    //       s"${s}:${r.toString}"
    //     }).toList.mkString(";")
    //     s"${sc1},${sc2},${r}"
    //   }).toList
    // }).toList.flatten

    // writeListOfString2File(
    //   content = allSim,
    //   to = s"${workd}/out/varOfState/all.chromHMMStateConservation.Jaccard.csv",
    //   head = "sc1,sc2,conservation"
    // )

  } // end of main
}   // end of object
