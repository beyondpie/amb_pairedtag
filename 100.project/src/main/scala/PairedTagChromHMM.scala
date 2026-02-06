package PairedTagChromHMM

import GRange.GenomicRange
import ChromHMM.PeakChromHMMStateAnnot
import Bed.BedElement4
import MetaData.TSCCMeta.projd
import SZUtils.readTable

val stateAnnotfnm =
  s"$projd/meta/chromHMM18statesAnnotation.csv"

val chromHMMStates = List("Chr-A", "Chr-B", "Chr-R", "Chr-P", "Chr-O",
    "Hc-P", "Hc-H", "ND")

val cellMarkFileTable =
  s"$projd/data/chromHMM/subclass_peak/cellmarkfiletable.tsv"

val CREAnnotd =
  s"$projd/06.ChromHMM/out/CREAnnot250429"

val genomicRegionAnnotd =
  s"$projd/06.ChromHMM/out/model_bypeak_b200_s18"

case class ChromHMMStateSum(size: Int, name: String) {
  override def toString: String = {
    s"${size.toString}:${name}"
  }
}

case class PeakChromHMMStateSum(p: GenomicRange,
  a: Vector[ChromHMMStateSum])


case class ChromHMMStateAnnot(
  stateid: Int, groupid: Int, groupName: String, fullGroupName: String,
  color: String
)

/**
 * Load 18 ChromHMM states to their annotations.
 *
 * @param fnm
 * @return
 */
def loadChromHMMStateAnnots(
  fnm: String = stateAnnotfnm): Vector[ChromHMMStateAnnot] = {
  readTable(fnm, sep = ",", head = true)
    .map(x =>
      new ChromHMMStateAnnot(
          stateid = x(0).toInt,
          groupid = x(1).toInt,
          groupName = x(2),
          fullGroupName = x(3),
          color = x(4)
      ))
    .toVector
}

/**
  * Read CRE's ChromHMM annotation on the summit for a given subclass.
  *
  * @param sc
  * @return
  */
def readCREChromHMMAnnot(sc: String): Vector[BedElement4] = {
  val fnm = s"$CREAnnotd/$sc.CRE.annotBySummit.bed"
  BedElement4.readBed4(fnm).toVector
}
