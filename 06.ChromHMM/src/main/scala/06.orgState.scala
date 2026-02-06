import os._
import SZUtils.readTable
import SZUtils.writeStrings2File
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*

object OrgChromHMMState {
  val projd  = "/projects/ps-renlab2/szu/projects/amb_pairedtag"
  val workd  = s"${projd}/06.ChromHMM"
  val modeld = s"${workd}/out/model_bypeak_b200_s18"
  val stated = s"${workd}/out/m-bpeak-s18_pd-obs"
  val chrs   = 1.to(19).map(x => s"chr${x}").toList
  val sep    = ","
  val homerOrdElements = List("3UTR", "5UTR", "CpG", "Exon", "GAP",
      "Intron", "Intergenic", "TES", "TSS", "LINE", "LTR", "SINE",
      "Satellite", "SimpleRepeat", "rRNA", "tRNA")

  val genomeBinHomerAnnotd = s"${projd}/data/genome/mm10GenomeBinAnnot"
  val homerAnnotOutd       = s"${workd}/out/m-bpeak-s18_homer_bannot"

  case class ChromHMMState(stateId: Int, pd: String, groupId: Int,
    name: String, note: String)

  case class ChromBinStateAnnot(chr: String, startFrom: String,
    endTo: String, binId: Int, relaTSS: String, state: Int, bATAC: Int,
    bH3K27ac: Int, bH3K27me3: Int, bH3K4me1: Int, bH3K9me3: Int)

  case class GroupSum(groupId: Int, nATAC: Int, nH3K27ac: Int,
    nH3K27me3: Int, nH3K4me1: Int, nH3K9me3: Int, n: Int) {
    def toString(sep: String = ","): String = {
      List(groupId, nATAC, nH3K27ac, nH3K27me3, nH3K4me1, nH3K9me3, n)
        .mkString(sep)
    }
  }

  def headOfGroupSum(sep: String): String = {
    List("groupId", "nATAC", "nH3K27ac", "nH3K27me3", "nH3K4me1",
        "nH3K9me3", "totalCount").mkString(sep)
  }

  def loadChromHMMState(fnm: String): List[ChromHMMState] = {
    readTable(fnm, sep = ",", head = true).map(x =>
      new ChromHMMState(
          stateId = x(0).toInt,
          pd = x(1),
          groupId = x(2).toInt,
          name = x(3),
          note = x(4)
      ))
  }

  /**
    * Use stateId + pd as groupId instead of the groupId in the file.
    *
    * Here "P"(promoter) will be mapped to 100, and "D"(distal)
    *  will be mapped to 0.
    * @param fnm
    * @return
    */
  def trivalLoadChromHMMState(nstate: Int = 18): List[ChromHMMState] = {
    1.to(nstate).map(i => {
      List(
        new ChromHMMState(
          stateId = i,
          pd = "P",
          groupId = 100 + i,
          name = s"${i}_P",
          note = "CRE_Promoter"
        ),
        new ChromHMMState(
          stateId = i,
          pd = "D",
          groupId = 0 + i,
          name = s"${i}_D",
          note = "CRE_Promoter"
        )
      )
    }).flatten.toList
  }

  /**
   * Read ChromBinStateAnnot file for one chromosome.
   *
   * @param fnm
   * @return
   */
  def loadChromBinStateAnnot(fnm: String): List[ChromBinStateAnnot] = {
    readTable(fnm, sep = ",", head = true).map(x =>
      new ChromBinStateAnnot(chr = x(0), startFrom = x(1), endTo = x(2),
          binId = x(3).toInt, relaTSS = x(4), state = x(5).toInt,
          bATAC = x(6).toInt, bH3K27ac = x(7).toInt,
          bH3K27me3 = x(8).toInt, bH3K4me1 = x(9).toInt,
          bH3K9me3 = x(10).toInt))
  }

  def getGroupDataSumOnChrom(i: List[ChromBinStateAnnot], state2Group: Map[(Int, String), Int]): List[GroupSum] = {
    i.map(x => (state2Group((x.state, x.relaTSS)), x))
      .groupBy((k, x) => k)
      .map((k, la) => {
        new GroupSum(
            groupId = k,
            nATAC = la.map((_, a) => a.bATAC).sum,
            nH3K27ac = la.map((_, a) => a.bH3K27ac).sum,
            nH3K27me3 = la.map((_, a) => a.bH3K27me3).sum,
            nH3K4me1 = la.map((_, a) => a.bH3K4me1).sum,
            nH3K9me3 = la.map((_, a) => a.bH3K9me3).sum,
            n = la.length
        )
      })
      .toList
  }

  def mapGroupSums2GroupSums(a: List[GroupSum]): List[GroupSum] = {
    a.groupBy(x => x.groupId)
      .map((k, v) =>
        new GroupSum(
            groupId = k,
            nATAC = v.map(x => x.nATAC).sum,
            nH3K27ac = v.map(x => x.nH3K27ac).sum,
            nH3K27me3 = v.map(x => x.nH3K27me3).sum,
            nH3K4me1 = v.map(x => x.nH3K4me1).sum,
            nH3K9me3 = v.map(x => x.nH3K9me3).sum,
            n = v.map(x => x.n).sum
        ))
      .toList
  }

  /**
   * Load mm10 Genome Bin Annotation for one chromsome
   *
   * @param fnm
   *   for one chromosome
   * @return
   */
  def loadmm10GenomeBinAnnot(fnm: String): (List[String], Vector[Vector[Int]]) = {
    val h = os.read.lines
      .stream(os.Path(fnm))
      .slice(0, 2)
      .head
      .split(",")
      .toList
    val v = readTable(fnm, sep = ",", head = true)
      .map(x => x.map(y => y.toInt).toVector)
      .toVector
    (h, v)
  }

  def main(args: Array[String]) = {

    // val state2group =
    //   loadChromHMMState(s"${modeld}/org.chromHMM18states.csv")
    //     .map(x => ((x.stateId, x.pd), x.groupId))
    //     .toMap

    val state2group =
      trivalLoadChromHMMState(nstate = 18)
        .map(x => ((x.stateId, x.pd), x.groupId))
        .toMap

    // all the subclasses
    val scs: List[String] =
      os.list(os.Path(stated))
        .map(x => x.baseName)
        .filter(x => x.contains("chr1_s18-200bin-pd-obs"))
        .distinct
        .map(x => x.replace("_chr1_s18-200bin-pd-obs", ""))
        .toList

    // groupSums (merge all the chrs) for each subclass
    // Under full paral, needs 30 mins to finish
    val scsGroupSums: Map[String, List[GroupSum]] =
      scs.par
        .map(sc => {
          println(s"Get groupSum for ${sc} ... ")
          val aa =
            chrs
              .map(chr => {
                println(s"Get groupSum for ${sc} on chr ${chr}.")
                val f = s"${stated}/${sc}_${chr}_s18-200bin-pd-obs.csv"
                val i = loadChromBinStateAnnot(f)
                val r = getGroupDataSumOnChrom(i, state2group)
                println(s"Get groupSum for ${sc} on chr ${chr}: done.")
                r
              })
              .flatten
          println(s"Get groupSum for ${sc}: done.")
          val r    = (sc, mapGroupSums2GroupSums(aa))
          val outf = s"${stated}/${sc}.stateOrgSum.csv"
          writeStrings2File(
              content = r._2.map(g => g.toString(sep)),
              to = outf,
              overwrite = true,
              head = headOfGroupSum(sep)
          )
          r
        })
        .toList
        .toMap

    // then merge all the subclasses' results.
    val groupSums: List[GroupSum] =
      mapGroupSums2GroupSums(
          scsGroupSums.map((_, v) => v).toList.flatten)
    // save both groupSum result
    writeStrings2File(
        content = groupSums.map(g => g.toString(sep)),
        to = s"${stated}/all.stateOrgSum.csv",
        overwrite = true,
        head = headOfGroupSum(sep)
    )

    // get genome bin annotations
    val chr2BinHomerAnnot: Map[String, Vector[Vector[Int]]] =
      chrs
      .map(chr => {
        val annotfnm =
          s"${genomeBinHomerAnnotd}/${chr}.binary.annot.csv"
        val a   = loadmm10GenomeBinAnnot(annotfnm)
        val aa  = a._1.zip(a._2.transpose).toMap
        val aaa = homerOrdElements.map(x => aa(x)).toVector
        (chr, aaa.transpose)
      })
      .toMap
    // val annotfnm = s"${genomeBinHomerAnnotd}/${chr}.binary.annot.csv"

    // get the homer annotations for each states
    // val scsGroupHomerAnnotSums: Map[String, List[Int]] =
    val scsGroupHomerAnnotSums =
      scs.par
        .map(sc => {
          // groupId, counts of homerOrdElements across chrs, size
          val a: List[(Int, Vector[Int], Int)] =
            chrs
              .map(chr => {
              println(s"Get HomerAnnotSum for ${sc} on ${chr}.")
              val f = s"${stated}/${sc}_${chr}_s18-200bin-pd-obs.csv"
              val t = readTable(f, sep = ",",head = true).map(
                x => (x(5).toInt, x(4)))
              t.zipWithIndex
                .groupBy(x => x._1)
                .map((k, v) => {
                  val index: Vector[Int] = v.map(i => i._2).toVector
                  val r: Vector[Vector[Int]] =
                    index.map(chr2BinHomerAnnot(chr))
                  // size of homerOrdElements.
                  val rr: Vector[Int] = r.transpose.map(x => x.sum)
                  (state2group(k), rr, index.length)
                })
                .toList
            })
            .flatten
            .groupBy(x => x._1)
            .map(
              // v: (groupId, #chrs x #homerOrdElement)
              (k, v) => {
                (k,
                  v.map(x => x._2)
                    .toVector
                    .transpose
                    .map(x => x.sum), v.map(x => x._3).sum)
              })
            .toList
          println(s"Done for ${sc}...")
          // save a to fnm
          val outf =
            s"${homerAnnotOutd}/${sc}.groupId2HomeGenomeAnnot.csv"
          writeStrings2File(
              content = a.map(x =>
                (x._1 +: x._2.toList :+ x._3).mkString(",")),
              to = outf,
              overwrite = true,
              head = ("groupId" +: homerOrdElements :+ "nTotal")
                .mkString(",")
          )
          a
        })
        .flatten
    // now merge sc results
    val groupHomerAnnotSums =
      scsGroupHomerAnnotSums
      .groupBy(x => x._1)
      .map((k, v) => {
        val n = v.map(x => x._3).sum
        val s = v.map(x => x._2).toVector.transpose.map(x => x.sum)
        (k +: s.toList :+ n).mkString(",")
      })
      .toList
    // output
    writeStrings2File(
        content = groupHomerAnnotSums,
        to = s"${homerAnnotOutd}/chromHMM.groupId2HomerGenomeAnnot.csv",
        overwrite = true,
        head = ("groupId" +: homerOrdElements :+ "nTotal").mkString(",")
    )

  } // end of main

} // end of object with main
