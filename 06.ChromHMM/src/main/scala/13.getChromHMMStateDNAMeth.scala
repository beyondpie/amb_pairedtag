import os._
import scala.collection.mutable.ArrayBuffer
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*
import SZUtils.writeStrings2File
import SZUtils.readTable

object GetChromHMMStateDNAMeth {

  def main(args: Array[String]) = {
    val projd =
      "/tscc/projects/ps-renlab2/szu/projects/amb_pairedtag"

    val mC              = "mCH"
    val workd           = s"${projd}/06.ChromHMM"
    val DNAmed          = s"${workd}/out/DNAMeth_binSize200"
    val genomeBinStated = s"${workd}/out/m-bpeak-s18_pd-obs"
    val outd            = s"${workd}/out/chromHMM18to8_DNAMe"
    val state18to8: Map[Int, String] =
      readTable(fnm = s"${projd}/meta/chromHMM18statesAnnotation.csv",
          sep = ",", head = true)
        .map(x => (x(0).toInt, x(2)))
        .toMap

    val sc  = args(0)
    val chr = args(1)
    // scs.par.foreach { sc =>
    //   {
    println(s"extract DNA ${mC} for ${sc} ...")
    //    chrs.foreach { chr =>
    //      {
    val binState =
      readTable(fnm =
            s"${genomeBinStated}/${sc}_${chr}_s18-200bin-pd-obs.csv",
          sep = ",", head = true)
        .map(x => state18to8(x(5).toInt))
        .toArray

    val binId2DNAme =
      readTable(
          fnm = s"${DNAmed}/${sc}_${chr}_${mC}_b200.avgByCount.csv",
          sep = ",", head = true)
        .map(x => (x(0).toInt, x(1).toDouble))
        .toList

    val r2 = binId2DNAme.map((binId, s) => {
      (binState(binId), s)
    })

    val f =
      s"${DNAmed}/${sc}.${chr}.8state2${mC}.b200.avgByCount.csv"

    val h = s"chr,state,${mC}"
    writeStrings2File(
        content = r2.toList
          .map(x => List(chr, x._1, x._2.toString).mkString(",")),
        to = f,
        head = h,
        overwrite = true
    )
    //     }
    //    }
    //  }
  }

  def stat: Unit = {
    val projd = "/projects/ps-renlab2/szu/projects/amb_pairedtag"

    val mC              = "mCH"
    val workd           = s"${projd}/06.ChromHMM"
    val DNAmed          = s"${workd}/out/DNAMeth_binSize200"
    val genomeBinStated = s"${workd}/out/m-bpeak-s18_pd-obs"
    val outd            = s"${workd}/out/chromHMM18to8_DNAMe"

    val scs: List[String] =
      os.list(os.Path(DNAmed))
        .map(x => x.toString)
        .filter(x => x.contains("mCG_b200.avgByCount.csv"))
        .map(x =>
          x.toString
            .split("/")
            .last
            .split("\\.")
            .head
            .replaceAll("_chr\\d+_mCG_b200", ""))
        .distinct
        .toList

    val chrs = 1.to(19).map(i => s"chr${i}")

    val sc2stateDNAme =
      scs
        .map(sc => {
          chrs.par
            .map(chr => {
              val f =
                s"${DNAmed}/${sc}.${chr}.8state2${mC}.b200.avgByCount.csv"

              readTable(fnm = f, head = true, sep = ",")
                .map(x => (x(1), x(2).toDouble))
                .groupBy(x => x._1)
                .map((state, v) =>
                  (state, v.map(x => x._2).sum / v.length))
            })
            .flatten
            .groupBy(x => x._1)
            .map((state, v) =>
              (sc, state, v.map(x => x._2).sum / v.length))
            .toList
        })
        .flatten

    // save to file
    writeStrings2File(
        content = sc2stateDNAme.map(x => s"${x._1},${x._2},${x._3}"),
        to = s"${outd}/subclass.chromHMM18to8.${mC}.avgByCount.csv",
        overwrite = true,
        head = s"subclass,state,${mC}"
    )

    val state2DNAme =
      sc2stateDNAme
        .groupBy(_._2)
        .map((state, v) => (state, v.map(_._3).sum / v.length))

  }
}
