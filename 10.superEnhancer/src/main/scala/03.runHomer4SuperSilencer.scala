import os._
import Homer2.MotifFinderByHomer2
import SZUtils.{str2path, path2str, fromString2PathRedirect}
import SZUtils.fromString2Option
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*
import MetaData.readSubclassMeta


@main def Homer4SuperSilencer(): Unit = {
  val projd = MetaData.TSCCMeta.projd
  val workd = os.Path(projd)/"10.superEnhancer"
  val silencerd = workd/"out"/"ROSE_H3K27me3_bed"
  val outd = workd/"out"/"homer4superSilencer"
  val logd = workd/"log"
  val flagd = workd/"flag"
  val ptscMeta       = readSubclassMeta().filter(_.atac)
  val ptscs = ptscMeta.map(_.name.PairedTagName).filter(
    sc => os.exists(silencerd/s"$sc.superEhancers.H3K27me3.bed")
  )

  val sc2t: Map[String, MotifFinderByHomer2] = ptscs.map(sc => {
    val name = s"$sc.superSilencer"
    val t = new MotifFinderByHomer2(
      inputfnm = silencerd/s"$sc.superEhancers.H3K27me3.bed",
      flagfnm = flagd/s"$name.homer.done",
      logfnm = logd/s"$name.homer.log",
      name = name,
      outd = outd/name,
      skip = false,
      check = true,
      homer2 = "/projects/ps-renlab2/szu/softwares/homer/bin",
      seqsize = -1,
      bgfnm = None,
      denovo = "-nomotif",
      mask = "-mask",
      hyperGenomic = "",
      genome = "mm10",
      ncore = 1
    )
    (sc, t)
  }).toMap

  sc2t.par.foreach((sc, t) => {
    println(s"run homer for superSilencer of $sc.")
    t.run()
    println(s"run homer for superSilencer of $sc done.")
  })

  // * test
  // val sc = ptscs(0)
  // val t = sc2t(sc)
  // t.run()
}


