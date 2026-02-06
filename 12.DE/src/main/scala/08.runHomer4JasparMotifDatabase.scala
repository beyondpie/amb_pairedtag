import os._
import Homer2.MotifFinderByHomer2
import SZUtils.{str2path, path2str, fromString2PathRedirect}
import SZUtils.fromString2Option
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*

@main def RunJasparMotifHomer(): Unit = {
  // for DAR CREs
  // * meta
  val projd = MetaData.TSCCMeta.projd
  val workd = projd + "/12.DE"
  val flagd = workd + "/flag"
  val logd  = workd + "/log"
  val outd  = os.Path(workd) / "out" / "JasparMotifDB"
  if (!os.exists(outd)) {
    os.makeDir(outd)
  }

  val chromStates = Vector("ChrA", "ChrO", "ChrP")
  val ptscMeta = MetaData.readSubclassMeta().filter(x => x.atac)
  val rptscs = MetaData.representPaiedTagSubclass.toSet
  val ptscs = ptscMeta.map(x => x.name.PairedTagName)
    .filter(x => !rptscs.contains(x))


  val scState2CREfnm: Map[(String, String), String] = ptscs.flatMap(sc => {
    chromStates.map(s => {
      val fnm = os.Path(workd)/"out"/s"DAR$s"/s"$sc.DAR$s.bed"
      ((sc, s), fnm.toString)
    })
  }).toMap

  val scState2HomerTask: Map[(String, String), MotifFinderByHomer2] = scState2CREfnm.map(
    (k, fnm) => {
      val sc = k._1
      val s = k._2
      val name = s"$sc.DAR.$s"
      val skip: Boolean = if(os.exists(outd/name/"knownResults.txt")) {
        println(s"$name DONOT need to run homer.")
        true
      } else {
        println(s"$name NEED to run homer.")
        false
      }
      val t = new MotifFinderByHomer2(
        inputfnm = fnm,
        flagfnm = os.Path(flagd)/s"$name.jaspar.done",
        logfnm = os.Path(logd)/s"$name.jaspar.log",
        name = name,
        outd = outd/name,
        skip = skip,
        check = false,
        homer2 =
          "/projects/ps-renlab2/szu/softwares/homer/bin",
        seqsize = -1,
        bgfnm = None,
        denovo = "-nomotif",
        mask = "-mask",
        hyperGenomic = "",
        genome = "mm10",
        mknown =
          (os.Path(projd)/"meta"/"jaspar2022_homer.motif").toString,
        ncore = 1
      )
      println(s"$name homer DONE.")
      (k, t)
    })

  scState2HomerTask.par.foreach((k, v) => {
    println(
        s"start ${k._1} ${k._2} homer using jaspar motif database.")
    v.run()
    println(s"${k._1} ${k._2} homer done.")
  })

  println("done.")
}
