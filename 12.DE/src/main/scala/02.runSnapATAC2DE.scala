import os._
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*
import SZUtils.{str2path, path2str}
import SZUtils.Extensions.*
import SZUtils.SimpleCommandTask
import SZUtils.readTable
import MetaData.TSCCMeta.{projd, metad}
import AllenMetaData.Allen.subclassNamePairedTag

object RunSnapATAC2DE {
  val minpct    = 0.005
  val minlog2fc = 0.25
  val dsfold    = 2
  val groupby   = "annot.sc"

  val workd = s"${projd}/12.DE"
  val logd  = s"${workd}/log"
  val flagd = s"${workd}/flag"
  val outd  = s"${workd}/out"

  val python = "/home/szu/mambaforge/envs/sa2stable/bin/python"
  val script = s"${workd}/src/main/python/02.runSnapATAC2DE.py"
  val hs     = List("H3K9me3")

  val scs: List[String] =
    readTable(
        fnm = (os.Path(metad) / "pt.subclass2count.csv").toString,
        sep = ",",
        skipEmptyField = true,
        head = true
    ).map(x => List(x(0), x(4)))
      .filter(x => x(1).toInt > 0)
      .map(x => x(0))

  def main(args: Array[String]) = {
    val DETasks: List[SimpleCommandTask] =
      scs
        .map(g => {
          hs.map(h => {
            val gstr    = subclassNamePairedTag(g)
            val suffix  = s"pct-${minpct}_log2fc-${minlog2fc}-sa2DE"
            val name    = s"${gstr}.${h}.${suffix}"
            val logfnm  = s"${logd}/${name}.log"
            val flagfnm = s"${flagd}/${name}.done"
            val outfnm  = s"${outd}/${name}.tsv"
            SimpleCommandTask(
                commands = List(python, script, "--annfnm",
                    s"${outd}/ann.${h}.pmat.h5ad", "--g1", g,
                    "--groupby", "annot.sc", "--minpct",
                    minpct.toString(), "--minlog2fc",
                    minlog2fc.toString(), "--dsfold", dsfold.toString(),
                    "--outfnm", outfnm),
                logfnm = logfnm,
                flagfnm = flagfnm,
                skip = true,
                check = true,
                verbose = true,
                name = name
            )
          })
        })
        .flatten
        .toList

    DETasks.par.foreach { t =>
      t.run()
    }

    // check job status
    // val t =
    //   scs
    //     .map(g => {
    //       val h      = "H3K9me3"
    //       val gstr   = subclassNamePairedTag(g)
    //       val suffix = s"pct-${minpct}_log2fc-${minlog2fc}-sa2DE"
    //       val name   = s"${gstr}.${h}.${suffix}"
    //       s"${flagd}/${name}.done"
    //     })
    //     .map(x => os.exists(x))
    //     .map(t => if(t) 1 else 0).sum()

    // * for debug
    // val g       = "001 CLA-EPd-CTX Car3 Glut"
    // val h       = "H3K9me3"
    // // val m       = "t-test"
    // val gstr    = subclassNamePairedTag(g)
    // val suffix  = s"pct-${minpct}_log2fc-${minlog2fc}-sa2DE"
    // val name    = s"${gstr}.${h}.${suffix}"
    // val logfnm  = s"${logd}/${name}.log"
    // val flagfnm = s"${flagd}/${name}.done"
    // val outfnm  = s"${outd}/${name}.tsv"
    // val a = SimpleCommandTask(commands = List(python, script, "--annfnm",
    //         s"${outd}/ann.${h}.pmat.h5ad", "--g1", g, "--groupby",
    //         "annot.sc", "--minpct", minpct.toString(), "--minlog2fc",
    //         minlog2fc.toString(), "--dsfold", dsfold.toString(),
    //         "--outfnm", outfnm), logfnm = logfnm,
    //     flagfnm = flagfnm, skip = false, check = true, verbose = true,
    //     name = name)
    // a.run()
  } // end of main

}
