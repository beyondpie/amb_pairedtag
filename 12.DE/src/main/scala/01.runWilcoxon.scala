import os._
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*
import SZUtils.{str2path, path2str}
import SZUtils.Extensions.*
import SZUtils.SimpleCommandTask
import SZUtils.readTable
import MetaData.TSCCMeta.{projd, metad}
import AllenMetaData.Allen.subclassNamePairedTag

object RunSimpleDE {
  val minpct    = 0.005
  val minlog2fc = 0.1
  val dsfold    = 3
  val groupby   = "annot.sc"

  val workd = s"${projd}/12.DE"
  val logd  = s"${workd}/log"
  val flagd = s"${workd}/flag"
  val outd  = s"${workd}/out"

  val python  = "/home/szu/mambaforge/envs/sa2stable/bin/python"
  val script  = s"${workd}/src/main/python/01.runSimpleDE.py"
  val hs      = List("H3K27ac", "H3K27me3", "H3K4me1", "H3K9me3")
  val methods = List("t-test", "wilcoxon")

  val scs: List[String] =
    readTable(
        fnm = (os.Path(metad) / "pt.subclass2count.csv").toString,
        sep = ",",
        skipEmptyField = true,
        head = false
    ).map(x => List(x(0), x(1)))
      // .filter(x => x(1).toInt > 100)
      .map(x => x(0))

  def main(args: Array[String]) = {
    val DETasks: List[SimpleCommandTask] =
      scs
        .map(g =>
          {
            val r = hs.map(h => {
              methods.map(m => {
                val gstr    = subclassNamePairedTag(g)
                val suffix  = s"pct-${minpct}_log2fc-${minlog2fc}-${m}"
                val name    = s"${gstr}.${h}.${suffix}"
                val logfnm  = s"${logd}/${name}.log"
                val flagfnm = s"${flagd}/${name}.done"
                val outfnm  = s"${outd}/${name}.tsv"
                SimpleCommandTask(
                    commands = List(python, script, "--annfnm",
                        s"${outd}/ann.${h}.pmat.h5ad", "--g1", g,
                        "--groupby", "annot.sc", "--minpct",
                        minpct.toString(), "--minlog2fc",
                        minlog2fc.toString(), "--dsfold",
                        dsfold.toString(), "--outfnm", outfnm,
                        "--method", m),
                    logfnm = logfnm,
                    flagfnm = flagfnm,
                    skip = true,
                    check = true,
                    verbose = true,
                    name = name
                )
              })
            })
            r
          }.flatten)
        .flatten

    // DETasks.par.foreach { t =>
    //   t.run()
    // }

    val g = "001 CLA-EPd-CTX Car3 Glut"
    val h = "H3K27ac"
    // val m       = "t-test"
    val m       = "wilcoxon"
    val gstr    = subclassNamePairedTag(g)
    val suffix  = s"pct-${minpct}_log2fc-${minlog2fc}-${m}"
    val name    = s"${gstr}.${h}.${suffix}"
    val logfnm  = s"${logd}/${name}.log"
    val flagfnm = s"${flagd}/${name}.done"
    val outfnm  = s"${outd}/${name}.tsv"
    val a = SimpleCommandTask(commands = List(python, script, "--annfnm",
            s"${outd}/ann.${h}.pmat.h5ad", "--g1", g, "--groupby",
            "annot.sc", "--minpct", minpct.toString(), "--minlog2fc",
            minlog2fc.toString(), "--dsfold", dsfold.toString(),
            "--outfnm", outfnm, "--method", m), logfnm = logfnm,
        flagfnm = flagfnm, skip = false, check = true, verbose = true,
        name = name)
    a.run()

  } // end of main

}
