import os._
import scala.concurrent.Future
import scala.concurrent.ExecutionContext
import java.util.concurrent.Executors
import org.jsoup._
import scala.jdk.CollectionConverters._
import scala.concurrent.Await
import scala.concurrent.duration.Duration

object DownloadDNAMeth {
  def main(args: Array[String]) = {
    val outd =
      os.Path(MetaData.TSCCMeta.projd) / "data" / "snmC_snm3C"
    val outCHNd = s"${outd.toString}/mCH_bigwig"

    val DNAMethHyperLink =
      "https://data.nemoarchive.org/biccn/grant/u19_cemba" + "/" +
        "ecker/epigenome/cellgroup/mCseq3/mouse/processed/other"
    val doc = Jsoup
      .connect(DNAMethHyperLink)
      .get()

    // * digest the content
    val headlines = doc.select("td").asScala
    val contents =
      headlines.map(x => x.select("a")).filter(x => x.size > 0)
    // val prefix = contents(0).attr("href")
    val CGNbws =
      contents.tail
        .map(x => x.attr("href"))
        .filter(x => x.contains("CGN"))
    val CHNbws =
      contents.tail
        .map(x => x.attr("href"))
        .filter(x => x.contains("CHN"))

    // * jobs
    val ec =
      ExecutionContext.fromExecutor(Executors.newFixedThreadPool(8))
    // this will run immediately
    val getCHNbws =
      CHNbws.toList
        .map(x =>
          Future {
            os.proc("wget", DNAMethHyperLink + "/" + x, "-O",
                s"${outCHNd}/${x}")
              .call(check = true)
          }(using ec))
    // check jobs are all done
    getCHNbws.map(x => x.isCompleted).filter(x => !x)
    // check jobs exits normally
    val a =
      getCHNbws
        .map(x =>
          Await.result(x, atMost = Duration(Int.MaxValue, "second")))
        .map(x => x.asInstanceOf[CommandResult])
        .map(x => x.exitCode)
        .filter(x => x == 0)
  }
}
