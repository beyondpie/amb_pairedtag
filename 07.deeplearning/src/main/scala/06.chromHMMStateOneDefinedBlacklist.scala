import os._
import MetaData.readSubclassMeta
import PairedTagChromHMM.genomicRegionAnnotd

object genChromHMMStateOneDefinedBlacklist {
  val blacklistState: String = "1"
  type LightGenomicRange = (chrom: String, startFrom: Int, endTo: Int)
  def lightLoadDenseBed(fnm: String): Vector[LightGenomicRange] = {
    os.read.lines.stream(os.Path(fnm))
      .map(x => x.split("\t").take(4))
      .filter(x => x(3) == blacklistState)
      .map(x => (chrom=x(0), startFrom=x(1).toInt, endTo=x(2).toInt))
      .toVector
  }
  
  def main(args: Array[String]) = {
    val ptscMeta = readSubclassMeta().filter(_.atac)

    // 1. load ChromHMM State One annotated Genomic Ranges
    val sc2blstate = ptscMeta.map(x => {
      val fnm = s"$genomicRegionAnnotd/${x.name.PairedTagName}"
      lightLoadDenseBed(fnm)
    })
    
    // 2. As long as a region happens in at least 80% of the subclasses
    //    mark it as blacklist.
    
    // 3. merge all the blacklist into dense bed
    // 4. output the result
  }
}
