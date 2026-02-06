import MergePeak.BlacklistFilter
import MergePeak.GenomeFilter
import MetaData.GenomeMeta.mm10BlacklistBed
import MetaData.GenomeMeta.mm10fafnm
import MetaData.readSubclassMeta
import SZUtils.ifelse
import SZUtils.writeStrings2File
import bioscala.LightCoord.Bed.LightBed
import bioscala.LightCoord.GenomeCoord.*
import bioscala.LightCoord.GenomeCoord.genomeCoordIgnoreStrandOrd
import os.*
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*

@main def getUnifiedPeakPerSubclass: Unit = {
  // * meta
  val projd       = MetaData.TSCCMeta.projd
  val peakd       = projd + "/data/pairedtag_peak"
  val mergedPeakd = peakd + "/merge_peak"
  val scpeakd     = peakd + "/subclass_peak"
  val ptscMeta    = readSubclassMeta().filter(_.atac)
  val ptscs       = ptscMeta.map(x => x.name.PairedTagName)

  val Histone2UnifedPeakfnm: Map[String, String] = Map(
      "H3K27ac" -> mergedPeakd.+(
          "/H3K27ac.merged.all.blv2.me.bestspm.peak"),
      "H3K27me3" -> mergedPeakd.+(
          "/H3K27me3.merged.all.blv2.me.peak"),
      "H3K4me1" -> mergedPeakd.+(
          "/H3K4me1.merged.all.blv2.me.peak"),
      "H3K9me3" -> mergedPeakd.+(
          "/H3K9me3.merged.all.blv2.me.peak")
  )

  val outd = mergedPeakd + "/subclassOnUnifiedPeak"

  // * functions
  def loadPeak(fnm: String): LightBed = {
    try {
      os.read.lines
        .stream(os.Path(fnm))
        .drop(1)
        .map(x => x.strip().split("\t"))
        .filter(x => x.length >= 4) // some files has one more empty line.
        .map(x =>
          (
              g = (chrom = x(0), coord = (x(1).toInt, x(2).toInt),
                  strand = "."),
              name = x(3),
              score = x.last.toDouble
          ))
        .toVector
    } catch {
      case e: java.lang.ArrayIndexOutOfBoundsException => {
        println(s"$fnm has no peaks.")
        Vector((
                g = (chrom = "chr1", coord = (0, 0),
                    strand = "."),
                name = "empty",
                score = 0.0
            ))
      }
    }
  }

  def loadSubclassModPeak(sc: String, mod: String): LightBed = {
    val fnm = mod match {
      case "H3K27ac" => scpeakd + s"/${sc}-${mod}.BestSPM.peak"
      case _ => scpeakd + s"/${sc}-${mod}.bedtoolmerge.peak"
    }
    loadPeak(fnm)
  }

  def loadUnifiedPeak(mod: String): LightBed = {
    loadPeak(Histone2UnifedPeakfnm(mod))
  }

  def map2UnifiedPeak(query: LightBed,
    subject: LightBed): GenomeCoords = {
    val r = findSubjectOvlpGenomeCoordIgnoreStrand(
        query.map(x => x.g), subject.map(x => x.g))
      .filter(x => x.coord.endTo != 0)
    if (r.nonEmpty) {
      r.distinct.sorted(using genomeCoordIgnoreStrandOrd)
    } else {
      r
    }
  }

  def save(x: GenomeCoords, sc: String, mod: String): Unit = {
    val r = x.map(x =>
      mkStringGenomeCoord(x, sep = "\t", useStrand = false))
    val fnm = outd + s"/$sc.$mod.upeak.bed"
    writeStrings2File(content = r, to = fnm, overwrite = true,
        head = "")
  }

  // * process
  Histone2UnifedPeakfnm.foreach((mod, ufnm) => {
    val subject = loadUnifiedPeak(mod)
    ptscs.par.foreach(sc => {
      val fnm = outd + s"/$sc.$mod.upeak.bed"
      if (!os.exists(os.Path(fnm))) {
        println(s"working on: $sc $mod")
        val query =
          loadSubclassModPeak(sc = sc, mod = mod)
            .filter(x =>
              x.g.chrom != "chrX" && x.g.chrom != "chrY")
            .filter(x => x.g.coord.endTo > 0)
        if (query.isEmpty) {
          println(s"$sc $mod has no peak when loading, skip.")
        } else {
          val r = map2UnifiedPeak(query, subject)
          if (r.isEmpty) {
            println(
                s"$sc $mod has no result after map2unifiedPeak, skip.")
          } else {
            save(r, sc, mod)
          }
        }
        println(s"$sc $mod done.")
      }
    })
  })

  // * test
  // val mod = "H3K27ac"
  // val ufnm = Histone2UnifedPeakfnm(mod)
  // val sc = ptscs.head
}
