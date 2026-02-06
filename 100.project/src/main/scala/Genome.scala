package MetaData

import SZUtils.{ifelse, str2path, path2str}
import os._
import bioscala.LightCoord.GenomeCoord.isOverlapIgnoreStrand
import bioscala.LightCoord.GenomeCoord.{GenomeCoord, GenomeCoords}
import bioscala.LightCoord.GenomeCoord.fromBed

object GenomeMeta {
  val projd   = TSCCMeta.projd
  val genomed = "/projects/ps-renlab2/szu/genome"

  val mm10chromSizefnm = s"$projd/meta/mm10.chrom.sizes.lite"
  val mm10autosomeChromSizefnm =
    s"$projd/meta/mm10.autosome.chromsize.tsv"

  val mm10BlacklistBed = s"${projd}/meta/mm10-blacklist.v2.bed"

  val mm10fafnm   = s"${projd}/genome/GRCm38.p6.genome.fa"
  val mm10gff3fnm = s"${genomed}/gencode.vM25.annotation.gff3"
  val mm10rDHSfnm = s"${genomed}/mm10-rDHS.unfiltered.v2.bed"
  val mm10CREfnm  = s"${genomed}/mm10-cCREs.SCREEN.v4.bed"
  val CEMBATACfnm = os.Path(projd) / "data" / "snATAC" / "wmb.snATAC.peak.srt.bed"

  val phastConsfnm = s"${genomed}/mm10.60way.phastCons.bw"

  /**
    * Generate background or Null genomic regions using
    *  bedtool's shuffleBed.
    *
    * @param shuffleBed
    * @param fgbedfnm
    * @param chromsizefnm
    * @param exclbedfnm
    * @param outfnm
    * @param keepChrom
    * @param novlp
    * @param seed
    * @return exitCode 0 if sucess else not.
    */
  def genNullGenomicRanges(shuffleBed: String, fgbedfnm: String,
    chromsizefnm: String, exclbedfnm: String, outfnm: String,
    keepChrom: Boolean = true, novlp: Boolean = true,
    seed: Int = 0): Int = {
    val cmds = List(shuffleBed, "-excl", exclbedfnm,
        ifelse(keepChrom, "-chrom", ""),
        ifelse(novlp, "-noOverlapping", ""), "-i", fgbedfnm, "-g",
        chromsizefnm)
      .filter(x => x.size >= 1)
      .map(x => os.Shellable(List(x)))
    val r = os
      .proc(cmds*)
      .call(check = true, stdout = os.Path(outfnm))
    r.exitCode
  }

  def filterMouseSignalRanges(raw: GenomeCoords): GenomeCoords = {
    // filter by blacklist region
    val bl = fromBed(mm10BlacklistBed)
    val rawInbl = isOverlapIgnoreStrand(query = raw, subject = bl)

    // filter by mm10CREs
    val mm10CRE = fromBed(mm10CREfnm)
    val rawInmm10CRE = isOverlapIgnoreStrand(query = raw, subject = mm10CRE)

    // filter by rDHS
    val mm10rDHS = fromBed(mm10BlacklistBed)
    val rawInmm10rDHS = isOverlapIgnoreStrand(query = raw, subject = mm10rDHS)

    // filter by CEMBA ATAC CREs
    val CEMBACRE = fromBed(CEMBATACfnm)
    val rawInCEMBACRE = isOverlapIgnoreStrand(query = raw, subject = CEMBACRE)

    raw.zipWithIndex
      .filterNot((x, i) => rawInbl(i) || rawInmm10CRE(i) ||
        rawInmm10rDHS(i) || rawInCEMBACRE(i))
      .map(_._1)
    
  }
  
}
