package CEMBATAC

import GRange.{GenomicRange, mouseGenomicRangeOrd}
import Bed.{BedElement4, BedPE10Element}
import Bed.BedPE10Element.readBedPE
import MetaData.TSCCMeta.projd
import SZUtils.readTable
import Peak.Peak.readBed
import Bed.BedElement4.readBed4
import DAR.DAR
import bioscala.LightCoord.GenomeCoord.{GenomeCoord, GenomeCoords}

val ATACPeakd =
  s"$projd/data/chromHMM/subclass_peak"
val ATACPeakChromHMMAnnotd =
  s"$projd/06.ChromHMM/out/CREAnnot250429"
val proximalDistalConnd =
  s"$projd/data/snATAC/sa2pdc_bedpe"
val ppdcfnm = s"$projd/data/snATAC/mba.whole.sa2subclass.pos.pdc.uniq.coacc.bedpe"
val CPMfnm = s"$projd/data/snATAC/cpm_peakBysubclass.csv"
val sa2DEd = s"$projd/12.DE/out/ATAC_sa2LRT_DE"

val wmbATACPeakfnm = s"${projd}/data/snATAC/wmb.snATAC.peak.srt.bed"

/**
  * Read CEMBA ATAC-seq CREs as List of GenomicRange
  *
  * @return List of CREs as GenomicRange
  */
def getATACPeaks: List[GenomicRange] = {
  os.read.lines
    .stream(os.Path(wmbATACPeakfnm))
    .slice(0, Int.MaxValue)
    .map(s => GenomicRange.fromStr(s, sep = "\t"))
    .toList
}

def loadAllCRE(f: String): GenomeCoords = {
  os.read.lines
    .stream(os.Path(f))
    .map(x => x.split("\t"))
    .filter(x => x.length >= 3)
    .map(x => (chrom = x(0),
      coord = (startFrom = x(1).toInt,
        endTo = x(2).toInt), strand = "+"))
    .toVector
}

/**
 * Load and extend ATAC Peaks.
 *   1. The extension will keep the original ranges if it match the
 *      original extension size, for example, 500bp. 2. The results will
 *      return a Map, chromosomes are keys, and coordinates inside the
 *      same chromosome will be sorted.
 * @param scfnm
 * @param ext
 * @return
 */
def loadandExtendATACPeak(scfnm: String,
  ext: Int = 250): Map[String, List[GenomicRange]] = {
  readTable(fnm = scfnm, sep = "\t", head = false)
    .map(x => {
      val l = x(1).toInt
      val r = x(2).toInt
      val s = l + math.floor((r - l) / 2).toInt
      new GenomicRange(chrom = x(0), startFrom = s - ext + 1,
          endTo = s + ext)
    })
    .groupBy(_.chrom)
    .map((k, v) => (k, v.sorted.toList))
}

/**
 * Load subclass-specific ATAC peaks.
 * - From data preprocessed for ChromHMM
 * - No ChrX and ChrY
 * @param sc
 *   PairedTag-format subclass name
 * @return
 * e.g, BedElement4, eg., GenomicRange("chr1", 0, 100), "chr1:0-100" 
 */
def loadscATAC(sc: String): List[BedElement4] = {
  readBed4(fnm = s"$ATACPeakd/${sc}.ATAC.peak.bed")
}

/**
  * Load ChromHMM-annotated subclass-specific ATAC peaks.
  * - No ChrX and ChrY
  * - Update 250501
  * @param sc PairedTag-format subclass name
  * @return
  * e.g., BedElement4 o GenomicRange("chr1", 0, 100), "Chr-O"
  */
def loadAnnotscATAC(sc: String): List[BedElement4] = {
  readBed4(fnm = s"$ATACPeakChromHMMAnnotd/$sc.CRE.annotBySummit.bed")
}

/**
  * Load subclass-specific cicero results of proximal-distal connections.
  * All of them pass the filtering by shuffled background. 
  * The scores are the co-accessible scores of GroupLASSO from cicero.
  * 
  * @param sc CEMBA ATAC-seq subclass name without id.
  * e.g., ABC_NN, Astro-CB_NN and so on.
  * @return
  */
def loadscpdc(sc: String): List[BedPE10Element] = {
  readBedPE(s"$proximalDistalConnd/sa2subclass.${sc}.pdc.bedpe")
}

/**
  * Get global positive proximal-distal connections from cicero.
  * 
  * The positive pearson correlations are calcualted using the ATAC-seq signals
  * of distal elements and the RNA-seq signals of putative target genes.
  * The scores are the connection scores from cicero.
  * 
  * @return
  */
def globalppdc: List[BedPE10Element] = {
  readBedPE(fnm = ppdcfnm, head = false)
}

/**
  * Load subclass-specific differential ATAC peaks using SnapATAC2. 
  * 
  * @param sc Paired-Tag-format subclass name
  * @return
  */
def loadscATACsa2DE(sc: String): List[DAR] = {
  readTable(fnm = s"$sa2DEd/${sc}_ATAC_sa2DE.tsv", sep = "\t", head = true)
    .map(x => {
      val t = x(0).split(":|-")
      val g = new GenomicRange(chrom = t(0), startFrom = t(1).toInt,
        endTo = t(2).toInt)
      new DAR(g, x(1).toDouble, x(2).toDouble, x(3).toDouble)
    })
}
