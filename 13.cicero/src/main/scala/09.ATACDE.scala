// analyze ATAC-seq DE results with chromHMM annot.
/*
 - Overlap between DE ATAC-peaks with chromHMM annots.
   - Draw heatmap of chr-A DE peaks.
 - for each subclass, we perform state-specific motif analysis.
   - for the TFs, check their expressions using Paired-Tag data.
 - using our previous cicero results, we check
   - if Chr-A is significant enriched for cicero results
   - The heatmap from cicero and correlated genes (using Paired-Tag).
 */
import os._

object VisualATACDE{
  case class EpiDE(name: String, log2fd: Double, pval: Double, qval: Double)
  def loadATACde(fnm:String) = ???

  def main(args: Array[String]) = {
    val projd = "/projects/ps-renlab2/szu/projects/amb_pairedtag"
    val workd = s"${projd}/13.cicero"
    val ATACded = s"${projd}/12.DE/out/ATAC_sa2LRT_DE"

    println("hi")
  }
}
