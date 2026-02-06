/**
 * Perform Sequence conservations for different ChromHMM States. Only
 * autosomes are considered.
 */

import os._
import CEMBATAC.getATACPeaks
import MetaData.GenomeMeta.genNullGenomicRanges
import MetaData.readSubclassMeta
import SZUtils.writeStrings2File
import PairedTagChromHMM.readCREChromHMMAnnot
import GRange.mouseGenomicRangeOrd
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*

object ChromHMMSequenceConservation {
  def deepToolComputeMatrix(bin: String, regionfnm: String, bwfnm: String, gzmatfnm: String, blacklistfnm: String, ncore: Int = 1): Int = {
    val cmds = List(
        bin,
        "scale-regions",
        "--regionsFileName",
        regionfnm,
        "--scoreFileName",
        bwfnm,
        "-m",
        "500",
        "--missingDataAsZero",
        "--outFileName",
        gzmatfnm,
        "-p",
        ncore.toString
    ).map(x => os.Shellable(List(x)))
    val r = os
      .proc(cmds*)
      .call(check = true)
    r.exitCode
  }

  def main(args: Array[String]) = {
    val projd = MetaData.TSCCMeta.projd
    val workd = projd + "/" + "06.ChromHMM/"
    val outd  = s"$workd/out"
    // 1. generate random sequences as background.
    val shuffleBed =
      "/home/szu/mambaforge/envs/sa2stable/bin/shuffleBed"
    val rDHSfnm        = MetaData.GenomeMeta.mm10rDHSfnm
    val ATACREfnm      = CEMBATAC.wmbATACPeakfnm
    val autosomeCREfnm = s"$outd/ATAC.autosome.CRE.bed"
    val autosomeCREs =
      getATACPeaks
        .filter(x => x.chrom != "chrX" && (x.chrom != "chrY"))
        .map(x => x.mkString("\t"))
    writeStrings2File(
        content = autosomeCREs,
        to = autosomeCREfnm,
        head = ""
    )

    val bgfnm = s"$outd/mm10.null.background.bed"
    val gfnm  = MetaData.GenomeMeta.mm10autosomeChromSizefnm
    genNullGenomicRanges(
        shuffleBed = shuffleBed,
        fgbedfnm = autosomeCREfnm,
        chromsizefnm = gfnm,
        exclbedfnm = rDHSfnm,
        outfnm = bgfnm,
        keepChrom = true,
        novlp = true,
        seed = 0
    )

    // 2. per subclass, run bigwig mappings using phastcons bigwig
    //    and save results.
    // - load all the scs
    val ptscMeta = readSubclassMeta().filter(x => x.atac)

    // - load the chromHMM annotations of different CREs for each scs
    val sc2CRE =
      ptscMeta
        .map(x => {
          val sc = x.name.PairedTagName
          (sc, readCREChromHMMAnnot(sc))
        })
        .toMap
    // - output the CRE annotations per subclasses
    val scCREAnnotd = s"$outd/CREAnnot250429"
    // record the output filename
    val scstate2fnm =
      sc2CRE
        .map((sc, CREs) => {
          CREs
            .groupBy(x => x.name)
            .map((state, annotCREs) => {
              val fnm  = s"$scCREAnnotd/$sc.CRE.$state.annotBySummit.bed"
              ((sc, state), fnm)
            }) // end of CREs's map
        })     // end of sc2CRE's map
        .flatten
        .toMap

    // - for each scs, run computeMatrix
    val computeMatrix = "/home/szu/mambaforge/bin/computeMatrix"
    val phastConsfnm = MetaData.GenomeMeta.phastConsfnm
    val outgzmatd = s"$outd/subclassCREStateComputeMatrix"
    val blfnm = MetaData.GenomeMeta.mm10BlacklistBed

    // - for null background
    deepToolComputeMatrix(
      bin = computeMatrix,
      regionfnm = bgfnm,
      bwfnm = phastConsfnm,
      gzmatfnm = s"$outd/mm10.backgroundRanges.computeMatrix.mat.gz",
      blacklistfnm = blfnm,
      ncore = 10
    )

    scstate2fnm.par.foreach((sc_state, fnm) => {
      val (sc, state) = sc_state
      deepToolComputeMatrix(
        bin = computeMatrix,
        regionfnm = fnm,
        bwfnm = phastConsfnm,
        gzmatfnm = s"$outgzmatd/$sc.$state.mat.gz",
        blacklistfnm = blfnm,
        ncore = 1
      )
      println(s"Finish computeMatrx for $sc $state.")
    })

    // 3. per subclass, get all the ATAC PhastCons scores
    val sc2autoCREd = projd + "/data/chromHMM/subclass_peak"
    val sc2autoCREfnm: Map[String, String] = ptscMeta.map(x => {
      val sc = x.name.PairedTagName
      val fnm = s"$sc2autoCREd/$sc.ATAC.peak.bed"
      (sc, fnm)
    }).toMap
    sc2autoCREfnm.par.foreach((sc, fnm) => {
      deepToolComputeMatrix(
        bin = computeMatrix,
        regionfnm = fnm,
        bwfnm = phastConsfnm,
        gzmatfnm = s"$outgzmatd/$sc.autosome.CRE.mat.gz",
        blacklistfnm = blfnm,
        ncore = 1
      )
    }) // end of sc2autoCRE foreach
  }
}
