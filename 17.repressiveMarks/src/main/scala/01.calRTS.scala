import os._
import bioscala.LightCoord.Bed._
import bioscala.LightCoord.GenCode._
import bioscala.LightCoord.GenomeCoord._
import bioscala.LightCoord.GeneBody3to5
import bioscala.LightCoord.GeneBody5to3
import Genome.MouseGenome.chr2size
import bioscala.LightCoord.getCenter
import bioscala.TRIAGE.listAllUniqueElement
import MetaData.readSubclassMeta
import bioscala.TRIAGE.getTopDomainWidth
import bioscala.TRIAGE.getRTS
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*
import bioscala.TRIAGE.assignBroadDomain2Gene
import SZUtils.{writeStrings2File, ifelse}

@main def calRTS(): Unit = {
  // * meta
  val projd = MetaData.TSCCMeta.projd
  val workd = s"$projd/17.repressiveMarks"
  val outd  = s"$workd/out"
  val figd  = s"$workd/figure"
  val peakd =
    s"$projd/data/pairedtag_peak/merge_peak/subclassOnUnifiedPeak"
  val gencodefnm =
    "/mnt/tscc2/szu/genome/gencode.vM25.annotation.gff3"
  val ptscMeta =
    readSubclassMeta()
      .filter(x => x.atac)
  val ptscs: Vector[String] =
    ptscMeta.map(x => x.name.PairedTagName)

  val mods         = Vector("H3K27me3", "H3K9me3")
  val qtop         = 0.05
  val defaultScore = -1.0
  val namefn       = (x: GFF3Element) => x.gname
  // * functions
  def loadMergedPeak(fnm: String): LightBed = {
    os.read.lines
      .stream(os.Path(fnm))
      .drop(1)
      .map(x => x.strip().split("\t"))
      .map(x =>
        (
            g = (chrom = x(0),
                coord = (startFrom = x(1).toInt,
                    endTo = x(2).toInt), strand = "."),
            name = s"${x(0)}:${x(1)}-${x(2)}",
            score = x(2).toDouble - x(1).toDouble
        ))
      .toVector
  }

  // * process
  val genes: GFF3 = readGeneFromGenCodeGFF(gencodefnm)
      .filter(x => x.feaType == "gene")
      .filter(x => x.g.chrom != "chrM")
      .filter(x => x.g.chrom != "chrX" && x.g.chrom != "chrY")
      .filter(x => x.bioType == "protein_coding")

  val proximalGenes =
    aroundGene(g = genes, upStream = 2500, downStream = 0,
        feaType = "proximal", crd3to5 = GeneBody3to5,
        crd5to3 = GeneBody5to3)(using chr2size)

  mods.foreach(mod => {
    val outscMapDomaind = outd + s"/subclass${mod}OnGene"
    if (!os.exists(os.Path(outscMapDomaind))) {
      os.makeDir(os.Path(outscMapDomaind))
    }

    // test for one subclass
    // val sc    = ptscs(0)
    // val brdmn = loadMergedPeak(peakd + d"/$sc.H3K27me3.upeak.bed")
    // val ctrdmn = brdmn.map(x => {
    //   val g = x.g.toTuple
    //     .copy(_2 = getCenter(x.g.coord))
    //     .asInstanceOf[GenomeCoord]
    //   x.toTuple.copy(_1 = g).asInstanceOf[LightBedElement]
    // })
    // val assgnbd = assignBroadDomain2Gene(genes = proximalGenes,
    //     d = ctrdmn, rd = brdmn)
    // val a  = assgnbd.filter(x => x.score > 0)
    // val b  = proximalGenes.filter(x => x.gname == "Rp1")
    // val b0 = genes.filter(x => x.gname == "Rp1")

    // // output
    // writeStrings2File(
    //     content = assgnbd.map(x =>
    //       mkStringLightBedElement(x, sep = "\t")),
    //     to = outscMapDomaind + d"/$sc.H3K27me3.broadDomain.bed",
    //     overwrite = true,
    //     head = ""
    // )

    // run for all subclasses
    // ptscs.par
    //   .foreach(sc => {
    //     val brdmn =
    //       loadMergedPeak(peakd + d"/$sc.${mod}.upeak.bed")
    //     val ctrdmn = brdmn.map(x => {
    //       val g = x.g.toTuple
    //         .copy(_2 = getCenter(x.g.coord))
    //         .asInstanceOf[GenomeCoord]
    //       x.toTuple.copy(_1 = g).asInstanceOf[LightBedElement]
    //     })
    //     val assgnbd = assignBroadDomain2Gene(
    //         genes = proximalGenes,
    //         d = ctrdmn,
    //         rd = brdmn,
    //         namefn = namefn,
    //         defaultScore = defaultScore
    //     )
    //     writeStrings2File(
    //         content = assgnbd.map(x =>
    //           mkStringLightBedElement(x, sep = "\t")),
    //         to = outscMapDomaind + d"/$sc.${mod}.broadDomain.bed",
    //         overwrite = true,
    //         head = ""
    //     )
    //     println(d"map ${mod} broad domain for $sc done.")
    //   }) // end of ptscs foreach of mapping broad domain

    val assgnbds: Vector[LightBed] = ptscs.map(sc =>
      readLightBed(
          fnm = outscMapDomaind + s"/$sc.${mod}.broadDomain.bed"))
    val nGene     = assgnbds(0).length
    val nCelltype = assgnbds.length
    println(
        s"Before filtering low overlapping genes, $nGene genes in total.")
    val rmgid: Seq[Int] = Seq
      .range(0, nGene)
      .map(i => {
        val r =
          assgnbds.map(x => (x(i).score - defaultScore).abs).sum
        ifelse(r < 1e-10, i.toInt, -1.toInt)
      })
      .filter(x => x > 0)

    // filter broad domains
    val bd = if (rmgid.nonEmpty) {
      println(
          s"${rmgid.length} genes to move due to low overlapping.")
      assgnbds.map(x =>
        x.zipWithIndex
          .filter((x, i) => !rmgid.contains(i))
          .map(_._1))
    } else {
      assgnbds
    }

    // filter genes
    val filterGenes = if (rmgid.nonEmpty) {
      Seq
        .range(0, nGene)
        .filter(i => !rmgid.contains(i))
        .map(i => genes(i))
    } else {
      genes
    }

    val width = getTopDomainWidth(x = bd, qtop = qtop,
        defaultScore = defaultScore)
    println(s"top ${qtop} domain width is: $width.")
    val rts: Vector[Double] =
      getRTS(x = bd, width, defaultScore = defaultScore)

    val r: LightBed = filterGenes.toVector
      .zip(rts)
      .map((x, s) =>
        (
            (g = x.g,
                name = namefn(x),
                score = s)
        ))
    val topGenes = r.sortBy(x => -x.score)

    // save results
    writeStrings2File(
        content = topGenes.map(x =>
          mkStringLightBedElement(x, sep = "\t")),
        to = outd + s"/gene2RTS.${mod}.top${qtop}.all.bed",
        overwrite = true,
        head = ""
    )

    // record process output
    /** Before filtering low overlapping genes, 20741
      * genes in total. 2079 genes to move due to low
      * overlapping. top 0.05 domain width is: 15237.0.
      * Acroo 151 celltypes: meanDomainWidth:
      * 15303.708299952994. meanStd: 19648.481648325975
      * done for H3K27me3.
      *
      * Before filtering low overlapping genes, 20741
      * genes in total. 2702 genes to move due to low
      * overlapping. top 0.05 domain width is:
      * 5592.0.Acroo 151 celltypes meanDomainWidth:
      * 6274.4058562249065 meanStd: 9925.542628878999
      * done for H3K9me3.
      */

    println(s"done for $mod.")
  }) // end of mods foreach
}
