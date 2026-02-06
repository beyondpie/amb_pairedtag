import Math.log2
import MetaData.readSubclassMeta
import SZUtils.{ifelse, path2str, writeStrings2File}
import bioscala.LightCoord.BedGraph.{getBigWigReader, getWigValue}
import bioscala.LightCoord.GenCode.{GFF3, readGeneFromGenCodeGFF}
import bioscala.LightCoord.GenomeCoord.{GenomeCoord, GenomeCoords}
import bioscala.LightCoord.getClosedCoord
import org.apache.commons.math3.stat.inference.TTest
import org.broad.igv.bbfile.BBFileReader
import os.*
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*

object SexDiffRegionCorrelatedGene {
  type DiffCRE = (g: GenomeCoord, log2fc: Double, p: Double,
    adjustP: Double, sc: String)
  type StatCompare = (p: Double, log2fc: Double)

  val projd =
    os.Path("/projects/ps-renlab2/szu/projects/amb_pairedtag")
  val workd    = projd / "18.sexDE"
  val outd     = workd / "out"
  val sexDEd   = outd / "20250619_sex_diff_peaks"
  val t        = TTest()
  val RNAbwd   = projd / "data" / "ptRNAbam" / "bigwig"
  val ptscMeta = readSubclassMeta().filter(_.atac)

  val allensc2ptscs =
    ptscMeta
      .map(x => (x.name.AllenIdName, x.name.PairedTagName))
      .toMap

  val mods  = Vector("H3K27me3", "H3K9me3", "H3K27ac", "H3K4me1")
  val sexes = Vector("Male", "Female")
  val gencodefnm =
    "/projects/ps-renlab2/szu/genome/gencode.vM25.annotation.gff3"

  def main(args: Array[String]): Unit = {
    // subclass, modality, sex
    val RNAbwReaders: Map[(String, String, String), BBFileReader] =
      ptscMeta
        .map(_.name.AllenIdName)
        .flatMap(allensc => {
          val ptsc = allensc2ptscs(allensc)
          mods
            .flatMap(mod => {
              sexes.map(sex => {
                val fnm: String =
                  RNAbwd / s"${ptsc}.${sex}.RPKM.s300.bw"
                ((allensc, mod, sex), getBigWigReader(fnm))
              })
            })
        })
        .toMap

    val mouseGenes: GFF3 =
      readGeneFromGenCodeGFF(gencodefnm)
        .filter(x => x.feaType == "gene")
        .filter(x => x.g.chrom != "chrM")
        .filter(x => x.g.chrom != "chrX" && x.g.chrom != "chrY")
        .filter(x => x.bioType == "protein_coding")

    val modsex2diff: Map[
        (String, String), Map[String, Vector[DiffCRE]]] =
      mods
        .flatMap(mod => {
          sexes.map(sex => {
            val r =
              loadSexDE(mod.toLowerCase, sex.toLowerCase)
            ((mod, sex), r.groupBy(_.sc))
          })
        })
        .toMap

    // get the closest genes for them
    val modsex2gene: Map[
        (String, String), Map[String, GenomeCoords]] =
      modsex2diff.map((k, v) => {
        (k,
            v.map((sc, gs) =>
              (
                  sc,
                  gs.flatMap(p => getNearGene(p.g, mouseGenes))
                    .distinct
              ))) // end of map k,v
      })
    // check if overall their expression are higher or lower than
    // the control group
    val modsex2DEcorGene: Map[
        (String, String), Map[String, StatCompare]] =
      modsex2gene.map((k, v) => {
        val sex    = k._2
        val oppSex = ifelse(sex == "Male", "Female", "Male")
        val r      =
          v.filter((sc, gs) => gs.length > 1)
            .filter((sc, gs) => {
              RNAbwReaders.contains((sc, k._1, sex)) &&
              RNAbwReaders.contains((sc, k._1, oppSex))
            })
            .map((sc, gs) => {
              // TODO: use average by window size
              val x = gs.map(g =>
                // val ptsc = allensc2ptscs(sc)
                getWigValue(RNAbwReaders((sc, k._1, sex)), g)
                  .map(x => x.s)
                  .sum)
              val y = gs.map(g =>
                getWigValue(RNAbwReaders((sc, k._1, oppSex)), g)
                  .map(x => x.s)
                  .sum)
              val p = compareUsingPairedTtest(x, y, t)
              (sc, p)
            })
        (k, r)
      })

    modsex2DEcorGene.foreach((modsex, sc2stats) => {
      val mod                     = modsex._1
      val sex                     = modsex._2
      val content: Vector[String] = ptscMeta
        .map(_.name.AllenIdName)
        .filter(sc => sc2stats.contains(sc))
        .map(sc => {
          sc + "," + sc2stats(sc).p.toString + "," + sc2stats(
              sc).log2fc.toString
        })
      writeStrings2File(
          content = content,
          to = outd / s"$mod.$sex.closestGene.PairedTtest.csv",
          head = Vector("AllenIdSubclass", "pValue",
              "log2foldchange").mkString(","),
          overwrite = true
      )
    }) // end of foreach

  } // end of main

  def loadSexDE(mod: String, sex: String): Vector[DiffCRE] = {
    os.read.lines
      .stream(sexDEd /
        s"20250609_${mod.toLowerCase}_${sex.toLowerCase}_enriched_regions.txt")
      .drop(1)
      .map(x => x.strip().split(","))
      .filter(_.nonEmpty)
      .map(x => {
        val g      = x(1)
        val chrom  = g.split(":").head
        val region =
          g.split(":").last.split("-").map(x => x.toInt)
        (
            g = (chrom = chrom,
                coord = (startFrom = region.head,
                    endTo = region.last), strand = ","),
            log2fc = x(2).toDouble,
            p = x(3).toDouble,
            adjustP = x(4).toDouble,
            sc = x.last
        )
      })
      .toVector
  }

  def getNearGene(x: GenomeCoord,
    gs: GFF3): Vector[GenomeCoord] = {
    val y = gs.filter(g => g.g.chrom == x.chrom)
    if (y.isEmpty) {
      Vector[GenomeCoord]()
    } else {
      val c = getClosedCoord(x.coord, y.map(_.g.coord))
      if (c.isEmpty) {
        Vector[GenomeCoord]()
      } else {
        c.map(t => (chrom = x.chrom, coord = t, strand = "."))
      }
    }
  }

  def compareUsingPairedTtest(x: Vector[Double],
    y: Vector[Double], t: TTest): StatCompare = {
    val p = t.pairedTTest(x.map(t => log2(t + 1e-6)).toArray,
        y.map(t => log2(t + 1e-6)).toArray)
    val log2fc =
      x
        .zip(y)
        .map((i, j) => log2((i + 1e-6) / (j + 1e-6)))
        .sum / x.length.toDouble
    (p, log2fc)
  }

} // end of object
