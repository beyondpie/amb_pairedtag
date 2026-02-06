import MetaData.readSubclassMeta
import SZUtils.{ifelse, path2str, writeStrings2File}
import os.*
import smile.math.MathEx
import scala.util.Random
import scala.math.{pow, sqrt}
import Math.d2
import scala.collection.parallel.*
import scala.collection.parallel.CollectionConverters.*

@main def CalscSimUsingSingleCellEmbed(): Unit = {
  type SingleCellEmbed = (
    barcode: String,
    sc: String,
    annotC: String,
    spectralEmbed: Vector[Double],
    umap: Vector[Double]
  )

  val projd =
    os.Path("/projects/ps-renlab2/szu/projects/amb_pairedtag")
  val workd        = projd / "17.repressiveMarks"
  val outd         = workd / "out" / "allGenomicRange"
  val embed2       = outd
  val nsample: Int = 5000
  val ptscMeta     = readSubclassMeta(f = projd / "meta" /
    "PairedTagSubclassMetaFromCellMetaChromHMM.csv").filter(x =>
    x.atac)
  val allenscs = ptscMeta.map(x => x.name.AllenIdName)

  def readSingleCellEmbed(mod: String): Vector[SingleCellEmbed] = {
    val fnm = outd / s"pt.singlecell.embed.$mod.csv"
    os.read.lines
      .stream(fnm)
      .drop(1)
      .map(x => x.split(","))
      .filter(_.nonEmpty)
      .map(x => {
        val se =
          x.drop(3).dropRight(2).map(y => y.toDouble).toVector
        val umap = x.takeRight(2).map(y => y.toDouble)
        (barcode = x.head, sc = x(1), annotC = x(2),
            spectralEmbed = se, umap = Vector(umap.head, umap.last))
      })
      .toVector
  }

  def getGroupEmbed(x: Vector[SingleCellEmbed], groupBy: String = "sc"): Vector[SingleCellEmbed] = {
    x.groupBy(x => {
      ifelse(groupBy == "sc", x.sc, x.annotC)
    }).map((k, xs) => {
      val se = xs
        .map(x => x.spectralEmbed)
        .transpose
        .map(x => x.sum / x.length.toDouble)
      val umap1 = xs.map(x => x.umap.head).sum / xs.length.toDouble
      val umap2 = xs.map(x => x.umap.last).sum / xs.length.toDouble
      (barcode = k, sc = xs.head.sc, annotC = xs.head.annotC,
          spectralEmbed = se, umap = Vector(umap1, umap2))
    }).toVector
  }

  def getGroupPairwiseSim(ge: Vector[SingleCellEmbed], usingEmbed: String = "umap"): Vector[Vector[Double]] = {
    val scs  = ge.map(_.barcode)
    val sc2v = ge.map(x => {
      if (usingEmbed == "umap") {
        Vector(x.umap.head, x.umap.last)
      } else {
        x.spectralEmbed
      }
    })
    val r = scs.indices
      .combinations(2)
      .map(xs => {
        val i = xs.head
        val j = xs.last
        val s = MathEx.cor(sc2v(i).toArray, sc2v(j).toArray)
        ((scs(i), scs(j)), s)
      })
      .toMap

    scs.map(sc1 => {
      scs.map(sc2 => {
        if (sc1 == sc2) {
          1.0
        } else if (r.contains((sc1, sc2))) {
          r((sc1, sc2))
        } else {
          r((sc2, sc1))
        }
      })
    })
  }

  def getWithinGroupDist(se: Vector[SingleCellEmbed], usingEmbed: String = "Umap", groupBy: String = "sc", nsample: Int = 1000): Vector[(String, Double)] = {
    se.groupBy(x => { ifelse(groupBy == "sc", x.sc, x.annotC) })
      .map((g, xs) => {
        val index    = Random.shuffle(xs.indices).toVector.combinations(2)
        val rdmIndex =
          index.take(nsample).toVector
        val xy = if (usingEmbed == "Umap") {
          rdmIndex.map(t =>
            (
                (xs(t.head).umap,
                    xs(t.last).umap)
            ))
        } else {
          rdmIndex.map(t => {
            (xs(t.head).spectralEmbed, xs(t.last).spectralEmbed)
          })
        }

        val ds: Vector[Double] =
          xy
            .toVector
            .map(lr => {
              d2(lr._1, lr._2)
            })
            .toVector
        
        (g, ds.sum / ds.length.toDouble)
      })
      .toVector
  }

  def getBetweenGroupDist(se: Vector[SingleCellEmbed], groups: Vector[String], usingEmbed: String = "Umap", groupBy: String = "sc", nsample: Int = 1000): Vector[(String, String, Double)] = {
    val allG =
      se.groupBy(x => ifelse(groupBy == "sc", x.sc, x.annotC))
    val mygroups = if (groups.nonEmpty){
      groups
    } else {
      allG.keys.toVector
    }
    mygroups.indices
      .combinations(2)
      .map(x => {
        val g1 = mygroups(x.head)
        val g2 = mygroups(x.last)
        val se1 = if (usingEmbed == "Umap") {
          Random
            .shuffle(allG(g1))
            .take(nsample)
            .map(x => x.umap)
            .toVector
        } else {
          Random
            .shuffle(allG(g1))
            .take(nsample)
            .map(x => x.spectralEmbed)
            .toVector

        }
        val se2 = if (usingEmbed == "Umap") {
          Random
            .shuffle(allG(g2))
            .take(nsample)
            .map(x => x.umap)
            .toVector
        } else {
          Random
            .shuffle(allG(g2))
            .take(nsample)
            .map(x => x.spectralEmbed)
            .toVector

        }
        val ds = se1.zip(se2).map((x, y) => d2(x, y))
        (g1, g2, ds.sum / ds.length.toDouble)
      })
      .toVector
  }

  // 1. calc similarity using pseduco-bulk
  /*  List("sc", "class").map(group => {
    List("H3K27me3", "H3K9me3").map(mod => {
      List("SpectralEmbed").map(embed => {

        // val group = "sc"
        // val mod = "H3K27me3"
        // val embed = "SpectralEmbed"

        val singleCellEmbed =
          readSingleCellEmbed(mod = mod)
        val groupEmbed =
          getGroupEmbed(singleCellEmbed, groupBy = group)
        val scSim =
          getGroupPairwiseSim(groupEmbed, usingEmbed = embed)
        // save result
        writeStrings2File(
            content = scSim.map(x => x.mkString(",")),
            to = outd / s"pt.$group.$mod.cor.$embed.csv",
            overwrite = true,
            head = groupEmbed.map(_.barcode).mkString(sep = ",")
        )
      })
    })
  })*/

  // 2. calc similarity using single-cell level
  val seK27me3 = readSingleCellEmbed("H3K27me3")
  val seK9me3  = readSingleCellEmbed("H3K9me3")
  List("sc", "class").par.foreach(group => {
    List("H3K27me3", "H3K9me3").par.foreach(mod => {
      List("Umap", "SpectralEmbed").par.foreach(embed => {
        // val group = "sc"
        // val mod   = "H3K27me3"
        // val embed = "SpectralEmbed"

        println(s"cal between distances for $group $mod $embed.")
        val groups = if(group == "sc"){
          allenscs
        } else {
          Vector[String]()
        }
        val distBetweenGroup =
          if (mod == "H3K27me3") {
            getBetweenGroupDist(seK27me3, groups, embed, group,
                nsample)
          } else {
            getBetweenGroupDist(seK9me3, groups, embed, group,
                nsample)
          }
        writeStrings2File(
            content = distBetweenGroup.map(x =>
              x._1 + "," + x._2 + "," + x._3.toString),
            to =
              outd / s"pt.$group.$mod.betweengroup.dist.$embed.csv",
            overwrite = true,
            head = ""
        )

        println(s"cal within distances for $group $mod $embed.")
        val distWithinGroup =
          if (mod == "H3K27me3") {
            getWithinGroupDist(seK27me3, embed, group, nsample)
          } else {
            getWithinGroupDist(seK9me3, embed, group, nsample)
          }
        writeStrings2File(
            content = distWithinGroup.map(x =>
              x._1 + "," + x._2.toString),
            to =
              outd / s"pt.$group.$mod.withingroup.dist.$embed.csv",
            overwrite = true,
            head = ""
        )

        println(
            s"Done: cal between distances for $group $mod $embed.")
      })
    })
  })
} // end of main function
