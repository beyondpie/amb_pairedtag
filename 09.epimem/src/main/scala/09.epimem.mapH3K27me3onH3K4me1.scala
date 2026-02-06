import os._
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*

object MapH3K27me3onH3Kme1forEpiMem {
  val projd = "/projects/ps-renlab2/szu/projects/amb_pairedtag"
  val ptDNAbwd = s"${projd}/data/ptDNAbam/bigwig"

  val srcd = s"${projd}/09.epimem/src/main"
  val python = "/home/szu/mambaforge/bin/python"
  val mapSignalScript = s"${srcd}/python/02.mappingBigWigOnBed.py"
  val rscript = "/home/szu/mambaforge/envs/seurat/bin/Rscript"
  val getEpiMemRegionScript = s"${srcd}/R/04.getEpiMemRegion.R"

  val repressiveMod = "H3K27me3"
  val sm = if (repressiveMod == "H3K9me3") {
    "1000"
  } else { "300" }
  val activeMod = "H3K4me1"

  val scs: List[String] =
    os.list(os.Path(ptDNAbwd))
      .map(x => x.baseName)
      .map(x => x.split("\\.")(0))
      .distinct
      .filter(x => x.split("_")(0).toInt < 500)
      .filter(x =>
        os.exists(
          os.Path(s"${ptDNAbwd}/${x}.${repressiveMod}.e100.bs100.sm${sm}.bw")
        )
      )
      .filter(x =>
        os.exists(os.Path(s"${ptDNAbwd}/${x}.${activeMod}.e100.bs100.sm300.bw"))
      )
      .filter(x => x != "021_L4_RSP_ACA_Glut")
      .toList

  def main(args: Array[String]) = {
    scs.par.foreach { sc =>
      {
        os.proc(
          python,
          mapSignalScript,
          sc,
          repressiveMod,
          activeMod
        ).call(check = true)
        println(s"Finish ${repressiveMod} on ${activeMod} for ${sc}.")
      }
    }

    scs.par.foreach { sc =>
      {
        os.proc(
          python,
          mapSignalScript,
          sc,
          activeMod,
          repressiveMod
        ).call(check = true)
        println(s"Finish ${activeMod} on ${repressiveMod} for ${sc}.")
      }
    }

    // Here we intersect inactive on active and active on inactive.
    scs.par.foreach { sc => {
        os.proc(
          rscript,
          getEpiMemRegionScript,
          sc,
          s"${repressiveMod},${activeMod}",
          s"${activeMod},${repressiveMod}",
          "cross_interact",
          "sequence"
        ).call(check = true)
        println(
          s"Finish epiMem for ${sc} using ${activeMod} and ${repressiveMod}."
        )
      }
    }
  }
}
