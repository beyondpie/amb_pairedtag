package MetaData

/*
 TODO:
 - add ATAC-seq subclass number into meta.
   check [projd]/meta/snATAC.subclass2cnt.csv for details.
 - move the codes about generating the meta files into 00.dataprocess
 
 */

import os._
import AllenMetaData.Allen.subclassIdNamePairedTag2ATAC
import AllenMetaData.Allen.subclassIdNamePairedTag2DNAMeth
import SZUtils.{readTable, writeStrings2File}
import MetaData.TSCCMeta
import PairedTagChromHMM.cellMarkFileTable
import AllenMetaData.Allen.subclassNameCEMBATAC
import AllenMetaData.Allen.subclassIdNamePairedTag2raw
import AllenMetaData.Allen.subclassNamePairedTag
import AllenMetaData.Allen.subclassNameDNAMeth

def bool2str(x: Boolean): String = x match {case true  => "1"
  case false => "0"
}

case class NumOfCell(name: String, male: Int, female: Int,
  maleRepA: Int, maleRepB: Int, femaleRepA: Int, femaleRepB: Int) {
  def mkString(sep: String = "-"): String = {
    name + sep +
      s"male${male}:${maleRepA}:${maleRepB}" + "_" +
      s"female${female}:${femaleRepA}:${femaleRepB}"
  }
}

def genNumOfCell(rep2num: Map[String, Int], name: String): NumOfCell = {
  val n =
    List("MaleA", "MaleB", "FemaleA", "FemaleB").map(x =>
      rep2num.getOrElse(x, 0))
  new NumOfCell(
      name = name,
      male = n(0) + n(1),
      female = n(2) + n(3),
      maleRepA = n(0),
      maleRepB = n(1),
      femaleRepA = n(2),
      femaleRepB = n(3)
  )
}

/**
 * Construct NumOfCell object from a formatted string.
 *
 * @param s
 *   e.g., "H3K27me3-male204:107:97_female152:105:47"
 * @param sep
 *   default "-", which seprates the name from the body.
 * @return
 */
def str2NumOfCell(s: String, sep: String = "-"): NumOfCell = {
  val r  = s.split(sep)
  val ss = r(1).split("_")
  val m  = ss(0).replace("male", "").split(":")
  val f  = ss(1).replace("female", "").split(":")
  new NumOfCell(
      name = r(0),
      male = m(0).toInt,
      female = f(0).toInt,
      maleRepA = m(1).toInt,
      maleRepB = m(2).toInt,
      femaleRepA = f(1).toInt,
      femaleRepB = f(2).toInt
  )
}

case class ClusterName(PairedTagName: String, ATACName: String,
  DNAmethName: String, AllenIdName: String) {
  def mkString(sep: String = ","): String = {
    List(PairedTagName, ATACName, DNAmethName, AllenIdName).mkString(
        sep)
  }
  override def toString(): String = {
    mkString(sep = ",")
  }
}

def getClusterNameHeader(sep: String = ","): String = {
  List("PairedTagName", "ATACName", "DNAmethName", "AllenIdName")
    .mkString(sep)
}

case class HasHistoneModTag(ifany: Boolean, ifmale: Boolean,
  ifemale: Boolean) {
  def mkString(sep: String = "-"): String = {
    List(bool2str(ifany), bool2str(ifmale), bool2str(ifemale))
      .mkString(sep)
  }
  override def toString(): String = {
    mkString(",")
  }
}

def toHasHistoneModTag(x: String,
  sep: String = "-"): HasHistoneModTag = {
  val r = x.split(sep).map(x => (x.toInt > 0))
  new HasHistoneModTag(ifany = r(0), ifmale = r(1), ifemale = r(2))
}

case class PairedTagClusterChromHMM(
  name: ClusterName, hasH3K27ac: HasHistoneModTag,
  hasH3K27me3: HasHistoneModTag, hasH3K4me1: HasHistoneModTag,
  hasH3K9me3: HasHistoneModTag, hasATAC: Boolean, hasDNAMeth: Boolean,
  hasHiC: Boolean
) {
  def mkString(sep: String = ",", insep: String = "-"): String = {
    List(name.mkString(sep), hasH3K27ac.mkString(insep),
        hasH3K27me3.mkString(insep), hasH3K4me1.mkString(insep),
        hasH3K9me3.mkString(insep), bool2str(hasATAC),
        bool2str(hasDNAMeth), bool2str(hasHiC)).mkString(sep)
  }
}

case class PairedTagCluster(
  name: ClusterName, all: NumOfCell, h3K4me1: NumOfCell,
  h3K9me3: NumOfCell, h3K27ac: NumOfCell, h3K27me3: NumOfCell,
  atac: Boolean, meDNA: Boolean, hiC: Boolean
)

def getPairedTagClusterHeader(sep: String = ","): String = {
  List(getClusterNameHeader(sep), "H3K27ac", "H3K27me3", "H3K4me1",
      "H3K9me3", "ATAC", "DNAMeth", "HiC").mkString(sep)
}

def initSubclassMetaFromChromHMMCellMarkFileTable() = {
  val projd = MetaData.TSCCMeta.projd
  val sep   = ","

  val scIdNamePt2atac = subclassIdNamePairedTag2ATAC
  val scIdNamePt2me   = subclassIdNamePairedTag2DNAMeth
  val scIdNamePt2raw  = subclassIdNamePairedTag2raw

  val r =
    readTable(cellMarkFileTable, sep = "\t", head = false)
      .map(x => x(0))
      .distinct
      .map(x =>
        new PairedTagClusterChromHMM(
            name = new ClusterName(x, scIdNamePt2atac(x),
                scIdNamePt2me(x), scIdNamePt2raw(x)),
            hasH3K27ac = HasHistoneModTag(true, true, true),
            hasH3K27me3 = HasHistoneModTag(true, true, true),
            hasH3K4me1 = HasHistoneModTag(true, true, true),
            hasH3K9me3 = HasHistoneModTag(true, true, true),
            hasATAC = true,
            hasDNAMeth = true,
            hasHiC = true
        ))

  // save to meta
  val content = r.map(x => x.mkString(sep))
  writeStrings2File(
      content = content,
      to = s"${projd}/meta/PairedTagSubclassMeta.csv",
      overwrite = true,
      head = getPairedTagClusterHeader(sep)
  )
} // end of getSubclassMeta

/**
 * Initiate the number of cells inforation for all the subclasses
 * annotated in PairedTag Cell Meta.
 *
 * @param sep1
 * @param sep2
 * @return
 *   Tuple of subclasses (ordered) and the number of cells The latter
 *   one is like: all-male400:201:199_female400:200:200,H3K27ac-male....
 */
def initSubclassMetaFromCellMeta(sep1: String = ",",
  sep2: String = "-"): Vector[(String, String)] = {
  val mods = "all" +: TSCCMeta.modality
  val cellMeta =
    MetaData.PairedTagBarcodeMeta
      .readPairedTagCellMetaFile()
      .filter(x => x.annotQuality == "Good")
  val sc2numofCell =
    cellMeta
      .groupBy(x => x.annot.sc.get)
      .map((sc, cells) => {
        val t = mods
          .map(h => {
            val rep2num = cells
              .filter(x => (h == "all") || (x.exp.modularity == h))
              .groupMapReduce(_.sample.rep)(x => 1)((x, y) => x + y)
            val r = genNumOfCell(rep2num, name = h)
            (h, r)
          })
          .toMap
        (sc, t)
      })
  val ordscs =
    sc2numofCell.keys.toVector.sortBy(x => x.split(" ").head.toInt)
  ordscs.map(sc => {
    val r   = sc2numofCell(sc)
    val str = mods.map(m => (r(m).mkString(sep2))).mkString(sep1)
    (sc, str)
  })
}

def readSubclassMetaFromChromHMM: List[PairedTagClusterChromHMM] = {
  val fnm =
    s"${TSCCMeta.projd}/meta/PairedTagSubclassMetaFromChromHMM.csv"
  readTable(fnm, sep = ",", head = true)
    .map(x =>
      new PairedTagClusterChromHMM(
          name = ClusterName(x(0), x(1), x(2), x(3)),
          hasH3K27ac = toHasHistoneModTag(x(4)),
          hasH3K27me3 = toHasHistoneModTag(x(5)),
          hasH3K4me1 = toHasHistoneModTag(x(6)),
          hasH3K9me3 = toHasHistoneModTag(x(7)),
          hasATAC = x(8).toInt > 0,
          hasDNAMeth = x(9).toInt > 0,
          hasHiC = x(10).toInt > 0
      ))
} // end of readSubclassMeta

// this is used to add more information about subclass meta.
def addNumOfCellsToSubclassMetaFromChromHMM(outf: String,
  sep: String = ","): Unit = {
  val scMetaFromChromHMM =
    readSubclassMetaFromChromHMM.map(x => (x.name.AllenIdName, x)).toMap
  val scMetaFromCells =
    initSubclassMetaFromCellMeta(sep1 = ",", sep2 = "-")
  val scMeta = scMetaFromCells.map((sc, numofCells) => {
    if (scMetaFromChromHMM.contains(sc)) {
      val r = scMetaFromChromHMM(sc)
      (r.name, numofCells, r.hasATAC, r.hasDNAMeth, r.hasHiC)
    } else {
      (new ClusterName(
              subclassNamePairedTag(sc),
              subclassNameCEMBATAC(sc),
              subclassNameDNAMeth(sc),
              sc
          ), numofCells, false, false, false)
    }
  })
  val scMetaStr = scMeta.map(x =>
    List(x._1.mkString(sep), x._2, bool2str(x._3), bool2str(x._4),
        bool2str(x._5))
      .mkString(sep))
  writeStrings2File(
      content = scMetaStr.toList,
      to = outf,
      overwrite = true,
      head = getClusterNameHeader(sep) + sep +
        ("all" +: TSCCMeta.modality).mkString(sep) + sep +
        List("ATAC", "DNAMeth", "HiC").mkString(sep)
  )
}

/**
 * Read Paired-Tag SubclassMeta file.
 *
 * This version is to read the one we generated by merging both
 * Single-Cell Meta and Subclass meta from ChromHMM meta data. This
 * function may change later if subclass meta changes.
 *
 * @return
 */
def readSubclassMeta(f: String = TSCCMeta.ptscMetafnm): Vector[PairedTagCluster] = {
  os.read.lines.stream(os.Path(f))
    .drop(1)
    .map(x => x.strip().split(","))
    .map(x => {
      new PairedTagCluster(name = new ClusterName(
              PairedTagName = x(0),
              ATACName = x(1),
              DNAmethName = x(2),
              AllenIdName = x(3)
          ), all = str2NumOfCell(x(4), sep = "-"),
          h3K4me1 = str2NumOfCell(x(5), sep = "-"),
          h3K9me3 = str2NumOfCell(x(6), sep = "-"),
          h3K27ac = str2NumOfCell(x(7), sep = "-"),
          h3K27me3 = str2NumOfCell(x(8), sep = "-"),
          atac = x(9).toInt > 0, meDNA = x(10).toInt > 0,
          hiC = x(11).toInt > 0)
    })
    .toVector
}


// FIXME: typo
val representPaiedTagSubclass: Vector[String] = Vector(
    "004_L6_IT_CTX_Glut",
    "008_L2_3_IT_ENT_Glut",
    "016_CA1_ProS_Glut",
    "017_CA3_Glut",
    "037_DG_Glut",
    "046_Vip_Gaba",
    "049_Lamp5_Gaba",
    "052_Pvalb_Gaba",
    "053_Sst_Gaba",
    "061_STR_D1_Gaba",
    "062_STR_D2_Gaba",
    "129_VMH_Nr5a1_Glut",
    "156_MB_ant_ve_Dmrta2_Glut",
    "319_Astro_TE_NN",
    "326_OPC_NN",
    "327_Oligo_NN",
    "334_Microglia_NN"
)

@main def run(): Unit = {
  // initSubclassMetaFromChromHMMCellMarkFileTable()
  addNumOfCellsToSubclassMetaFromChromHMM(outf = TSCCMeta.ptscMetafnm,
      sep = ",")
}
