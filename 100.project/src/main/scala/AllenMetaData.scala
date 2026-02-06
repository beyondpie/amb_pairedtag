/*
 Allen whole mouse brain Meta Data
 */
package AllenMetaData

import os._
import MetaData.TSCCMeta.AllenClusterMetaFile
import MetaData.toStringOfAllenIdLabel
import SZUtils.readTable

case class AllenCluster(
  id: Int, label: String, id_label: String
)

case class AllenSex(
  female_ratio: Float, male_ratio: Float, sex_bias: String
)

case class AllenNeuronType(
  nt_type_label: String, nt_type_combo_label: String
)

case class AllenDarkLight(
  dark_ratio: Float, light_ratio: Float
)

case class AllenSequencingSize(
  v3_size: Int, v2_size: Int, multiome_size: Int
)

case class AllenGeneMarker(
  cluster_markers: List[String], merfish_markers: List[String]
)

case class AllenAnatomy(
  annot: String, note: String, ccf_broad_freq: String,
  ccf_acronym_freq: String
)

case class AllenCellType(
  cl: Int, cluster: AllenCluster, supertype: AllenCluster,
  subclass: AllenCluster, allenClass: AllenCluster,
  anatomy: AllenAnatomy, size: AllenSequencingSize,
  marker: AllenGeneMarker, sex: AllenSex, darklight: AllenDarkLight,
  ntype: AllenNeuronType
)

val sc2spNN: Map[String, List[String]] =
  List(
    ("318_Astro_NT_NN", List("1160_Astro_NT_NN_2")),
    (
      "319_Astro_TE_NN",
      List(
        "1161_Astro_TE_NN_1",
        "1162_Astro_TE_NN_2",
        "1163_Astro_TE_NN_3",
        "1165_Astro_TE_NN_5"
      )
    ),
    ("322_Tanycyte_NN", List("1173_Tanycyte_NN_2")),
    ("323_Ependymal_NN", List("1175_Ependymal_NN_1")),
    ("325_CHOR_NN", List("1178_CHOR_NN_1")),
    ("326_OPC_NN", List("1179_OPC_NN_1", "1180_OPC_NN_2")),
    (
      "327_Oligo_NN",
      List(
        "1181_COP_NN_1",
        "1182_NFOL_NN_2",
        "1183_MFOL_NN_3",
        "1184_MOL_NN_4"
      )
    ),
    ("329_ABC_NN", List("1186_ABC_NN_1")),
    (
      "330_VLMC_NN",
      List("1187_VLMC_NN_1", "1188_VLMC_NN_2", "1189_VLMC_NN_3")
    ),
    ("331_Peri_NN", List("1191_Peri_NN_1")),
    ("332_SMC_NN", List("1192_SMC_NN_1")),
    ("333_Endo_NN", List("1193_Endo_NN_1")),
    ("334_Microglia_NN", List("1194_Microglia_NN_1")),
    ("335_BAM_NN", List("1195_BAM_NN_1"))
  ).toMap


object Allen {
  def AllenClMeta: List[AllenCellType] = {
    val rawTables =
      readTable(fnm = AllenClusterMetaFile,
        sep = "\t", head = true, skipEmptyField = true)

    rawTables.map(x =>
      AllenCellType(
          // the key commonly used in Allen
          cl = x(0).toInt,
          // Allen has no labels for cluster
          cluster = AllenCluster(
              id = x(1).toInt,
              label = "",
              id_label = x(2)
          ),
          supertype = AllenCluster(
              id = x(3).toInt,
              label = x(4),
              id_label = x(5)
          ),
          subclass = AllenCluster(
              id = x(6).toInt,
              label = x(7),
              id_label = x(8)
          ),
          // Allen use different order for class
          allenClass = AllenCluster(
              id = x(10).toInt,
              label = x(9),
              id_label = x(11)
          ),
          anatomy = AllenAnatomy(
              annot = x(12),
              note = x(13),
              ccf_broad_freq = x(14),
              ccf_acronym_freq = x(15)
          ),
          size = AllenSequencingSize(
              v3_size = x(16).toInt,
              v2_size = x(17).toInt,
              multiome_size = x(18).toInt
          ),
          marker = AllenGeneMarker(
              cluster_markers = x(19).split(",").toList,
              merfish_markers = x(20).split(",").toList
          ),
          sex = AllenSex(
              female_ratio = x(21).toFloat,
              male_ratio = x(22).toFloat,
              sex_bias = x(23)
          ),
          darklight = AllenDarkLight(
              dark_ratio = x(24).toFloat,
              light_ratio = x(25).toFloat
          ),
          ntype = AllenNeuronType(
              nt_type_label = x(26),
              nt_type_combo_label = x(27)
          )
      ))
  } // end of AllenClMeta

  /*
   Transform raw subclass label or subclass id label
     to subclass name recorded in CEMBA ATAC.
   */
  def subclassNameCEMBATAC(sc: String): String = {
    sc.replaceAll("^\\d+_", "")
      .replace("  ", " ")
      .replace(" ", "_")
      .replace("/", "-")
  }

  /*
   Transform raw subclass label or subclass id label
     to subclass name recorded in PairedTag.
   */
  def subclassNamePairedTag(sc: String): String = {
    sc.replace("  ", " ")
      .replace(" ", "_")
      .replace("/", "_")
      .replace("-", "_")
  }

  /*
   Transform raw subclass label or subclass id label
   to subclass name recorded in CEMBA DNA Methylation.
   In CEMBA DNA Methylation subclass name,
   All the spaces, "/" are changed into "_", but keep "-" as "-".
  */
  def subclassNameDNAMeth(sc: String): String = {
    sc.replace("  ", " ")
      .replace(" ", "_")
      .replace("/", "_")
  }

  def subclassIdNamePairedTag2raw: Map[String, String] = {
    AllenClMeta
      .map(x => x.subclass.id_label)
      .map(x => (subclassNamePairedTag(x), x))
      .toMap
  }

  def subclassIdNamePairedTag2ATAC: Map[String, String] = {
    AllenClMeta
      .map(x => x.subclass.id_label)
      .map(x => (subclassNamePairedTag(x), subclassNameCEMBATAC(x)))
      .toMap
  }

  def subclassIdNamePairedTag2DNAMeth: Map[String, String] = {
    AllenClMeta
      .map(x => x.subclass.id_label)
      .map(x => (subclassNamePairedTag(x), subclassNameDNAMeth(x)))
      .toMap
  }

} // end of Allen
