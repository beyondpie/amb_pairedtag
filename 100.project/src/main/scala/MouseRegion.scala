package MouseRegion

import scala.collection.immutable.HashMap

object RegionMapping {

  val PairedTagRegion2MajorRegion = HashMap(
    "CPU" -> "CPU",
    "HYP" -> "HYP",
    "HCa" -> "HIP",
    "HCp" -> "HIP",
    "ERC" -> "ERC",
    "AMY" -> "AMY",
    "NAC" -> "NAC",
    "VTA_SnR" -> "VTA",
    "PFC" -> "PFC"
  )

  val AllenRegion2MajorRegion = HashMap(
    "STR_STRd" -> "CPU",
    "HY_HYa1" -> "HYP",
    "HY_HYa2" -> "HYP",
    "HY_HYml" -> "HYP",
    "HY_HYm2" -> "HYP",
    "HY_HYpl" -> "HYP",
    "HY_HYpm" -> "HYP",
    "HY_LZ" -> "HYP",
    "HY_MEZ_PVZ_PVR" -> "HYP",
    "HIP_CA" -> "HIP",
    "HIP" -> "HIP",
    "ENT" -> "ERC",
    "STR_sAMY" -> "AMY",
    "STR_STRv" -> "NAC",
    "MB_VTA_SN" -> "VTA",
    "ACA" -> "PFC",
    "PL_ILA_ORB" -> "PFC"
  )
}
