package MetaData

object PairedTagQCMeta {
  val tfoutd = TSCCMeta.tfoutd
  val unCoEmbedLQFiles: List[String] =
    List(
      "AMY.to.be.deleted.L5r.txt",
      "CPU.to.be.deleted.3.txt",
      "HIP.to.be.deleted.1.txt",
      "HIP.to.be.deleted.2.txt",
      "HIP.to.be.deleted.3.txt",
      "HYP.to.be.deleted.txt",
      "NAC.to.be.deleted.txt",
      "PFC.to.be.deleted.2.txt",
      "VTA.to.be.deleted.1.txt",
      "VTA.to.be.deleted.2.txt"
    ).map(x => s"${tfoutd}/un_coembed/$x")

  def unCoEmbedHQFiles: List[String] =
    List(
      "CPU.enrich.L5r.with.annot.txt",
      "ERC.enrich.L5r.with.annot.txt"
    ).map(x => s"${tfoutd}/un_coembed/$x")
}
