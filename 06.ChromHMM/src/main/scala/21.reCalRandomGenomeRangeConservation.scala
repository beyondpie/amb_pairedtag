import os._
import MetaData.GenomeMeta.{genNullGenomicRanges, filterMouseSignalRanges}
import SZUtils.{writeStrings2File, str2path, path2str}
import ChromHMMSequenceConservation.deepToolComputeMatrix
import bioscala.LightCoord.GenomeCoord.{GenomeCoord, GenomeCoords}
import bioscala.LightCoord.GenomeCoord.fromBed
import bioscala.LightCoord.GenomeCoord.mkStringGenomeCoord
 
object CalRandomGenomeRangeConservation {
  val projd = os.Path("/projects/ps-renlab2/szu/projects/amb_pairedtag")
  val workd = projd / "06.ChromHMM"
  val outd = workd / "out"
  val shuffleBedProg = "/home/szu/miniforge3/bin/shuffleBed"
  val computeMatrixProg = "/home/szu/miniforge3/bin/computeMatrix"
  val phastConsfnm = MetaData.GenomeMeta.phastConsfnm
  val rdmfnm = outd / "mm10.null2.background.bed"
  val blfnm = MetaData.GenomeMeta.mm10BlacklistBed

  val rawRdn: GenomeCoords = fromBed(outd / "mm10.null.background.bed")
  val fltRdn = filterMouseSignalRanges(rawRdn)

  writeStrings2File(
    content = fltRdn.map(x => mkStringGenomeCoord(x, "\t", false)),
    to = rdmfnm,
    head = ""
  )

  deepToolComputeMatrix(
    bin = computeMatrixProg,
    regionfnm = rdmfnm,
    bwfnm = phastConsfnm,
    gzmatfnm = s"$outd/mm10.random2.bg.computeMatrix.mat.gz",
    blacklistfnm = blfnm,
    ncore = 10
  )

}



