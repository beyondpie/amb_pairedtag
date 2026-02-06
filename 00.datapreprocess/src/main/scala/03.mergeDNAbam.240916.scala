import os._
import scala.collection.parallel.CollectionConverters.*
import MetaData.TSCCMeta.{ptscDNAbamd, ptclusterDNAbamd, ptCellMetaFile}
import MetaData.PairedTagBarcodeMeta
import MetaData.toStringOfAllenIdLabel
import MetaData.TSCCTools.samtools
import MetaData.TSCCTools.MergeBamUsingSamtools
import scala.collection.immutable.HashMap
import scala.util.Random
import SZUtils.readLines
import SZUtils.writeStrings2File
import MetaData.TSCCMeta.getBarcode2Bam
import MetaData.PairedTagBarcodeMeta.readPairedTagCellMetaFile

object MergeDNAbam240916 {
  val projd = "/projects/ps-renlab2/szu/projects/amb_pairedtag"
  // in encoder.
  val samtools = "/home/szu/mambaforge/bin/samtools"
  val barcodeListDir =
    os.Path("/projects/ps-renlab2/zhw063/99.MouseBrainPairedTag") /
      "20240916_analysis" /
      "cell_barcodes_output"

  val todir = s"${projd}/data/mergeDNABam240916"

  val cellMeta = readPairedTagCellMetaFile()
  val barcode2bam: Map[String, String] = getBarcode2Bam(cellMeta)

  val listOfBarcodesfiles: Map[String, String] =
    os.list(barcodeListDir)
      .filter(x => x.toString.contains("txt"))
      .map(fullpath => (fullpath.baseName, fullpath.toString))
      .toMap

  val listOfBamsfiles: Map[String, String] =
    listOfBarcodesfiles.keys
    .map(k => (k, s"${todir}/log/${k}.bamfnms.list"))
    .toMap

  listOfBarcodesfiles.foreach { (k, f) =>
    {
      // load from f
      val barcodes: List[String] =
        readLines(listOfBarcodesfiles(k), head = false)
      // map to the corresponding bam fnms
      val barcodefnms: List[String] = barcodes.map(x => barcode2bam(x))
      // write barcodefnms to listOfBamsfiles(k)
      writeStrings2File(
        content = barcodefnms,
        to = listOfBamsfiles(k),
        overwrite = true
      )
    }
  }

  val flagfnms: Map[String, String] =
    listOfBarcodesfiles.keys
    .map(k => (k, s"${todir}/flag/${k}.done"))
    .toMap
  val logfnms: Map[String, String] =
    listOfBarcodesfiles.keys
    .map(k => (k, s"${todir}/log/${k}.log"))
    .toMap
  val outMergedBamfnms: Map[String, String] =
    listOfBarcodesfiles.keys
    .map(k => (k, s"${todir}/out/${k}.bam"))
    .toMap

  def main(args: Array[String]) = {
    val mergeBamTasks: List[MergeBamUsingSamtools] =
      listOfBamsfiles
      .map((k, f) => {
        new MergeBamUsingSamtools(
          samtools = samtools,
          flagfnm = flagfnms(k),
          logfnm = logfnms(k),
          fnmofListOfBams = listOfBamsfiles(k),
          toBamfile = outMergedBamfnms(k),
          filterMAPQ = false,
          skip = true
        )
      })
        .toList

    mergeBamTasks.par.foreach { x =>
      x.run()
    }

  } // end of main
} // end of object containing main
