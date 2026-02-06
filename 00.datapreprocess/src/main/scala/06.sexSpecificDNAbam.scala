import os._
import MetaData.TSCCMeta
import MetaData.TSCCTools.MergeBamUsingSamtools
import scala.collection.parallel._
import scala.collection.parallel.CollectionConverters.*
import MetaData.PairedTagBarcodeMeta.readPairedTagCellMetaFile
import MetaData.TSCCMeta.getBarcode2Bam
import MetaData.readSubclassMeta
import SZUtils.readLines
import SZUtils.writeStrings2File
import AllenMetaData.Allen.subclassNamePairedTag
// non-neuronal cells
// four histones
// split by sex
@main def getSexSpecificDNAbam: Unit = {
  val projd         = TSCCMeta.projd
  val workd         = projd + "/00.datapreprocess"
  val flagd: String = workd + "/flag"
  val logd: String  = workd + "/log"
  val outMergeBamd  = projd + "/data/ptDNAbam/sexbam"
  val samtools      = "/home/szu/mambaforge/bin/samtools"
  val sexs          = Vector("Male", "Female")
  val hs = Vector("H3K27me3", "H3K27ac", "H3K4me1", "H3K9me3")
  // load scMeta to get barcodes
  val cellMeta =
    readPairedTagCellMetaFile()
      .filter(x => x.annotQuality == "Good")

  // set up subclass and sex pairs
  val ptscMeta = readSubclassMeta().filter(_.atac)

  val scsexh: Vector[(String, String, String)] =
    ptscMeta
      .map(_.name.AllenIdName)
      .filter(x => x.endsWith("NN"))
      .flatMap(x => sexs.map(y => (x, y)))
      .flatMap(x => hs.map(h => (x._1, x._2, h)))

  // load barcode2bam
  val b2b = getBarcode2Bam(cellMeta)
  // - make sure the files are sorted and index before merging
  // generate list of bams file
  val listOfBamsfiles: Map[(String, String, String), String] =
    scsexh
      .map((sc, sex, h) => {
        ((sc, sex, h),
            logd +
              s"/${subclassNamePairedTag(sc)}.$sex.$h.bams.list")
      })
      .toMap

  listOfBamsfiles.foreach((sc_sex_h, fnm) => {
    val sc       = sc_sex_h._1
    val sex      = sc_sex_h._2
    val h        = sc_sex_h._3
    val barcodes =
      cellMeta
        .filter(x =>
          x.sample.sex == sex &&
            x.annot.sc.get == sc &&
            x.exp.modularity == h)
        .map(x => x.barcode)
    val fnms = barcodes.map(b => b2b(b))
    if (fnms.nonEmpty) {
      writeStrings2File(
          content = fnms,
          to = fnm,
          overwrite = true,
          head = ""
      )
    } else {
      println(s"$sc $sex $h has no barcodes.")
    }
  })
  // generate the tasks
  val t =
    scsexh
      .filter((sc, sex, h) => {
        os.exists(os.Path(logd +
          s"/${subclassNamePairedTag(sc)}.$sex.$h.bams.list"))
      })
      .map((sc, sex, h) => {
        val ptsc = subclassNamePairedTag(sc)
        val r    = new MergeBamUsingSamtools(
          samtools = samtools,
          flagfnm = flagd + s"/$ptsc.$sex.$h.mergeBam.done",
          logfnm = logd + s"/$ptsc.$sex.$h.mergeBam.log",
          toBamfile = outMergeBamd + s"/$ptsc.$h.$sex.bam",
          fnmofListOfBams = listOfBamsfiles((sc, sex, h)),
          filterMAPQ = false,
          skip = true
        )
        ((sc, sex, h), r)
      })
      .toMap


  // run mergeBam
  t.par.foreach((sc_sex_h, a) => {
    val sc  = sc_sex_h._1
    val sex = sc_sex_h._2
    val h   = sc_sex_h._3
    println(s"mergeBam for $sc $sex $h.")
    a.run()
    println(s"mergeBam for $sc $sex $h done.")
  })

}
