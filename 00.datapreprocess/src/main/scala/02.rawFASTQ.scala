import os._
import monocle.syntax.all._

import SZUtils.readTable
import ZSTCompress.gz2zst
import MetaData.PairedTagSublib.{SublibRawFASTQ, readSublib}
import MetaData.PairedTagSublib.Sublib
import SZUtils.writeStrings2File

object OrgRawFASTQs {
  val old = "/projects/ps-renlab"
  val d = "/projects/ps-renlab2"

  val projd    = s"${d}/szu/projects/amb_pairedtag"
  val sublibf  = s"${projd}/meta/amb_pairedtag_sublibrary.csv"
  val fastqd   = s"${d}/szu/shared_data/wmb_EP/fastqRgz"
  val demultid = s"${d}/szu/shared_data/wmb_EP/fromNASguac"

  def getFileLink(fnm: String): String = os.Path(fnm) match {
    case p if os.isLink(p) =>
      os.readLink.absolute(p).toString match {
        case r1 if os.exists(os.Path(r1)) => r1
        case r1 if r1.contains(s"${old}/zhw063") =>
          r1.replace(s"${old}/zhw063", s"${d}/ps-renlab/zhw063")
        // case r1 if r1.contains(s"${d}/zhw063/demultiplex") =>
        //   r1.replace(s"${d}/zhw063/demultiplex",
        //       s"${d}/zhw063/demultiplex/backup")
        case r1
            if r1.contains(s"${d}/zhw063/demultiplex") || r1.contains(
                "/home/zhw063/renlab2/demultiplex") => {
          val localfnm = os.Path(r1).baseName + ".gz"
          s"${demultid}/${localfnm}"
        }
        case r1
            if r1.contains(
                s"${d}/zhw063/15.MouseBrainPTexp8B_20230508") =>
          r1.replace(s"${d}/zhw063/15.MouseBrainPTexp8B_20230508",
              s"${d}/zhw063/99.MouseBrainPairedTag/08.MouseBrainExp8" +
                "/03.MouseBrainPTexp8B_20230508")
        case r1
            if r1.contains(
                s"${d}/zhw063/17.MouseBrainPTexp10_230521") =>
          r1.replace(s"${d}/zhw063/17.MouseBrainPTexp10_230521",
              s"${d}/zhw063/99.MouseBrainPairedTag/10.MouseBrainExp10" +
                "/04.MouseBrainPTexp10_230521")
        case r1 =>
          r1.replace(s"${d}/zhw063/99.MouseBrainPairedTag", fastqd)
      }
    case p if p.toString.contains(s"${old}/zhw063") =>
      p.toString.replace(old, d)
    case _ => fnm
  }

  def modifyFile(fnm: String): String = ???

  def main(args: Array[String]) = {
    // * add all the directories containing raw fastq files
    //   data copy from szu NAS.
    // RNA_name / DNA_name 2 fastq.gz map.
    // eg, ZW331 to /somewhere/ZW331_middle_R1(and R2)_suffix.fastq.gz
    // val rawFASTQgzd: Map[String, (String, String)] = ???

    // * check and update, for each sublib, if R1 and R2 exists in above files.
    val origSublibs: List[Sublib] =
      readSublib(sublibf)
        .map(x =>
          x.focus(_.workd)
            .modify(y =>
              y.replace(s"${old}/zhw063/renlab2", s"${d}/zhw063")))

    val origSublibFastqs: List[SublibRawFASTQ] =
      origSublibs
        .map(x =>
          List(
              new SublibRawFASTQ(
                  sublibid = x.sublibid,
                  tag = x.idDNA,
                  R1 = s"${x.workd}/01.rawdata/${x.idDNA}_R1.fq.gz",
                  R2 = s"${x.workd}/01.rawdata/${x.idDNA}_R2.fq.gz"
              ),
              new SublibRawFASTQ(
                  sublibid = x.sublibid,
                  tag = x.idRNA,
                  R1 = s"${x.workd}/01.rawdata/${x.idRNA}_R1.fq.gz",
                  R2 = s"${x.workd}/01.rawdata/${x.idRNA}_R2.fq.gz"
              )
          ))
        .flatten

    val sublibFastqs = origSublibFastqs.map(x =>
      x.focus(_.R1)
        .modify(f => getFileLink(f))
        .focus(_.R2)
        .modify(f => getFileLink(f)))

    // * check files all exist
    val notExistSublibR1Fastqs =
      sublibFastqs.filter(x => !os.exists(os.Path(x.R1)))

    val notExistSublibR2Fastqs =
      sublibFastqs.filter(x => !os.exists(os.Path(x.R2)))

    // output
    // use /tscc/projects for general access
    val r1 = sublibFastqs.map(x => List(
      x.tag, x.R1.replace("/projects/", "/tscc/projects/")))
    val r2 = sublibFastqs.map(x => List(x.tag,
      x.R2.replace("/projects/", "/tscc/projects/")))

    writeStrings2File(
      content = r1.map(x => x.mkString(",")),
      to = s"${projd}/meta/pairedtag.Raw.FASTQgz.R1.csv",
      head = ""
    )

    writeStrings2File(
      content = r2.map(x => x.mkString(",")),
      to = s"${projd}/meta/pairedtag.Raw.FASTQgz.R2.csv",
      head = ""
    )

  } // end of main
}
