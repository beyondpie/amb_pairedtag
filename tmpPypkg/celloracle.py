import re
from typing import List, Dict
import pandas as pd
import celloracle as co
from celloracle.network_analysis.links_object import Links
from celloracle.motif_analysis.tfinfo_core import load_TFinfo, TFinfo

from tmpPypkg.globalvar import ppdcWithCorfnm, ppdcWithConnfnm

def read_pdc_bedpe(pdcfnm: str) -> pd.DataFrame:
    """
    Read proximal-distal connections saved as bedpe file.
    """
    pdc: pd.DataFrame = pd.read_csv(pdcfnm, sep="\t", header=None)
    pdc.columns = ["pchr", "pstart", "pend",
                   "dchr", "dstart", "dend",
                   "name", "conn", "strand1", "strans2"]
    index = pdc.groupby("name")["conn"].idxmax()
    pdc = pdc.loc[index]
    pdc.set_index("name", drop=False, inplace=True,
                  verify_integrity=True)
    pdc.insert(loc=pdc.shape[1], column="CRE",
               value=pdc.apply(
                   lambda x: f"{x.dchr}:{x.dstart}-{x.dend}", axis=1))
    pdc.insert(loc=pdc.shape[1], column="geneGR",
               value=pdc.apply(
                   lambda x: f"{x.dchr}:{x.pstart}-{x.pend}", axis=1))
    pdc.insert(loc=pdc.shape[1], column="gene",
               value=pdc.apply(lambda x: x.name.split("|")[0], axis=1))
    return pdc


def load_ppdc_with_conn() -> pd.DataFrame:
    r = read_pdc_bedpe(ppdcWithConnfnm)
    return r

def load_ppdc_with_cor() -> pd.DataFrame:
    r: pd.DataFrame = load_ppdc_with_conn()
    x = pd.read_csv(ppdcWithCorfnm, sep="\t", header=0)
    x.set_index("conns", drop=False, inplace=True,
                verify_integrity=True)
    r.insert(column="cor", loc=r.shape[1],
             value=x.loc[r.index, "pcc"])
    r.insert(column="corFDR", loc=r.shape[1],
             value=x.loc[r.index, "FDR"])
    return r

def load_ppdcSum() -> pd.DataFrame:
    r = load_ppdc_with_cor()
    x = pd.DataFrame(
        {
            "CRE": r.CRE,
            "geneGR": r.geneGR,
            "gene": r.gene,
            "cor": r.cor,
            "conn": r.conn
        }
    )
    x.reset_index(inplace=True, drop=True)
    name = x.apply(lambda x: f"{x.gene}|{x.CRE}", axis=1)
    x.index = name
    return x

def getTF2gene(linkfnm: str, ptsc: str, allensc: str,
               atacsc: str,
               allenRNAlogCPMscbyg: pd.DataFrame,
               ptRNAlogCPMscbyg: pd.DataFrame) -> pd.DataFrame:
    """
    From the CellOracle's GRN, load the final TFs to target genes.

    In whole mouse brain snATAC-seq data analaysis, we use default
    parameters for filtering the edges:
    p=0.001, weight='coef_abs', threshold_number=10000.

    - In the field of links_dict, it records all the edges from GRN.
    - In the field of filtered_links, it only records the edges after
      filtering.

    NOTE: in our previous analys, for each subclass, we run GRN model
    seperatedly, which might not be the best strategy for the model if it based
    on Bayesian to borrow information from other subclasses.

    Args:
      linkfnm: h5ad filename of links from CellOracle
      atacsc: subclass name for this analysis, use `None` to extract
          the only one value under filtered_links field

    Return:
      pd.DataFrame: TF to target genes with columns including
      - TF gene symbols
      - target gene symbols
      - Allen's RNAseq expression for TFs
      - Allen's RNAseq expression for target genes
      - PairedTag RNA logCPM for TFs
      - PairedTag RNA logCPM for tareget genes
      - and coef_mean and neglogp
    """

    grn: Links = co.utility.load_hdf5(linkfnm)
    if atacsc is None:
        _, r = grn.filtered_links.popitem()
    else:
        r = grn.filtered_links[atacsc]
    genes = allenRNAlogCPMscbyg.columns.intersection(
        ptRNAlogCPMscbyg.columns
    )
    index = r.source.isin(genes) & r.target.isin(genes)
    r = r[index]
    r.insert(loc=r.shape[1],
             column="source_ptRNAlogCPM",
             value=ptRNAlogCPMscbyg.loc[ptsc, r.source].to_list())
    r.insert(loc=r.shape[1],
             column="source_allenRNAlogCPM",
             value=allenRNAlogCPMscbyg.loc[allensc, r.source].to_list())
    r.insert(loc=r.shape[1],
             column="target_ptRNAlogCPM",
             value=ptRNAlogCPMscbyg.loc[ptsc, r.target].to_list())
    r.insert(loc=r.shape[1],
             column="target_allenRNAlogCPM",
             value=allenRNAlogCPMscbyg.loc[allensc, r.target].to_list())
    return r


def getCRE2tf(tfinfofnm: str,
              atacsc: str,
              atacCPMpbysc: pd.DataFrame) -> pd.DataFrame:
    """
    From the GimmeMotif's tfinfo file, load CRE to TFs.

    Args:
      tfinfo: TFinfo object from CellOracle

    Return:
      pd.DataFrame: CRE to TFs wtith columns including
      - CRE genomic strings, like chr1:1-500
      - TF gene symbols
      - ATACSeq CPM values for CRE
      - gimmemotif scores
    """
    # For filtering, we use the one shown in CellOracle tutorial
    tfinfo: TFinfo = load_TFinfo(tfinfofnm)
    tfinfo.reset_filtering()
    # parameters we used in snATAC-seq paper
    tfinfo.filter_motifs_by_score(
        threshold=10.0,
        method="cumulative_score")
    motif2tf: Dict[str, List[str]] = {
        k: list(set(v))
        for k, v in tfinfo.dic_motif2TFs.items()
        if len(v) > 0
    }

    t: pd.DataFrame = tfinfo.scanned_filtered[
        tfinfo.scanned_filtered.apply(
            lambda x: x.motif_id in motif2tf, axis=1)]
    t.set_index(["seqname", "motif_id"],
                drop=False, inplace=True)

    scandf = tfinfo.scanned_df
    scandf.set_index(["seqname", "motif_id"],
                     drop=False,
                     inplace=True)

    # only use motif shown in motif2tf
    scandf = scandf[scandf.index.isin(t.index)].reset_index(drop=True)[
        ["seqname", "motif_id", "factors_direct",
         "factors_indirect", "score", "pos", "strand"]
    ]
    scandf.insert(loc=scandf.shape[1],
                  column="tf",
                  value=scandf.motif_id.apply(lambda x: motif2tf[x]))
    # expand List[str] field of tf, so each tf has a record.
    y = scandf.explode("tf").drop(
        ["factors_direct", "factors_indirect"], axis=1)

    # insert targeted genes from pdc
    y = y.merge(tfinfo.peak_df, how="inner", left_on="seqname",
                right_on="peak_id")

    # reformat genomic range.
    seqname = y.seqname.apply(
        lambda x: re.sub("-", ":", x.replace("_", "-"), count=1))
    y["seqname"] = seqname

    # add atac signal if not None.
    if atacCPMpbysc is not None:
        y.insert(loc=y.shape[1], column="atacCPM",
                value=atacCPMpbysc.loc[
                    seqname, atacsc].reset_index(drop=True))

    return y

def getCRE2gene(pdcfnm: str, ppdcSum: pd.DataFrame) -> pd.DataFrame:
    """
    From the cicero and correlations, load CRE to target genes.

    Args:
      pdcfnm: subclass-specific proximal-distal connection fnm.
      ppdcSum: global positivie proximal-distal connections

    Return:
      pd.DataFrame: CRE to target genes with column s including
      - CRE genomic strings
      - target genes' genomeic strings
      - target genes' gene symbols
      - connection scores
      - positive direction scores, if not in ppdc, will be 0.0
    """
    pdc: pd.DataFrame = read_pdc_bedpe(pdcfnm)
    r = pd.DataFrame({
        "CRE": pdc.CRE,
        "gene": pdc.gene,
        "conn": pdc.conn,
        "cor": 0.0
    })
    name_has_cor = pdc.index[pdc.index.isin(ppdcSum.index)]
    r.loc[name_has_cor, "cor"] = ppdcSum.loc[name_has_cor, "cor"]
    r.reset_index(drop=True, inplace=True)
    return r

def geteRegulon(CRE2gene: pd.DataFrame,
                CRE2tf: pd.DataFrame,
                tf2gene: pd.DataFrame) -> pd.DataFrame:
    """
    Establish eRegulons by manually linking CRE2gene (cicero),
    CRE2tf (gimme) and tf2gene (CellOracle). 
    """
    x = CRE2tf.merge(CRE2gene, how="inner",
                left_on=["seqname", "gene_short_name"],
                right_on=["CRE", "gene"])

    x = x.merge(tf2gene, how="inner", left_on=["tf", "gene"],
                right_on=["source", "target"])

    x = x[["seqname", "tf", "target",
           "conn", "coef_mean", "cor", "atacCPM",
           "source_ptRNAlogCPM", "target_ptRNAlogCPM",
           "source_allenRNAlogCPM", "target_allenRNAlogCPM",
           "motif_id", "score", "pos", "strand"
           ]]
    x.rename(columns={
        "seqname": "CRE",
        "motif_id": "motif",
        "target": "gene",
        "score": "gimmiescore",
        "pos": "motifpos",
        "strand": "motifstrand",
        "coef_mean": "coefCellOracle",
        "cor": "corSubclass",
        "conn": "connCicero",
        "source_ptRNAlogCPM": "tfptRNAlogCPM",
        "target_ptRNAlogCPM": "geneptRNAlogCPM",
        "source_allenRNAlogCPM": "tfAllenlogCPM",
        "target_allenRNAlogCPM": "geneAllenlogCPM"
    }, inplace=True)
    return x
