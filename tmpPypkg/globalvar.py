import os
from typing import Dict, List
import re
import numpy as np
import pandas as pd

# * global variables
pt_projd = "/projects/ps-renlab2/szu/projects/amb_pairedtag"
atac_projd = "/projects/ps-renlab2/szu/projects/CEMBA2"

ptcellmetafnm = os.path.join(pt_projd, "meta",
                             "pairedtag.cell.meta.all.240626.csv")
ptscmetafnm = os.path.join(pt_projd, "meta",
                           "PairedTagSubclassMetaFromCellMetaChromHMM.csv")

ataccellmetafnm = os.path.join(
    pt_projd, "data", "snATAC", "snATAC_cellmeta.v9.8.240830.tsv")

atacsa2_anndataset_fnm = os.path.join(
    pt_projd, "data", "snATAC", "all.sa2v26.adset.h5ad")


pt2dissect_fnm = os.path.join(pt_projd, "meta",
                              "pairedtag_region2dissect.csv")

cellOracled = os.path.join(pt_projd, "data", "cellOracle")
pdcd = os.path.join(cellOracled, "cicero_pdc_bedpe")
tfscand = os.path.join(cellOracled, "tfscan")
GRNd = os.path.join(cellOracled, "GRN_subclass")
ppdcWithConnfnm = os.path.join(
    pt_projd, "data", "snATAC",
    "mba.whole.sa2subclass.pos.pdc.uniq.coacc.bedpe")
ppdcWithCorfnm = os.path.join(
    pt_projd, "data", "snATAC",
    "mba.whole.sa2subclass.pearson.pos.pdc.alignv1.tsv")

ATACPM_pbysc_fnm = os.path.join(pt_projd, "data", "snATAC",
                                "cpm_peakBysubclass.csv")
PairedTagRNACPM_scbyg_fnm = os.path.join(
    pt_projd, "data", "pairedtag_ann", "pt.RNAseq.scbyg.CPM.csv")
AllenRNAlogCOM_scbyg_fnm = os.path.join(
    pt_projd, "data", "allen", "sa2.allen.avg.logCPM.scbyg.csv")

# * functions

def load_ptscmeta() -> pd.DataFrame:
    """
    Load PairedTag Subclass-level meta data.
    - Add CEMBATACname column for ATAC-seq data usage.
    - Filter subclasses without ATAC-seq data.
    - Remove Hypendymal_NN since no CellOracle result for it.
    """
    r = pd.read_csv(ptscmetafnm, sep=",", header=0)
    # filter subclass without ATAC signals
    r = r[r.ATAC > 0]
    r['CEMBATACname'] = r.ATACName.apply(
        lambda x: re.sub("^\d+_", "", x))
    # filter the subclass without cell oracle result.
    r = r[~r.CEMBATACname.isin(["Hypendymal_NN"])]
    return r

def get_rawCEMBAdissect_in_ambPairedTag() -> Dict[str, List[str]]:
    with open(pt2dissect_fnm, "r") as f:
        lines: List[List[str]] = [line.strip().split(",") for line in f.readlines()]
    rr = {line[0]: line[1] for line in lines}
    r: Dict[str, List[str]] = {k: v.split(";") for k, v in rr.items()}
    return r


def get_snATACdissect_in_ambPairedTag(add_4D: bool = False) -> Dict[str, List[str]]:
    """
    Get CEMBA dissections for the amb PairedTag project.

    Option: Add '4D' if '4E' is included.
    See the discussion here:
    https://github.com/beyondpie/CEMBA_wmb_snATAC/discussions/20
    """
    raw = get_rawCEMBAdissect_in_ambPairedTag()
    r = {k: v if ["4E"] not in v else v + ["4D"] for k, v in raw.items()}
    return r

def getPairedTagCellMeta() -> pd.DataFrame:
    r = pd.read_csv(ptcellmetafnm, sep=",", header=0)
    r.set_index("barcode", drop=False, inplace=True)
    # filter LQ cells
    r = r[r["annotQuality"] == "Good"]
    return r

def transform_allenlabel(a: str) -> str:
    """Transform Allen official names to an unified string.
    """
    return re.sub(r" +|/|-", "_", a)

def reorder_scs(scs: List[str], sep: str = "_") -> List[str]:
    scs_index = [int(i.split(sep=sep)[0]) for i in scs]
    myord = [i[0] for i in sorted(enumerate(scs_index), key=lambda x: x[1])]
    return [scs[i] for i in myord]

def load_AllenRNA_logCPM_scbyg() -> pd.DataFrame:
    r = pd.read_csv(AllenRNAlogCOM_scbyg_fnm, sep=",",
                header=0, index_col=0)
    return r


def load_ATA_CPM_pbysc() -> pd.DataFrame:
    r = pd.read_csv(ATACPM_pbysc_fnm, sep=",", header=0, index_col=0)
    return r

def load_PairedTag_RNA_logCPM_scbyg() -> pd.DataFrame:
    r = pd.read_csv(PairedTagRNACPM_scbyg_fnm, sep=",", header=0,
                    index_col=0)
    r = np.log1p(r)
    return r
