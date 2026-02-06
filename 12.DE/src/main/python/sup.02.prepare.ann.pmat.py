import os
import snapatac2 as sa2

# * meta
projd = "/projects/ps-renlab2/szu/projects/amb_pairedtag"
sa2h5adir = os.path.join(
    "/projects/ps-renlab2/zhw063", "99.MouseBrainPairedTag", "20241202_for_Songpeng"
)
outd = os.path.join(projd, "12.DE", "out")


# * prepare H3K27ac pmat
raw_ann_h27ac = sa2.read(
    filename=os.path.join(sa2h5adir, "2024701.brain.snapatac2.h3k27ac.h5ad"),
    backed=None,
)
# repeat exp columns lead the copy failed.
raw_ann_h27ac.obs.drop(columns=["exp"], inplace=True)
ann_h27ac = raw_ann_h27ac[
    raw_ann_h27ac.obs["annotQuality"].astype(str) == "Good"
].copy()
ann_h27ac = sa2.pp.make_peak_matrix(
    ann_h27ac, peak_file=os.path.join(outd, "H3K27ac.all.bestspm.peak.bed")
)
ann_h27ac.write_h5ad(filename=os.path.join(outd, "ann.H3K27ac.pmat.h5ad"))

# * prepare H3K9me3 pmat
raw_ann_k9 = sa2.read(
    filename=os.path.join(sa2h5adir,
                          "2024701.brain.snapatac2.h3k9me3.50kbin.h5ad"),
    backed=None,
)

raw_ann_k9.obs.drop(columns=["exp"], inplace=True)
ann_k9 = raw_ann_k9[raw_ann_k9.obs["annotQuality"].astype(str) == "Good"].copy()
ann_k9 = sa2.pp.make_peak_matrix(
    ann_k9, peak_file=os.path.join(sa2h5adir, "H3K9me3.merged.all.blv2.me.peak.bed")
)
ann_k9.write_h5ad(filename=os.path.join(outd, "ann.H3K9me3.pmat.h5ad"))

# * check ATAC-seq snapATAC2 object
atac_sa2_fnm2 = os.path.join(
    projd, "data", "snATAC", "all.sa2v26.pmat.snATAC.nofrag.h5ad"
)
r2 = sa2.read(filename=atac_sa2_fnm2, backed=None)

# polars data.frame
sc2cnt = r2.obs["subclass"].value_counts()
sc2cnt.write_csv(
    file=os.path.join(projd, "meta", "snATAC.subclass2cnt.csv"),
    include_header=True,
    separator=","
)
 
