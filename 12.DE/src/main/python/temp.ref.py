import numpy as np
import snapatac2 as snap


brain_h3k4me1 = snap.read("./2024701.brain.snapatac2.h3k4me1.h5ad", backed=None)
brain_h3k4me1.obs.drop(columns=["cellular_barcode"], inplace=True)
brain_h3k4me1.obs.drop(columns=["exp"], inplace=True)
brain_h3k4me1_clean = brain_h3k4me1[brain_h3k4me1.obs["annotQuality"].astype(str) == "Good"].copy()

h3k4me1_peak_mat = snap.pp.make_peak_matrix(
    brain_h3k27me3_clean,
    peak_file="./H3K4me1.merged.all.blv2.me.peak.bed"
)
subclass_marker_peaks_h3k4me1 = snap.tl.marker_regions(
    h3k4me1_peak_mat, groupby="annot.sc", pvalue=0.01
)

# Initialize a dictionary to store results for each cell group
results_sc = {}

# Loop through each unique cell group
for cellgroup in np.unique(brain_h3k4me1_clean.obs["annot.sc"]):

    # Define the two groups: the current cell group and the rest
    group1 = brain_h3k4me1_clean.obs["annot.sc"] == cellgroup
    group2 = brain_h3k4me1_clean.obs["annot.sc"] != cellgroup

    # Get the indices of cells in group1 and group2
    indices_group1 = np.where(group1)[0]
    indices_group2 = np.where(group2)[0]

    # Downsample group2 to match the size of group1
    if len(indices_group1) < len(indices_group2):
        # Randomly select cells from group2 to match the number of cells in group1
        downsampled_indices_group2 = np.random.choice(
            indices_group2, size=len(indices_group1) * 2, replace=False
        )
    else:
        # If group2 has fewer or equal cells than group1, use all cells in group2
        downsampled_indices_group2 = indices_group2

    # Create downsampled boolean masks for group2
    downsampled_group2 = np.isin(np.arange(len(group2)),
                                 downsampled_indices_group2)

    # Perform differential peak analysis with downsampled group2
    diff_peaks = snap.tl.diff_test(
        h3k4me1_peak_mat,
        cell_group1=group1,
        cell_group2=downsampled_group2,
        direction="positive",
        min_pct=0.005,
    )

    # Save the results in the dictionary with the cell group as the key
    results_sc[cellgroup] = diff_peaks
