import os
import pandas as pd
from multiprocessing import Pool

import tmpPypkg.globalvar as gv
from tmpPypkg.utils import reimport

reimport("tmpPypkg")
reimport("tmpPypkg.globalvar")

bedzstd = os.path.join(
    gv.pt_projd, "data/snATAC", "subclass_bed_zst")

outd = os.path.join(
    gv.pt_projd, "data/snATAC", "subclassbed4ChromHMM")

# * load subclasses in pariedtag
ptcm = gv.getPairedTagCellMeta()
ptsc = ptcm["annot.sc"].unique()
# 165
ptsc = [gv.transform_allenlabel(i) for i in ptsc]

# * get list of subclass in snATAC subclass_bed
atacscs = [f.split(".")[0]
           for f in os.listdir(bedzstd)
           if os.path.isfile(os.path.join(bedzstd, f))]

# 152
jscs = [sc for sc in atacscs if sc in ptsc]
# * for each of the intersected subclass, perform the update.

def update_bed(sc) -> None:
    f = os.path.join(bedzstd, f"{sc}.bed.zst")
    a = pd.read_csv(f, sep="\t", header=None)
    a.insert(5, column='strand', value="+")
    a.to_csv(os.path.join(outd, f"{sc}.chromHMM.bed"),
            sep="\t", header=False, index=False)


# This needs about 2 hours in single-thread mode
# for sc in jscs:
#     print(f"Running transformation of bed for: {sc}")
#     update_bed(sc)

if __name__ == '__main__':
    with Pool(40) as p:
        p.map(update_bed, jscs)
