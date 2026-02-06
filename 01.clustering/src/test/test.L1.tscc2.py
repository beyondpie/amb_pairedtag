import os
import snapatac2 as sa2

a = sa2.read(
    os.path.join("/tscc/projects/ps-renlab2/szu/projects/amb_pairedtag/",
                 "01.clustering/out/tscc2/L1/anns",
                 "RNA_k8_npc50_k40_L0_c0.h5ad")
)
