import os
import sys
from typing import List
import snapatac2 as sa2


snATAC_sample = sys.argv[1]

snATAC_cellmeta_fnm = os.path.join(
    "/projects/ps-renlab2/szu/projects/CEMBA2",
    "meta/wmb.cellmeta.v9.7.tsv")

snATAC_samples = os.path.join(
    "/projects/ps-renlab2/szu/projects/amb_pairedtag",
    "data/snATAC",
    "snapatac2.6_samples")
tmpdir = "/projects/ps-renlab2/szu/tmpdir"
snATAC_fragd = os.path.join(
        "/projects/ps-renlab2/szu/projects/CEMBA2",
        "00.data.preprocess",
        "snapatac2_pp_out/raw_fragment_files")

snATAC_sample2barcode = os.path.join(
    "/projects/ps-renlab2/szu/projects/amb_pairedtag",
    "05.CRE/src/main/resource",
    "snATAC_sample2barcode")


if __name__ == '__main__':
    def get_barcodes(s: str) -> List[str]:
        with open(os.path.join(snATAC_sample2barcode,
                               f"{s}.sample.txt"), 'r') as f:
            barcodes = [b.strip() for b in f.readlines()]
        return(barcodes)
    def get_sa2obj(wmb_sample: str) -> None:
        fragfnm = os.path.join(snATAC_fragd, f"{wmb_sample}.sa2.frag.tsv")
        barcodes = get_barcodes(wmb_sample)
        print(f"process {fragfnm} on {len(barcodes)} barcodes ...")
        outf = os.path.join(snATAC_samples, f"{wmb_sample}.sa2v26.h5ad")
        if os.path.exists(outf):
            print(f"{outf} exists, and remove it.")
            os.remove(outf)
        # fragment has been shift before.
        a = sa2.pp.import_data(
            fragment_file=os.path.join(snATAC_fragd, fragfnm),
            whitelist=barcodes,
            file=outf,
            shift_left=0,
            shift_right=0,
            chrom_sizes=sa2.genome.mm10,
            sorted_by_barcode=False,
            tempdir=tmpdir,
            chrM=["chrM", "M"],
            n_jobs=1
        )
        a.close()
        print(f"generate {outf}.")
    get_sa2obj(snATAC_sample)
