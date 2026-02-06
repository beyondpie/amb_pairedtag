.PHONY: clean test

python := /home/szu/miniforge3/envs/sa2/bin/python
projd := /projects/ps-renlab2/szu/projects/amb_pairedtag
nmf_script := ${projd}/extratools/singlecell/nmf/nmf.py
nmfstat_script := ${projd}/extratools/singlecell/nmf/nmf.stat.py
nmf_outd := ${projd}/05.CRE/out/nmf

cpmd := ${projd}/data/ptHistoneCounts/ATACPeak_rds
CPM_pbyc_H3K27ac := ${cpmd}/cpm.scbyp.H3K27ac.h5
peakfnm := ${cpmd}/peak.sc.H3K27ac.txt
scfnm := ${cpmd}/subclass.H3K27ac.txt

nmfd := ${projd}/05.CRE/out/nmf

.PHONY: H3K27ac_nmf
H3K27ac_nmf:
	${python} ${nmf_script} --pbych5 ${CPM_pbyc_H3K27ac} \
      --seed 0 --rank 150 --outPrefix ${nmf_outd}/H3K27ac_nmf150

H3K27ac_nmfstat:
	${python} ${nmfstat_script} -p ${CPM_pbyc_H3K27ac} \
      -x ${peakfnm} -y ${scfnm} \
      --basis ${nmfd}/H3K27ac_nmf150.W.mx \
      --coef ${nmfd}/H3K27ac_nmf150.H.mx \
      -c 0.2 -o ${nmfd}/nmfPmat.H3K27ac.r150.n0


.PHONY: addpmat_snATACpeak addpmat_H3K27acpeak addpmat_test
coded := ${projd}/05.CRE/src/main/python
addpmat_script := ${coded}/sup.genpmat.sa2v26.snATAC.py
snATAC_peak := ${projd}/data/snATAC/wmb.snATAC.peak.srt.bed
H3K27ac_peak := ${projd}/04.peakcalling/out/mergePeak/H3K27ac.merge.all.blv2.me.peak.bed

addpmat_snATACpeak:
	${python} ${addpmat_script} ${snATAC_peak} snATAC

addpmat_H3K27acpeak:
	${python} ${addpmat_script} ${H3K27ac_peak} H3K27ac


test_peak := ${projd}/05.CRE/src/test/resource/test.peak.bed
addpmat_test:
	${python} ${addpmat_script} ${test_peak} test
