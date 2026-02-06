.PHONY: syncBioScala cleanScalaIndexing

bioscalaRoot := "/home/szu/git-recipes/bioscala"
syncBioscala:
	@rsync --progress -u .scalafmt.conf ${bioscalaRoot}/
	@rsync --progress -r -u bioscala/src/main/scala/* ${bioscalaRoot}/src/main/scala/
	@rsync --progress -r -u bioscala/src/test/scala/* ${bioscalaRoot}/src/test/scala/
	@rsync --progress -u project/build.properties ${bioscalaRoot}/build.properties

cleanScalaIndexing:
	-@rm -rf .bloop
	-@rm -rf .bsp
	-@rm -rf .metals
	-@rm -rf project/.bloop
	-@rm project/metals.sbt
	-@rm -rf target
	-@rm -rf project/project
	-@rm -rf project/target
	-@rm -rf bioscala/target
	-@rm -rf 100.project/target
	-@rm -rf 00.datapreprocess/target
	-@rm -rf 03.integration/target
	-@rm -rf 04.peakcalling/target
	-@rm -rf 05.CRE/target
	-@rm -rf 06.ChromHMM/target
	-@rm -rf 07.deeplearning/target
	-@rm -rf 09.epimem/target
	-@rm -rf 10.superEnhancer/target
	-@rm -rf 12.DE/target
	-@rm -rf 13.cicero/target
	-@rm -rf 17.repressiveMarks/target

projd := /projects/ps-renlab2/szu/projects/amb_pairedtag
google_drive := ${projd}/data/google_drive
h5ad_softln_path := ${google_drive}/h5ad

pt_h5ad_from := ${projd}/data/pairedtag_ann
pt_h5ad_to := ${h5ad_softln_path}/pairedtag_ann
add_pt_h5ad_softlink:
	ln -s ${pt_h5ad_from} ${pt_h5ad_to}

mouse_atac_h5ad_from := ${projd}/data/snATAC
mouse_atac_h5ad_to := ${h5ad_softln_path}/snATAC
add_mouseatac_h5ad_softlink:
	ln -s ${mouse_atac_h5ad_from}/all.sa2v26.adset.h5ad ${mouse_atac_h5ad_to}
	ln -s ${mouse_atac_h5ad_from}/all.sa2v26.ann.h5ad ${mouse_atac_h5ad_to}
	ln -s ${mouse_atac_h5ad_from}/all.sa2v26.ann.pmat.nofrag.h5ad ${mouse_atac_h5ad_to}
	ln -s ${mouse_atac_h5ad_from}/atac.ann.ptregion.h5ad ${mouse_atac_h5ad_to}
	ln -s ${mouse_atac_h5ad_from}/atac.sa2.pmat.113subclass.h5ad ${mouse_atac_h5ad_to}
	ln -s ${mouse_atac_h5ad_from}/snATAC.gmat.cnt.ann.h5ad ${mouse_atac_h5ad_to}
	ln -s ${mouse_atac_h5ad_from}/snATAC_gmat_obsv9.8_ann.h5ad ${mouse_atac_h5ad_to}
	ln -s ${mouse_atac_h5ad_from}/snATAC_gmat_obsv9.8_ptregion_ann.h5ad ${mouse_atac_h5ad_to}
	ln -s ${mouse_atac_h5ad_from}/snATAC_gmat_obsv9.8_ptregion_neu_ann.h5ad ${mouse_atac_h5ad_to}
	ln -s ${mouse_atac_h5ad_from}/snATAC_gmat_obsv9.8_ptregion_nn_ann.h5ad ${mouse_atac_h5ad_to}
  ln -s ${mouse_atac_h5ad_from}/snapatac2.6_samples ${mouse_atac_h5ad_to}









