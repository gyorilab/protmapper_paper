OUTPUT := output
PLOTS := plots
DATA := data
DEPLOY := ../protmapper_manuscript/figures/figure_panels

all: fig1 brca mouse

fig1: indra_sites \
      $(PLOTS)/site_stats_by_site.pdf \
      $(PLOTS)/psp_reader_site_overlap.pdf

indra_sites: $(OUTPUT)/indra_stmts_by_site.pkl

brca: $(OUTPUT)/brca_site_stats.txt

mouse: $(OUTPUT)/psp_relations_by_site.pkl
#$(OUTPUT)/mouse_kin_sub_count.txt

# MAKEFILE GRAPH
graph: makegraph.pdf

makegraph.pdf: makegraph.dot
	dot -T pdf makegraph.dot -o makegraph.pdf

makegraph.dot: Makefile
	make -Bnd | make2graph > makegraph.dot

deploy:
	rsync -av $(PLOTS)/*.pdf $(DEPLOY)

clean:
	cd $(OUTPUT); rm -rf *


# DATA -----------------------------------------------------------------------

#$(DATA)/PathwayCommons9.All.hgnc.txt:
#	wget -P $(DATA) http://www.pathwaycommons.org/archives/PC2/v9/PathwayCommons9.All.hgnc.txt.gz
#	gunzip $@
#
# Get phospho statements from INDRA DB/Reading -----------------------
$(OUTPUT)/indra_phos_stmts.pkl:
	python get_indra_sites.py get_phos_stmts $@

$(OUTPUT)/indra_phos_stmts_gmap_uniq_respos.pkl: $(OUTPUT)/indra_phos_stmts.pkl
	python get_indra_sites.py preprocess_stmts $< $@

$(OUTPUT)/indra_stmts_by_site.pkl: \
    $(OUTPUT)/indra_phos_stmts_gmap_uniq_respos.pkl
	python get_indra_sites.py stmts_by_site $< $@


# FIG1 --------------------------------------------------------------
# PC Sites/Biopax
$(OUTPUT)/%.sites.pkl: $(DATA)/biopax/%.owl
	python -m protmapper_paper.get_sites.biopax $< $@

#$(OUTPUT)/psp_kinase_substrate_tsv.pkl: \
#    $(DATA)/Kinase_Substrate_Dataset
#	python -m protmapper_paper.get_sites.psp $< $@


# Skipping the other sources for PSP and Reactome data
#$(OUTPUT)/biopax/PathwayCommons10.psp.BIOPAX.pkl
#$(OUTPUT)/biopax/Homo_sapiens.pkl

# BEL Sites
$(OUTPUT)/large_corpus_pybel.pkl: $(DATA)/large_corpus.bel
	python -m protmapper_paper.get_sites.bel parse_belscript $< $@

#$(OUTPUT)/bel_mod_agents.pkl: $(OUTPUT)/large_corpus_pybel.pkl
#	python -m protmapper_paper.get_sites.bel get_pybel_mod_agents

$(OUTPUT)/bel_large_corpus.sites.pkl: $(OUTPUT)/large_corpus_pybel.pkl
	python -m protmapper_paper.get_sites.bel get_pybel_stmts_by_site $< $@

$(OUTPUT)/signor.sites.pkl:
	python -m protmapper_paper.get_sites.signor $@

#$(OUTPUT)/bel_sites.pkl: $(OUTPUT)/bel_mod_agents.pkl
#	python sitemap_fig.py map_bel_sites > /dev/null


#$(OUTPUT)/biopax_sites_by_db.pkl: \
#    $(OUTPUT)/PathwayCommons10.hprd.BIOPAX.sites.pkl \
#    $(OUTPUT)/PathwayCommons10.kegg.BIOPAX.sites.pkl \
#    $(OUTPUT)/PathwayCommons10.panther.BIOPAX.sites.pkl \
#    $(OUTPUT)/PathwayCommons10.pid.BIOPAX.sites.pkl \
#    $(OUTPUT)/PathwayCommons10.reactome.BIOPAX.sites.pkl \
#    $(OUTPUT)/PathwayCommons10.wp.BIOPAX.sites.pkl \
#    $(OUTPUT)/psp_kinase_substrate_tsv.sites.pkl \
#    $(OUTPUT)/biopax/Kinase_substrates.sites.pkl
#	python map_pc_sites > /dev/null


# Reader Sites
$(OUTPUT)/reader_sites.pkl: $(OUTPUT)/indra_phos_stmts_gmap_uniq_respos.pkl
	python get_indra_sites.py reader_sites $<

# All sites, combined into a single dict
$(OUTPUT)/all_sites.pkl: \
    $(OUTPUT)/signor.sites.pkl \
    $(OUTPUT)/bel_large_corpus.sites.pkl \
    $(OUTPUT)/PathwayCommons10.hprd.BIOPAX.sites.pkl \
    $(OUTPUT)/PathwayCommons10.kegg.BIOPAX.sites.pkl \
    $(OUTPUT)/PathwayCommons10.panther.BIOPAX.sites.pkl \
    $(OUTPUT)/Kinase_substrates.sites.pkl
	python -m protmapper_paper.get_sites.combine $@ $(OUTPUT)/*.sites.pkl

#$(OUTPUT)/PathwayCommons10.reactome.BIOPAX.sites.pkl \
$(OUTPUT)/PathwayCommons10.pid.BIOPAX.sites.pkl \
$(OUTPUT)/PathwayCommons10.wp.BIOPAX.sites.pkl \


# All sites combined into a single dataframe
$(OUTPUT)/all_db_sites.csv: \
    $(OUTPUT)/bel_sites.pkl \
    $(OUTPUT)/biopax_sites_by_db.pkl \
    $(OUTPUT)/reader_sites.pkl
	python sitemap_fig.py create_site_csv

# Plots on correctness/mappability
$(PLOTS)/site_stats_by_site.pdf: $(OUTPUT)/all_db_sites.csv
	python sitemap_fig.py plot_site_stats

$(PLOTS)/psp_reader_site_overlap.pdf: \
    $(OUTPUT)/reader_sites.pkl \
    $(OUTPUT)/psp_kinase_substrate_tsv.pkl
	python psp_reading_venn.py

# BRCA data ----------------------------------------------------------
$(OUTPUT)/brca_up_mappings.txt: \
    $(DATA)/breast_phosphosites.txt \
    $(DATA)/HUMAN_9606_idmapping.dat
	python brca_data.py map_uniprot

$(OUTPUT)/brca_site_stats.txt: \
    $(OUTPUT)/indra_stmts_by_site.pkl \
    $(OUTPUT)/brca_up_mappings.txt
	python brca_data.py site_stats $< $@



# MOUSE data ---------------------------------------------------------

$(OUTPUT)/psp_relations_by_site.pkl: $(DATA)/Kinase_Substrate_Dataset
	python get_psp_tsv_sites.py

#$(OUTPUT)/mouse_kin_sub_count.txt: $(OUTPUT)/psp_relations_by_site.pkl

