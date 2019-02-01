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
	python get_db_sites.py get_phos_stmts $@

$(OUTPUT)/indra_phos_stmts_gmap_uniq_respos.pkl: $(OUTPUT)/indra_phos_stmts.pkl
	python get_db_sites.py preprocess_stmts $< $@

$(OUTPUT)/indra_stmts_by_site.pkl: \
    $(OUTPUT)/indra_phos_stmts_gmap_uniq_respos.pkl
	python get_db_sites.py stmts_by_site $< $@


# FIG1 --------------------------------------------------------------
# PC Sites
$(OUTPUT)/biopax/%.pkl: $(DATA)/biopax/%.owl
	python get_pc_sites.py $< $@

$(OUTPUT)/psp_kinase_substrate_tsv.pkl: \
    $(DATA)/Kinase_Substrate_Dataset
	python get_psp_tsv_sites.py

$(OUTPUT)/biopax_sites_by_db.pkl: \
    $(OUTPUT)/biopax/PathwayCommons10.hprd.BIOPAX.pkl \
    $(OUTPUT)/biopax/PathwayCommons10.kegg.BIOPAX.pkl \
    $(OUTPUT)/biopax/PathwayCommons10.panther.BIOPAX.pkl \
    $(OUTPUT)/biopax/PathwayCommons10.pid.BIOPAX.pkl \
    $(OUTPUT)/biopax/PathwayCommons10.psp.BIOPAX.pkl \
    $(OUTPUT)/biopax/PathwayCommons10.reactome.BIOPAX.pkl \
    $(OUTPUT)/biopax/PathwayCommons10.wp.BIOPAX.pkl \
    $(OUTPUT)/biopax/Kinase_substrates.pkl \
    $(OUTPUT)/psp_kinase_substrate_tsv.pkl \
    $(OUTPUT)/biopax/Homo_sapiens.pkl
	python sitemap_fig.py map_pc_sites > /dev/null

# BEL Sites
$(OUTPUT)/large_corpus_pybel.pkl: $(DATA)/large_corpus.bel
	python process_bel_large_corpus.py parse_belscript 2> /dev/null

$(OUTPUT)/bel_mod_agents.pkl: $(OUTPUT)/large_corpus_pybel.pkl
	python process_bel_large_corpus.py get_pybel_mod_agents

$(OUTPUT)/bel_sites.pkl: $(OUTPUT)/bel_mod_agents.pkl
	python sitemap_fig.py map_bel_sites > /dev/null

# Reader Sites
$(OUTPUT)/reader_sites.pkl: $(OUTPUT)/indra_phos_stmts_gmap_uniq_respos.pkl
	python get_db_sites.py reader_sites $<

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

