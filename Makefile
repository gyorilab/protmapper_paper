OUTPUT := output
PLOTS := plots
DATA := data
DEPLOY := ../sitemapper_manuscript/figures/figure_panels

all: indra_sites fig1 graph

indra_sites: $(OUTPUT)/indra_phos_stmts_gmap_uniq_respos.pkl

fig1: $(PLOTS)/site_stats_by_site.pdf

# MAKEFILE GRAPH
graph: makegraph.pdf

makegraph.pdf: makegraph.dot
	dot -T pdf makegraph.dot -o makegraph.pdf

makegraph.dot: Makefile
	make -Bnd | make2graph > makegraph.dot

deploy:
	rsync -av $(OUTPUT)/*.pdf $(DEPLOY)

clean:
	cd $(OUTPUT); rm -rf *


# DATA -----------------------------------------------------------------------

#$(DATA)/PathwayCommons9.All.hgnc.txt:
#	wget -P $(DATA) http://www.pathwaycommons.org/archives/PC2/v9/PathwayCommons9.All.hgnc.txt.gz
#	gunzip $@

# FIG1 --------------------------------------------------------------
# PC Sites
$(OUTPUT)/pc_psp_modified_agents.pkl: \
    $(DATA)/pc/PathwayCommons10.psp.BIOPAX.owl \
    $(DATA)/pc/PathwayCommons10.reactome.BIOPAX.owl \
    $(DATA)/pc/PathwayCommons10.pid.BIOPAX.owl
	python get_pc_sites.py

$(OUTPUT)/pc_sites_by_db.pkl: \
    $(OUTPUT)/pc_psp_modified_agents.pkl \
    $(OUTPUT)/pc_pid_modified_agents.pkl \
    $(OUTPUT)/pc_reactome_modified_agents.pkl
	python sitemap_fig.py map_pc_sites > /dev/null

# BEL Sites
$(OUTPUT)/large_corpus_pybel.pkl: $(DATA)/large_corpus.bel
	python process_bel_large_corpus.py parse_belscript 2> /dev/null

$(OUTPUT)/bel_mod_agents.pkl: $(OUTPUT)/large_corpus_pybel.pkl
	python process_bel_large_corpus.py get_pybel_mod_agents

$(OUTPUT)/bel_sites.pkl: $(OUTPUT)/bel_mod_agents.pkl
	python sitemap_fig.py map_bel_sites > /dev/null

$(OUTPUT)/all_db_sites.csv: \
    $(OUTPUT)/bel_sites.pkl \
    $(OUTPUT)/pc_sites_by_db.pkl
	python sitemap_fig.py create_site_csv

$(PLOTS)/site_stats_by_site.pdf: $(OUTPUT)/all_db_sites.csv
	python sitemap_fig.py plot_site_stats

# Get phospho statements from INDRA DB/Reading -----------------------
$(OUTPUT)/indra_phos_stmts.pkl:
	python get_db_sites.py get_phos_stmts $@

$(OUTPUT)/indra_phos_stmts_gmap_uniq_respos.pkl: $(OUTPUT)/indra_phos_stmts.pkl
	python get_db_sites.py preprocess_stmts $< $@


