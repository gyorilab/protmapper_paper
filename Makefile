OUTPUT := output
PLOTS := plots
DATA := data
DEPLOY := ../sitemapper_manuscript/figures/figure_panels

all: fig1 graph

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
$(OUTPUT)/pc_psp_modified_agents.pkl:
	python get_pc_sites.py

$(OUTPUT)/pc_sites_by_db.pkl: \
    $(OUTPUT)/pc_psp_modified_agents.pkl \
    $(OUTPUT)/pc_pid_modified_agents.pkl \
    $(OUTPUT)/pc_reactome_modified_agents.pkl
	python sitemap_fig.py map_pc_sites

# BEL Sites
$(OUTPUT)/large_corpus_pybel.pkl: $(DATA)/large_corpus.bel
	python process_bel_large_corpus.py parse_belscript

$(OUTPUT)/bel_mod_agents.pkl: $(OUTPUT)/large_corpus_pybel.pkl
	python process_bel_large_corpus.py get_pybel_mod_agents

$(OUTPUT)/bel_sites.pkl: $(OUTPUT)/bel_mod_agents.pkl
	python sitemap_fig.py map_bel_sites

$(OUTPUT)/all_db_sites.csv: \
    $(OUTPUT)/bel_sites.pkl \
    $(OUTPUT)/pc_sites_by_db.pkl
	python sitemap_fig.py create_site_csv


$(PLOTS)/site_stats_by_site.pdf: $(OUTPUT)/all_db_sites.csv
	python sitemap_fig.py plot_site_stats
