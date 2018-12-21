OUTPUT := output
DATA := data
DEPLOY := ../sitemapper_manuscript/figures/figure_panels

fig1: $(OUTPUT)/fig1_pct_incorrect_sites.pdf

deploy:
	rsync -av $(OUTPUT)/*.pdf $(DEPLOY)

clean:
	cd $(OUTPUT); rm -rf *


# DATA -----------------------------------------------------------------------

#$(DATA)/PathwayCommons9.All.hgnc.txt:
#	wget -P $(DATA) http://www.pathwaycommons.org/archives/PC2/v9/PathwayCommons9.All.hgnc.txt.gz
#	gunzip $@

# FIG1 --------------------------------------------------------------
$(OUTPUT)/pc_psp_modified_agents.pkl:
	python get_db_sites.py

$(OUTPUT)/pc_sites_by_db.pkl:
	python sitemap_fig.py map_pc_sites

$(OUTPUT)/large_corpus_stmts.pkl: $(DATA)/large_corpus.bel
	python process_bel_large_corpus.py $(DATA)/large_corpus.bel

$(OUTPUT)/fig1_pct_incorrect_sites.pdf: \
    $(OUTPUT)/large_corpus_stmts.pkl \
    $(OUTPUT)/pc_sites_by_db.pkl
	echo "Hello"
