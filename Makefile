OUTPUT := output
PLOTS := plots
DATA := data
DEPLOY := ../protmapper_manuscript/figures/figure_panels

all: fig1 brca mouse

fig1: $(PLOTS)/site_stats_by_site.pdf \
      $(PLOTS)/psp_reader_site_overlap.pdf


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


agent_mods: $(OUTPUT)/indra_sparser_agent_mod.sites.pkl \
            $(OUTPUT)/indra_reach_agent_mod.sites.pkl \
            $(OUTPUT)/indra_rlimsp_agent_mod.sites.pkl

# DATA -----------------------------------------------------------------------

#$(DATA)/PathwayCommons9.All.hgnc.txt:
#	wget -P $(DATA) http://www.pathwaycommons.org/archives/PC2/v9/PathwayCommons9.All.hgnc.txt.gz
#	gunzip $@
#
# Get phospho statements from INDRA DB/Reading -----------------------
$(OUTPUT)/indra_db_stmts.pkl:
	python -m protmapper_paper.get_sites.indra get_db_phos_stmts $@

$(OUTPUT)/indra_rlimsp_stmts.pkl: $(DATA)/rlims.medline.json $(DATA)/rlims.pmc.json
	python -m protmapper_paper.get_sites.indra get_rlimsp_phos_stmts $@ $(DATA)/rlims.medline.json $(DATA)/rlims.pmc.json

$(OUTPUT)/indra_all_stmts.pkl: $(OUTPUT)/indra_db_stmts.pkl $(OUTPUT)/indra_rlimsp_stmts.pkl
	python -m protmapper_paper.get_sites.indra get_all_indra_phos_stmts $@ $(OUTPUT)/indra_db_stmts.pkl $(OUTPUT)/indra_rlimsp_stmts.pkl

$(OUTPUT)/indra_phos_stmts_gmap_uniq_respos.pkl: $(OUTPUT)/indra_all_stmts.pkl
	python -m protmapper_paper.get_sites.indra preprocess_stmts $< $@ true

$(OUTPUT)/indra_reach.sites.pkl: \
    $(OUTPUT)/indra_phos_stmts_gmap_uniq_respos.pkl
	python -m protmapper_paper.get_sites.indra stmts_by_site $< reach $@

$(OUTPUT)/indra_rlimsp.sites.pkl: \
    $(OUTPUT)/indra_phos_stmts_gmap_uniq_respos.pkl
	python -m protmapper_paper.get_sites.indra stmts_by_site $< rlimsp $@

$(OUTPUT)/indra_sparser.sites.pkl: \
    $(OUTPUT)/indra_phos_stmts_gmap_uniq_respos.pkl
	python -m protmapper_paper.get_sites.indra stmts_by_site $< sparser $@

# Get modified Agent statements from INDRA DB/Reading -----------------------
$(OUTPUT)/indra_agent_mod_stmts.pkl:
	python -m protmapper_paper.get_sites.indra get_agent_mod_stmts $@

$(OUTPUT)/indra_agent_mod_stmts_gmap_uniq_respos.pkl: $(OUTPUT)/indra_agent_mod_stmts.pkl
	python -m protmapper_paper.get_sites.indra preprocess_stmts $< $@ false

$(OUTPUT)/indra_reach_agent_mod.sites.pkl: \
    $(OUTPUT)/indra_agent_mod_stmts_gmap_uniq_respos.pkl
	python -m protmapper_paper.get_sites.indra agent_mod_stmts_by_site $< reach $@

$(OUTPUT)/indra_rlimsp_agent_mod.sites.pkl: \
    $(OUTPUT)/indra_agent_mod_stmts_gmap_uniq_respos.pkl
	python -m protmapper_paper.get_sites.indra agent_mod_stmts_by_site $< rlimsp $@

$(OUTPUT)/indra_sparser_agent_mod.sites.pkl: \
    $(OUTPUT)/indra_agent_mod_stmts_gmap_uniq_respos.pkl
	python -m protmapper_paper.get_sites.indra agent_mod_stmts_by_site $< sparser $@


# COLLECT_SITES --------------------------------------------------------------
# PC Sites/Biopax
$(OUTPUT)/%.sites.pkl: $(DATA)/biopax/%.owl
	python -m protmapper_paper.get_sites.biopax $< $@

# BEL Sites
$(OUTPUT)/large_corpus_pybel.pkl: $(DATA)/large_corpus.bel
	python -m protmapper_paper.get_sites.bel parse_belscript $< $@

$(OUTPUT)/bel_large_corpus.sites.pkl: $(OUTPUT)/large_corpus_pybel.pkl
	python -m protmapper_paper.get_sites.bel get_pybel_stmts_by_site $< $@

# SIGNOR sites
$(OUTPUT)/signor.sites.pkl:
	python -m protmapper_paper.get_sites.signor $@

# Reader Sites
$(OUTPUT)/reader_sites.pkl: $(OUTPUT)/indra_phos_stmts_gmap_uniq_respos.pkl
	python get_indra_sites.py reader_sites $<

# Phosphosite from TSV
#$(OUTPUT)/psp_kinase_substrate_tsv.pkl: \
#    $(DATA)/Kinase_Substrate_Dataset
#	python -m protmapper_paper.get_sites.psp $< $@

# All sites, combined into a single dict
# Skipping the other sources for PSP and Reactome data
#$(OUTPUT)/biopax/PathwayCommons10.psp.BIOPAX.pkl
#$(OUTPUT)/biopax/Homo_sapiens.pkl
$(OUTPUT)/all_sites.pkl: \
    $(OUTPUT)/signor.sites.pkl \
    $(OUTPUT)/bel_large_corpus.sites.pkl \
    $(OUTPUT)/PathwayCommons10.hprd.BIOPAX.sites.pkl \
    $(OUTPUT)/PathwayCommons10.kegg.BIOPAX.sites.pkl \
    $(OUTPUT)/PathwayCommons10.panther.BIOPAX.sites.pkl \
    $(OUTPUT)/PathwayCommons10.pid.BIOPAX.sites.pkl \
    $(OUTPUT)/PathwayCommons10.reactome.BIOPAX.sites.pkl \
    $(OUTPUT)/PathwayCommons10.wp.BIOPAX.sites.pkl \
    $(OUTPUT)/Kinase_substrates.sites.pkl \
    $(OUTPUT)/indra_reach.sites.pkl \
    $(OUTPUT)/indra_sparser.sites.pkl \
    $(OUTPUT)/indra_rlimsp.sites.pkl \
    $(OUTPUT)/indra_reach_agent_mod.sites.pkl \
    $(OUTPUT)/indra_sparser_agent_mod.sites.pkl \
    $(OUTPUT)/indra_rlimsp_agent_mod.sites.pkl
	python -m protmapper_paper.get_sites.combine $@ $(OUTPUT)/*.sites.pkl

$(OUTPUT)/mapping_results.pkl: $(OUTPUT)/all_sites.pkl
	python -m protmapper_paper.map_sites $< $@

# All sites combined into a single dataframe
$(OUTPUT)/site_info.csv: \
    $(OUTPUT)/all_sites.pkl \
    $(OUTPUT)/mapping_results.pkl
	python -m protmapper_paper.analyze_sites create_site_csv $< $(word 2,$^) $@

# Plots on correctness/mappability
$(PLOTS)/site_stats_by_site.pdf: $(OUTPUT)/site_info.csv
	python -m protmapper_paper.analyze_sites plot_site_stats $< \
        $(PLOTS)/site_stats

#$(PLOTS)/psp_reader_site_overlap.pdf: \
#    $(OUTPUT)/reader_sites.pkl \
#    $(OUTPUT)/psp_kinase_substrate_tsv.pkl
#	python psp_reading_venn.py


# SITE SAMPLE ----------------------------------------------
# Plots on correctness/mappability
$(OUTPUT)/site_sample.csv: $(OUTPUT)/site_info.csv
	python -m protmapper_paper.analyze_sites site_samples $< $@


# BRCA data ----------------------------------------------------------
$(OUTPUT)/brca_up_mappings.txt: \
    $(DATA)/breast_phosphosites.txt \
    $(DATA)/HUMAN_9606_idmapping.dat
	python -m protmapper_paper.brca_data map_uniprot $< $(word 2,$^) $@


$(OUTPUT)/brca_site_stats.txt: \
    $(OUTPUT)/all_sites.pkl \
    $(OUTPUT)/brca_up_mappings.txt
	python -m protmapper_paper.brca_data site_stats $< $(word 2,$^) $@


# MOUSE data ---------------------------------------------------------

$(OUTPUT)/psp_relations_by_site.pkl: $(DATA)/Kinase_Substrate_Dataset
	python get_psp_tsv_sites.py

#$(OUTPUT)/mouse_kin_sub_count.txt: $(OUTPUT)/psp_relations_by_site.pkl

