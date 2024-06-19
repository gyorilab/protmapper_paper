OUTPUT := output
PLOTS := plots
DATA := data
DEPLOY := ../protmapper_manuscript/figures/figure_panels

all: figs cptac export

figs: $(PLOTS)/site_stats_by_site.pdf \
      $(PLOTS)/psp_db_reader_sites_overlap_distinct.pdf \
      $(PLOTS)/reader_sites_overlap_distinct.pdf \
      $(PLOTS)/psp_db_reader_annotation_overlap_distinct.pdf \
      $(PLOTS)/reader_annotation_overlap_distinct.pdf \
      $(PLOTS)/psp_db_reader_annotation_overlap_distinct_kinase.pdf \
      $(PLOTS)/reader_annotation_overlap_distinct_kinase.pdf \
      $(PLOTS)/psp_db_reader_annotation_overlap_distinct_kinase_nofamplex.pdf \
      $(PLOTS)/reader_annotation_overlap_distinct_kinase_nofamplex.pdf \
      $(PLOTS)/annotations_valid_counts.pdf

cptac: \
    $(OUTPUT)/brca_peptide_mapping_results.csv \
    $(OUTPUT)/ovca_peptide_mapping_results.csv \
    $(OUTPUT)/brca_annotation_counts.csv \
    $(OUTPUT)/ovca_annotation_counts.csv

export: $(OUTPUT)/export.csv

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

$(DATA)/PathwayCommons12.reactome.BIOPAX.owl.gz:
	wget -P $(DATA) https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.reactome.BIOPAX.owl.gz

$(DATA)/PathwayCommons12.pid.BIOPAX.owl.gz:
	wget -P $(DATA) https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.pid.BIOPAX.owl.gz

$(DATA)/HPRD_FLAT_FILES_041310.tar.gz:
    wget -P $(DATA) https://rescued.omnipathdb.org/HPRD_FLAT_FILES_041310.tar.gz

# Get phospho statements from INDRA DB/Reading -----------------------
$(OUTPUT)/indra_db_stmts.pkl $(OUTPUT)/indra_agent_mod_stmts.pkl:
	python -m protmapper_paper.get_sites.indra_sites get_db_phos_stmts \
		$(OUTPUT)/indra_db_stmts.pkl $(OUTPUT)/indra_agent_mod_stmts.pkl

$(OUTPUT)/indra_phos_stmts_gmap_uniq_respos.pkl: $(OUTPUT)/indra_db_stmts.pkl
	python -m protmapper_paper.get_sites.indra_sites preprocess_stmts $< $@

$(OUTPUT)/indra_reach.sites.pkl: \
    $(OUTPUT)/indra_phos_stmts_gmap_uniq_respos.pkl
	python -m protmapper_paper.get_sites.indra_sites stmts_by_site $< reach $@

$(OUTPUT)/indra_rlimsp.sites.pkl: \
    $(OUTPUT)/indra_phos_stmts_gmap_uniq_respos.pkl
	python -m protmapper_paper.get_sites.indra_sites stmts_by_site $< rlimsp $@

$(OUTPUT)/indra_sparser.sites.pkl: \
    $(OUTPUT)/indra_phos_stmts_gmap_uniq_respos.pkl
	python -m protmapper_paper.get_sites.indra_sites stmts_by_site $< sparser $@

# Get modified Agent statements from INDRA DB/Reading -----------------------
$(OUTPUT)/indra_agent_mod_stmts_gmap_uniq_respos.pkl: $(OUTPUT)/indra_agent_mod_stmts.pkl
	python -m protmapper_paper.get_sites.indra_sites preprocess_stmts $< $@

$(OUTPUT)/indra_reach_agent_mod.sites.pkl: \
    $(OUTPUT)/indra_agent_mod_stmts_gmap_uniq_respos.pkl
	python -m protmapper_paper.get_sites.indra_sites agent_mod_stmts_by_site $< reach $@

$(OUTPUT)/indra_rlimsp_agent_mod.sites.pkl: \
    $(OUTPUT)/indra_agent_mod_stmts_gmap_uniq_respos.pkl
	python -m protmapper_paper.get_sites.indra_sites agent_mod_stmts_by_site $< rlimsp $@

$(OUTPUT)/indra_sparser_agent_mod.sites.pkl: \
    $(OUTPUT)/indra_agent_mod_stmts_gmap_uniq_respos.pkl
	python -m protmapper_paper.get_sites.indra_sites agent_mod_stmts_by_site $< sparser $@


# COLLECT_SITES --------------------------------------------------------------
# PID sites
$(OUTPUT)/pid.sites.pkl: $(DATA)/PathwayCommons12.pid.BIOPAX.owl.gz
	python -m protmapper_paper.get_sites.biopax $< $@

# Reactome sites
$(OUTPUT)/reactome.sites.pkl: $(DATA)/PathwayCommons12.reactome.BIOPAX.owl.gz
	python -m protmapper_paper.get_sites.biopax $< $@

# PSP sites
$(OUTPUT)/psp.sites.pkl: $(DATA)/Kinase_substrates.owl.gz
	python -m protmapper_paper.get_sites.biopax $< $@

# BEL Sites
$(OUTPUT)/bel_large_corpus.sites.pkl: $(DATA)/large_corpus.bel
	python -m protmapper_paper.get_sites.bel $< $@

# SIGNOR sites
$(OUTPUT)/signor.sites.pkl:
	python -m protmapper_paper.get_sites.signor $@

# HPRD sites
$(OUTPUT)/hprd.sites.pkl: $(DATA)/HPRD_FLAT_FILES_041310.tar.gz
	python -m protmapper_paper.get_sites.hprd $< $@

# All sites, combined into a single dict
$(OUTPUT)/all_sites.pkl: \
    $(OUTPUT)/signor.sites.pkl \
    $(OUTPUT)/hprd.sites.pkl \
    $(OUTPUT)/bel_large_corpus.sites.pkl \
    $(OUTPUT)/pid.sites.pkl \
    $(OUTPUT)/reactome.sites.pkl \
    $(OUTPUT)/psp.sites.pkl \
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
	python -m protmapper_paper.analyze_sites create_site_csv \
        $< $(word 2,$^) $@ $(OUTPUT)/annotations.csv

$(OUTPUT)/annotations.csv: $(OUTPUT)/site_info.csv

# Plots on correctness/mappability
$(PLOTS)/site_stats_by_site.pdf: $(OUTPUT)/site_info.csv
	python -m protmapper_paper.analyze_sites plot_site_stats $< \
        $(PLOTS)/site_stats

$(PLOTS)/psp_db_reader_sites_overlap_distinct.pdf $(PLOTS)/reader_sites_overlap_distinct.pdf $(PLOTS)/psp_db_reader_annotation_overlap_distinct.pdf $(PLOTS)/reader_annotation_overlap_distinct.pdf $(PLOTS)/psp_db_reader_annotation_overlap_distinct_kinase.pdf $(PLOTS)/reader_annotation_overlap_distinct_kinase.pdf $(PLOTS)/psp_db_reader_annotation_overlap_distinct_kinase_nofamplex.pdf $(PLOTS)/reader_annotation_overlap_distinct_kinase_nofamplex.pdf: \
    $(OUTPUT)/site_info.csv \
    $(OUTPUT)/annotations.csv
	python -m protmapper_paper.psp_reading_venn

# Plots on correctness/mappability
$(PLOTS)/annotations_valid_counts.pdf: $(OUTPUT)/annotations.csv
	python -m protmapper_paper.analyze_sites plot_annot_stats $< \
        $(PLOTS)/annotations


# SITE SAMPLE ----------------------------------------------
# Sample of sites for curation
$(OUTPUT)/site_sample.csv: $(OUTPUT)/site_info.csv
	python -m protmapper_paper.analyze_sites site_samples $< $@


# BRCA and OVCA peptide mapping ---------------------------------------------
$(OUTPUT)/brca_peptide_mapping_results.csv: \
    $(DATA)/CPTAC2_Breast_Prospective_Collection_BI_Phosphoproteome.phosphosite.tmt10.tsv
	python -m protmapper_paper.mapping_stats $< $@

$(OUTPUT)/ovca_peptide_mapping_results.csv: \
    $(DATA)/TCGA_Ovarian_PNNL_Phosphoproteome.phosphosite.itraq.tsv
	python -m protmapper_paper.mapping_stats $< $@

$(OUTPUT)/brca_annotation_counts.csv: \
    $(OUTPUT)/annotations.csv \
    $(DATA)/CPTAC2_Breast_Prospective_Collection_BI_Phosphoproteome.phosphosite.tmt10.tsv
	python -m protmapper_paper.annotation_count $< $(word 2,$^) $@

$(OUTPUT)/ovca_annotation_counts.csv: \
    $(OUTPUT)/annotations.csv \
    $(DATA)/TCGA_Ovarian_PNNL_Phosphoproteome.phosphosite.itraq.tsv
	python -m protmapper_paper.annotation_count $< $(word 2,$^) $@

# Statements, beliefs and export
$(OUTPUT)/annotation_statements.pkl: \
    $(OUTPUT)/all_sites.pkl \
    $(OUTPUT)/mapping_results.pkl
	python -m protmapper_paper.annotation_stmts \
		$(OUTPUT)/all_sites.pkl $(OUTPUT)/mapping_results.pkl $@


$(OUTPUT)/stmt_beliefs.json: \
    $(OUTPUT)/annotation_statements.pkl \
    $(DATA)/protmapper_belief_training_corpus.json \
    $(DATA)/protmapper_belief_training_corpus_curations.json \
    $(DATA)/protmapper_belief_training_corpus_extra_evidences.pkl
	python -m protmapper_paper.belief \
		    $(OUTPUT)/annotation_statements.pkl \
		    $(DATA)/protmapper_belief_training_corpus.json \
		    $(DATA)/protmapper_belief_training_corpus_curations.json \
		    $(DATA)/protmapper_belief_training_corpus_extra_evidences.pkl $@

# Export of all annotated sites with evidence
$(OUTPUT)/export.csv: \
    $(OUTPUT)/annotation_statements.pkl
	python -m protmapper_paper.analyze_sites export \
        $(OUTPUT)/annotation_statements.pkl $(OUTPUT)/stmt_beliefs.json $@ $(OUTPUT)/evidences.csv

