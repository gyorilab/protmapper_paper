# Assembling a corpus of phosphoproteomic annotations using ProtMapper to normalize site information from databases and text mining

This repository contains scripts and other files for generating the results
of the manuscript. The results build on both the ProtMapper and INDRA
systems.

ProtMapper is available at: https://github.com/indralab/protmapper
INDRA is available at: https://github.com/sorgerlab/indra

The `Makefile` contains all the commands for generating results based on
input data. The scripts used for processing are inside the `protmapper_paper`
package.

## Collecting and mapping phosphosite data

The `protmapper_paper/get_sites` module contains scripts to extract
phosphosite annotation data from various sources including:
- Machine reading-based annotations aggregated by INDRA using the Reach,
  Sparser and RLIMS-P systems.
- Annotations from the SIGNOR, HPRD, Reactome, NCI-PID, Panther and
  PhopsphoSitePlus databases using source processors implemented in INDRA.

The `protmapper_paper/map_sites.py` script uses the ProtMapper to
map all "raw" sites collected from the above sources.

## Analyzing phosphosite data

The following scripts in the `protmapper_paper` module perform analysis on the
aggregated and mapped phosphosite annotation data:
- `analyze_sites.py`: Statistics on the number of sites and annotations
  from different sources.
- `annotation_count.py`: Comparing the number of sites and annotations
  from different sources and the effect of applying ProtMapper mapping
  in the context of analyzing phosphoproteomic data.
- `mapping_stats.py`: Statistics on the effect of applying ProtMapper mapping
  to phosphosite annotations aggregated from different sources.
- `psp_reading_venn`: Create Venn diagrams of overlap between sites and
  annotations from different sources.
- `sample_reading_sites`: Script to sample sites from machine reading for
  curation.
