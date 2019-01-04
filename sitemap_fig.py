from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import sys
import glob
import pickle
from os.path import join as pjoin
from collections import OrderedDict, Counter, namedtuple, defaultdict
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('agg')
import seaborn as sns
from matplotlib import pyplot as plt
from indra.tools import assemble_corpus as ac
from indra.util import write_unicode_csv, read_unicode_csv
from indra.util import plot_formatting as pf
from indra.tools import assemble_corpus as ac

from protmapper import ProtMapper

# Constants
CACHE_PATH = 'output/pc_site_cache.pkl'
BIOPAX_SITES_BY_DB = 'output/biopax_sites_by_db.pkl'
BEL_AGENTS = 'output/bel_mod_agents.pkl'
BEL_SITES = 'output/bel_sites.pkl'
READER_SITES = 'output/reader_sites.pkl'
ALL_SITES_CSV = 'output/all_db_sites.csv'


def map_statements(stmts, source, outfile=None):
    """Tabulate valid, invalid, and mapped sites from a set of Statements."""
    # Look for errors in database statements
    sm = SiteMapper(default_site_map)
    valid_stmts, mapped_stmts = sm.map_sites(stmts)
    # Collect stats from SiteMapper itself
    sites = []
    for site_key, mapping in sm._cache.items():
        gene, res, pos = site_key
        freq = sm._sitecount[site_key]
        if mapping == 'VALID':
            valid, mapped, mapped_res, mapped_pos, explanation = \
                                                      (1, 0, None, None, None)
        else:
            valid = 0
            # Not mapped
            if mapping is None:
                mapped, mapped_res, mapped_pos, explanation = \
                                                    (0, None, None, None)
            # Mapped!
            else:
                mapped_res, mapped_pos, explanation = mapping
                mapped = 1 if mapped_pos else 0
        si = SiteInfo(gene, res, pos, valid, mapped, mapped_res, mapped_pos,
                      explanation, freq, source)
        sites.append(si)
    # Write to CSV file
    if outfile:
        header = [[field.upper() for field in si._asdict().keys()]]
        rows = header + replace_nones(sites)
        write_unicode_csv(outfile, rows)
    return sites


def map_agents(mod_agents_file, pm, source, save_csv=True):
    """Tabulate valid, invalid, and mapped sites from a set of Agents."""
    # Load the agents
    with open(mod_agents_file, 'rb') as f:
        mod_agents = pickle.load(f)
    print("Mapping %s" % mod_agents_file)
    sites = []
    for ag_ix, ag in enumerate(mod_agents):
        if ag_ix % 1000 == 0:
            print('%d of %d' % (ag_ix, len(mod_agents)))
        if ag is None:
            print("Skipping None agent")
            continue
        try:
            up_id = ag.db_refs.get('UP')
            if len(ag.mods) == 1:
                ms = pm.map_to_human_ref(up_id, 'uniprot', ag.mods[0].residue,
                                     ag.mods[0].position)
                sites.append(ms)
            elif len(ag.mods) > 1:
                sitelist = [(up_id, 'uniprot', m.residue, m.position)
                             for m in ag.mods]
                ms_list = pm.map_site_list_to_human_ref(sitelist)
                sites.extend(ms_list)
        except Exception as e:
            print("Error: %s" % str(e))
            print("agent: %s, up_id: %s, res %s, pos %s, db_refs %s" %
                  (ag, up_id, ag.mods[0].residue, ag.mods[0].position,
                   str(ag.db_refs)))
    # Now that we've collected a list of all the sites, tabulate frequencies
    site_counter = Counter(sites)
    return list(site_counter.items())


# UTIL for tabulating and saving Site Info --------------------------------

SiteInfo = namedtuple('SiteInfo', ['source', 'gene_name', 'up_id',
                                   'error_code',
                                   'valid', 'orig_res', 'orig_pos',
                                   'mapped_res', 'mapped_pos',
                                   'description', 'freq'])


def ms_to_si(source, freq, ms):
    return SiteInfo(source=source, gene_name=ms.gene_name,
                  up_id=ms.up_id, error_code=ms.error_code,
                  valid=ms.valid, orig_res=ms.orig_res, orig_pos=ms.orig_pos,
                  mapped_res=ms.mapped_res, mapped_pos=ms.mapped_pos,
                  description=ms.description, freq=freq)


def replace_nones(rows):
    return [[cell if cell is not None else '' for cell in row]
             for row in rows]


def print_stats(site_df):
    """Print statistics about site validity."""
    def pct(n, d):
        return 100 * n / float(d)
    df = site_df[site_df.ERROR_CODE.isna()]
    results = {}
    # Group statistics by source
    sources = df.SOURCE.unique()
    for source in sources:
        print("Stats for %s -------------" % source)
        s = df[df.SOURCE == source]
        s_valid = s[s.VALID]
        s_map = s[s.MAPPED_POS.notnull()]
        n = len(s)
        n_val = len(s_valid)
        n_inv = n - n_val
        n_map = len(s_map)
        n_unmap = n_inv - n_map
        f = s.FREQ.sum()
        f_val = s_valid.FREQ.sum()
        f_inv = f - f_val
        f_map = s_map.FREQ.sum()
        f_unmap = f_inv - f_map
        db_results = {
                'Source': source,
                'Total Sites': n,
                'Valid Sites': n_val,
                'Valid Sites Pct.': pct(n_val, n),
                'Invalid Sites': n_inv,
                'Invalid Sites Pct.': pct(n_inv, n),
                'Mapped Sites': n_map,
                'Mapped Sites Pct.': pct(n_map, n_inv),
                'Mapped Sites Pct. Total': pct(n_map, n),
                'Unmapped Sites': n_unmap,
                'Unmapped Sites Pct.': pct(n_unmap, n_inv),
                'Unmapped Sites Pct. Total': pct(n_unmap, n),
                'Total Occurrences': f,
                'Valid Occ.': f_val,
                'Valid Occ. Pct.': pct(f_val, f),
                'Invalid Occ.': f_inv,
                'Invalid Occ. Pct.': pct(f_inv, f),
                'Mapped Occ.': f_map,
                'Mapped Occ. Pct.': pct(f_map, f_inv),
                'Mapped Occ. Pct. Total': pct(f_map, f),
                'Unmapped Occ.': f_unmap,
                'Unmapped Occ. Pct.': pct(f_unmap, f_inv),
                'Unmapped Occ. Pct. Total': pct(f_unmap, f),
        }
        results[source] = db_results

        print("Total sites: %d" % n)
        print("  Valid:   %d (%0.1f)" % (n_val, pct(n_val, n)))
        print("  Invalid: %d (%0.1f)" % (n_inv, pct(n_inv, n)))
        print("  Mapped:  %d (%0.1f)" % (n_map, pct(n_map, n)))
        print("%% Mapped:  %0.1f\n" % pct(n_map, n_inv))
        print("Total site occurrences: %d" % f)
        print("  Valid:   %d (%0.1f)" % (f_val, pct(f_val, f)))
        print("  Invalid: %d (%0.1f)" % (f_inv, pct(f_inv, f)))
        print("  Mapped:  %d (%0.1f)" % (f_map, pct(f_map, f)))
        print("Pct occurrences mapped: %0.1f\n" % pct(f_map, f_inv))
    # Sample 100 invalid-unmapped (by unique sites)
    # Sample 100 invalid-mapped (by unique sites)
    results_df = pd.DataFrame.from_dict(results, orient='index')
    return results_df


# -- PLOTTING ----------------------------------------------------------

def plot_pc_pe_mods(all_mods):
    dbs = Counter([row.source for row in all_mods])
    dbs = sorted([(k, v) for k, v in dbs.items()], key=lambda x: x[1],
                 reverse=True)
    plt.ion()
    width = 0.8
    ind = np.arange(len(dbs)) + (width / 2.)
    plt.figure(figsize=(2, 2), dpi=150)
    for db_ix, (db, db_freq) in enumerate(dbs):
        db_mods = [row for row in all_mods if row.source == db]
        valid = [row for row in db_mods if row.valid == 1]
        invalid = [row for row in db_mods if row.valid == 0]
        h_valid = plt.bar(db_ix, len(valid), width=width, color='g')
        h_invalid = plt.bar(db_ix, len(invalid), width=0.8, color='r',
                            bottom=len(valid))
    plt.xticks(ind, [db[0] for db in dbs])
    ax = plt.gca()
    pf.format_axis(ax)
    plt.show()



def make_bar_plot(site_info, num_genes=120):
    # Build a dict based on gene name
    # Get counts summed across gene names
    gene_counts = defaultdict(lambda: 0)
    site_counts = defaultdict(list)
    for site, freq in site_info:
        if not (site.mapped_res and site.mapped_pos):
            pass
            #continue
        gene_counts[site.gene_name] += freq
        site_counts[site.gene_name].append((freq, site.valid,
                                            site.mapped_res, site.mapped_pos,
                                            site.description))
    # Sort the individual site counts by frequency
    for gene, freq_list in site_counts.items():
        site_counts[gene] = sorted(freq_list, key=lambda x: x[0], reverse=True)
    gene_counts = sorted([(k, v) for k, v in gene_counts.items()],
                          key=lambda x: x[1], reverse=True)

    def plot_sites(gene_count_subset, figsize, subplot_params, do_legend=True):
        ind = np.array(range(len(gene_count_subset)))
        plt.figure(figsize=figsize, dpi=150)
        width = 0.8
        handle_dict = {}
        for ix, (gene, freq) in enumerate(gene_count_subset):
            # Plot the stacked bars
            bottom = 0
            for site_freq, valid, mapped_res, mapped_pos, explanation \
                    in site_counts[gene]:
                # Don't show sites that are valid
                mapped = True if mapped_res and mapped_pos else False
                if valid:
                    color = 'gray'
                    handle_key = 'Valid'
                elif mapped and \
                        explanation.startswith('INFERRED_METHIONINE_CLEAVAGE'):
                    color = 'b'
                    handle_key = 'Methionine'
                elif mapped and \
                        explanation.startswith('INFERRED_MOUSE_SITE'):
                    color = 'c'
                    handle_key = 'Mouse'
                elif mapped and \
                        explanation.startswith('INFERRED_RAT_SITE'):
                    color = 'purple'
                    handle_key = 'Rat'
                elif mapped and \
                        explanation.startswith('INFERRED_ALTERNATIVE_ISOFORM'):
                    color = 'orange'
                    handle_key = 'Alternative isoform'
                elif mapped:
                    color = 'g'
                    handle_key = 'Manually mapped'
                elif not mapped and explanation is not None \
                                            and explanation != 'VALID':
                    color = 'r'
                    handle_key = 'Curated as incorrect'
                elif not valid:
                    assert False
                else:
                    assert False # Make sure we handled all cases above
                handle_dict[handle_key] = \
                        plt.bar(ix + 0.4, site_freq, bottom=bottom,
                                color=color, linewidth=0.5, width=width)
                bottom += site_freq
        plt.xticks(ind + (width / 2.), [x[0] for x in gene_count_subset],
                   rotation='vertical')
        plt.ylabel('Stmts with invalid sites')
        plt.xlim((0, max(ind)+1))
        ax = plt.gca()
        pf.format_axis(ax)
        plt.subplots_adjust(**subplot_params)
        if do_legend:
            plt.legend(loc='upper right', handles=list(handle_dict.values()),
                       labels=list(handle_dict.keys()), fontsize=pf.fontsize,
                       frameon=False)
        plt.show()
    plot_sites(gene_counts[0:4], (0.23, 2),
               {'left': 0.24, 'right': 0.52, 'bottom': 0.31}, do_legend=False)
    plot_sites(gene_counts[4:num_genes], (11, 2),
               {'bottom': 0.31, 'left': 0.06, 'right':0.96})
    plot_sites(gene_counts, (11, 2),
               {'left': 0.24, 'right': 0.52, 'bottom': 0.31}, do_legend=False)
    return gene_counts

# ---------------------


def plot_site_count_dist(sites, num_sites=240):
    # Plot site frequencies, colored by validity
    sites.sort(key=lambda s: s[1], reverse=True)
    width = 0.8
    plt.figure(figsize=(11, 2), dpi=150)
    ind = np.arange(num_sites) + (width / 2.)
    for site_ix, site in enumerate(sites[:num_sites]):
        if site.valid:
            color = 'g'
        else:
            color = 'r'
        plt.bar(site_ix, site.freq, color=color)
    ax = plt.gca()
    pf.format_axis(ax)
    plt.show()


def map_pc_sites():
    pm = ProtMapper(use_cache=True, cache_path=CACHE_PATH)
    agent_files = {
        'hprd': 'output/biopax/PathwayCommons10.hprd.BIOPAX.pkl',
        'kegg': 'output/biopax/PathwayCommons10.kegg.BIOPAX.pkl',
        'panther': 'output/biopax/PathwayCommons10.panther.BIOPAX.pkl',
        'pid': 'output/biopax/PathwayCommons10.pid.BIOPAX.pkl',
        'psp_pc': 'output/biopax/PathwayCommons10.psp.BIOPAX.pkl',
        'pc_tsv': 'output/psp_kinase_substrate_tsv.pkl',
        'reactome': 'output/biopax/PathwayCommons10.reactome.BIOPAX.pkl',
        'wp': 'output/biopax/PathwayCommons10.wp.BIOPAX.pkl',
        'psp_biopax': 'output/biopax/Kinase_substrates.pkl',
        #'reactome_human': 'output/biopax/Homo_sapiens.pkl',
    }
    all_sites = {}
    for db_name, agent_file in agent_files.items():
        sites = map_agents(agent_file, pm, db_name)
        all_sites[db_name] = sites
    with open(BIOPAX_SITES_BY_DB, 'wb') as f:
        pickle.dump(all_sites, f)
    pm.save_cache()


def map_bel_sites():
    with open(BEL_AGENTS, 'rb') as f:
        bel_agents = pickle.load(f)
    pm = ProtMapper(use_cache=True, cache_path=CACHE_PATH)
    bel_sites = map_agents(BEL_AGENTS, pm, 'bel')
    with open(BEL_SITES, 'wb') as f:
        pickle.dump(bel_sites, f)
    pm.save_cache()


def create_site_csv():
    all_sites = []
    # Load Biopax sites
    with open(BIOPAX_SITES_BY_DB, 'rb') as f:
        pc_sites = pickle.load(f)
    for db, sites in pc_sites.items():
        for ms, freq in sites:
            all_sites.append(ms_to_si(db, freq, ms))
    # Load BEL sites
    with open(BEL_SITES, 'rb') as f:
        bel_sites = pickle.load(f)
    for ms, freq in bel_sites:
        all_sites.append(ms_to_si('bel', freq, ms))
    # Load sites from reading
    with open(READER_SITES, 'rb') as f:
        sites_by_reader = pickle.load(f)
    for reader, sites in sites_by_reader.items():
        for ms, freq in sites.items():
            all_sites.append(ms_to_si(reader, freq, ms))
    header = [[field.upper() for field in all_sites[0]._asdict().keys()]]
    rows = header + replace_nones(all_sites)
    write_unicode_csv(ALL_SITES_CSV, rows)


def plot_site_stats(csv_file):
    site_df = pd.read_csv(csv_file)
    # Drop rows with error_code (invalid gene names in BEL)
    site_df = site_df[site_df.ERROR_CODE.isna()]
    results = print_stats(site_df)
    # Now make figures for the sites
    by_site = results[['Valid Sites', 'Mapped Sites', 'Unmapped Sites']]
    by_site_pct = results[['Valid Sites Pct.', 'Mapped Sites Pct. Total',
                           'Unmapped Sites Pct. Total']]
    by_occ = results[['Valid Occ.', 'Mapped Occ.', 'Unmapped Occ.']]
    by_occ_pct = results[['Valid Occ. Pct.', 'Mapped Occ. Pct. Total',
                          'Unmapped Occ. Pct. Total']]
    for df, kind in ((by_site, 'by_site'), (by_site_pct, 'by_site_pct'),
                     (by_occ, 'by_occ'), (by_occ_pct, 'by_occ_pct')):
        plt.figure()
        df.plot(kind='bar', stacked=True)
        plt.subplots_adjust(bottom=0.2)
        plt.savefig('plots/site_stats_%s.pdf' % kind)


if __name__ == '__main__':
    # Map sites from Pathway Commons
    if sys.argv[1] == 'map_pc_sites':
        map_pc_sites()
    # Map sites from BEL large corpus
    elif sys.argv[1] == 'map_bel_sites':
        map_bel_sites()
    # Create a single CSV file containing information about all sites from
    # databases
    elif sys.argv[1] == 'create_site_csv':
        create_site_csv()
    # Load the CSV file and plot site statistics
    elif sys.argv[1] == 'plot_site_stats':
        plot_site_stats(ALL_SITES_CSV)
    else:
        pass
