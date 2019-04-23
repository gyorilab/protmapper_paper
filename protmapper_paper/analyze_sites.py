import csv
import sys
import pickle
from collections import Counter
import matplotlib
matplotlib.use('agg')
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from indra.util import plot_formatting as pf


def create_site_csv(site_dict, mapping_results, site_file, annot_file):
    site_header = ['SOURCE', 'GENE_NAME', 'UP_ID', 'ERROR_CODE', 'VALID',
              'ORIG_RES', 'ORIG_POS', 'MAPPED_ID', 'MAPPED_RES', 'MAPPED_POS',
              'DESCRIPTION', 'SIDE', 'HAS_SUBJECT', 'FREQ', 'TOTAL_CONTROLLERS',
              'PROTEIN_CONTROLLERS']
    ann_header = ['SOURCE', 'GENE_NAME', 'UP_ID', 'ERROR_CODE', 'VALID',
              'ORIG_RES', 'ORIG_POS', 'MAPPED_ID', 'MAPPED_RES', 'MAPPED_POS',
              'DESCRIPTION', 'SIDE', 'CTRL_NAME', 'CTRL_NS', 'CTRL_ID',
              'CTRL_IS_PROTEIN', 'CTRL_FREQ']
    all_sites = [site_header]
    annotations = [ann_header]
    for site in site_dict:
        ms = mapping_results[site]
        up_id, res, pos = site
        for side in ('lhs', 'rhs'):
            for source, stmts in site_dict[site][side].items():
                if len(stmts) == 0:
                    continue
                elif side == 'rhs':
                    none_enz = [s for s in stmts if s.agent_list()[0] is None]
                    with_enz = [s for s in stmts
                                if s.agent_list()[0] is not None]
                    # Count distinct controllers and protein controllers
                    total_controllers = []
                    protein_controllers = set()
                    for s in with_enz:
                        ctrl = s.agent_list()[0]
                        db_ns, db_id = ctrl.get_grounding()
                        if db_ns is None or db_id is None:
                            continue
                        total_controllers.append((ctrl.name, db_ns, db_id))
                        if db_ns in ('FPLX', 'HGNC', 'UP'):
                            protein_controllers.add((ctrl.name, db_ns, db_id))
                    ctrl_ctr = Counter(total_controllers)
                    for (cname, cns, cid), cfreq in ctrl_ctr.items():
                        is_protein = True if cns in ('FPLX', 'HGNC', 'UP') \
                                          else False
                        ann = [
                           source, ms.gene_name, ms.up_id, ms.error_code,
                           ms.valid, ms.orig_res, ms.orig_pos, ms.mapped_id,
                           ms.mapped_res, ms.mapped_pos, ms.description,
                           side, cname, cns, cid, is_protein, cfreq]
                        annotations.append(ann)
                    # Add count for stmts without subject
                    if len(none_enz) > 0:
                        all_sites.append([
                               source, ms.gene_name, ms.up_id, ms.error_code,
                               ms.valid, ms.orig_res, ms.orig_pos, ms.mapped_id,
                               ms.mapped_res, ms.mapped_pos, ms.description,
                               side, False, len(none_enz), 0, 0])
                    # Add count for stmts *with* subject
                    if len(with_enz) > 0:
                        all_sites.append([
                               source, ms.gene_name, ms.up_id, ms.error_code,
                               ms.valid, ms.orig_res, ms.orig_pos, ms.mapped_id,
                               ms.mapped_res, ms.mapped_pos, ms.description,
                               side, True, len(with_enz),
                               len(ctrl_ctr), len(protein_controllers)])
                else:
                    all_sites.append([
                           source, ms.gene_name, ms.up_id, ms.error_code,
                           ms.valid, ms.orig_res, ms.orig_pos, ms.mapped_id,
                           ms.mapped_res, ms.mapped_pos, ms.description, side,
                           True, len(stmts)])
    print("Saving %d site entries to %s" % (len(all_sites)-1, site_file))
    with open(site_file, 'wt') as f:
        csvwriter = csv.writer(f)
        csvwriter.writerows(all_sites)
    print("Saving %d annotation entries to %s" %
                        (len(annotations)-1, annot_file))
    with open(annot_file, 'wt') as f:
        csvwriter = csv.writer(f)
        csvwriter.writerows(annotations)


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


def plot_site_stats(csv_file, output_base):
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
    by_occ_corr = results[['Valid Occ.', 'Invalid Occ.']]
    by_occ_corr_pct = results[['Valid Occ. Pct.', 'Invalid Occ. Pct.']]

    for df, kind in ((by_site, 'by_site'), (by_site_pct, 'by_site_pct'),
                     (by_occ, 'by_occ'), (by_occ_pct, 'by_occ_pct')):
        plt.figure()
        df.plot(kind='bar', stacked=True, color=['blue', 'green', 'orange'])
        plt.subplots_adjust(bottom=0.2)
        plt.savefig('%s_%s.pdf' % (output_base, kind))
    for df, kind in ((by_occ_corr, 'by_occ_corr'),
                     (by_occ_corr_pct, 'by_occ_corr_pct')):
        plt.figure()
        df.plot(kind='bar', stacked=True, color=['blue', 'orange'])
        plt.subplots_adjust(bottom=0.2)
        plt.savefig('%s_%s.pdf' % (output_base, kind))


def plot_annot_stats(csv_file, output_base):
    site_df = pd.read_csv(csv_file, output_base)
    # Drop rows with error_code (invalid gene names in BEL)
    site_df = site_df[site_df.ERROR_CODE.isna()]
    


def get_sites_by_source(sites_dict, source, side):
    filt_dict = {}
    for site in sites_dict.keys():
        source_stmts = sites_dict[site][side].get(source)
        if source_stmts:
            filt_dict[site] = source_stmts
    return filt_dict


def site_sample(all_sites_file, output_file):
    all_sites = pd.read_csv(all_sites_file, dtype={'ORIG_POS': 'str'})
    return all_sites

if __name__ == '__main__':
    # Create a single CSV file containing information about all sites from
    # databases
    if sys.argv[1] == 'create_site_csv':
        site_dict_file = sys.argv[2]
        mapping_results_file = sys.argv[3]
        csv_file = sys.argv[4]
        annot_file = sys.argv[5]
        with open(site_dict_file, 'rb') as f:
            site_dict = pickle.load(f)
        with open(mapping_results_file, 'rb') as f:
            mapping_results = pickle.load(f)
        create_site_csv(site_dict, mapping_results, csv_file, annot_file)
    # Load the CSV file and plot site statistics
    elif sys.argv[1] == 'plot_site_stats':
        input_file = sys.argv[2]
        output_base = sys.argv[3]
        plot_site_stats(input_file, output_base)
    elif sys.argv[1] == 'site_samples':
        all_sites_file = sys.argv[2]
        output_file = sys.argv[3]
        sample = site_sample(all_sites_file, output_file)

