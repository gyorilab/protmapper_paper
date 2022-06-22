import csv
import sys
import pickle
from collections import Counter, defaultdict
import matplotlib
matplotlib.use('agg')
import pandas as pd
from matplotlib import pyplot as plt
from indra.util import plot_formatting as pf
from indra.ontology.bio import bio_ontology


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


def create_export(site_stmts, mapping_results, export_file, evs_file):
    from indra.statements import Agent
    from indra.databases import uniprot_client, hgnc_client

    # Make header for main export file
    export_header = ['ID',
                     'CTRL_NS', 'CTRL_ID', 'CTRL_GENE_NAME', 'CTRL_IS_KINASE',
                     'TARGET_UP_ID', 'TARGET_GENE_NAME', 'TARGET_RES',
                     'TARGET_POS','SOURCES'
                     ]
    # Make header for evidence export file
    evidence_header = ['ID', 'SOURCE', 'PMID', 'DBID', 'TEXT',
                       'DESCRIPTION',
                       'ORIG_UP_ID', 'ORIG_RES', 'ORIG_POS',
                       'MAPPED_UP_ID', 'MAPPED_RES', 'MAPPED_POS']

    site_info = {}
    site_evidence = defaultdict(list)
    idx = 0
    # site_stmts is a dict with structure like
    # site_stmts[('Q15438', 'T', '395')]['rhs']['signor'] ->
    # [Phosphorylation(PRKCD(), CYTH1(), T, 395)]
    for (orig_up_id, orig_res, orig_pos), stmt_dict in site_stmts.items():
        # We skip sites that are missing residue or position
        if not orig_res or not orig_pos:
            continue
        # Next, we construct keys for the *final* site (either valid to begin
        # with or mapped to be valid), and if there is no valid final
        # site, we skip the site
        ms = mapping_results[(orig_up_id, orig_res, orig_pos)]
        if ms.valid:
            final_site = [ms.up_id, ms.orig_res, ms.orig_pos]
        elif ms.mapped_res and ms.mapped_pos:
            final_site = [ms.mapped_id, ms.mapped_res, ms.mapped_pos]
        else:
            continue
        # Skip non-human substrates
        if not uniprot_client.is_human(final_site[0]):
            continue
        target_gene_name = uniprot_client.get_gene_name(final_site[0])
        final_site = [final_site[0], target_gene_name, final_site[1],
                      final_site[2]]

        # We now look at all the Statements where the given site
        # appears as a substrate and get controllers and evidences
        for source, stmts in stmt_dict['rhs'].items():
            for stmt in stmts:
                # If there is no controller, we skip the entry
                if stmt.enz is None:
                    continue
                # We next get the grounding for the controller and
                # if there is no grounding, we skip it
                ctrl_ns, ctrl_id = stmt.enz.get_grounding()
                if ctrl_ns not in ['UP', 'HGNC', 'FPLX'] or ctrl_id is None:
                    continue

                ctrl_gene_name = None
                ctrl_is_kinase = False
                # Get human gene name for UniProt entries
                if ctrl_ns == 'UP':
                    # Skip non-human protein controllers
                    if not uniprot_client.is_human(ctrl_id):
                        continue
                    ctrl_gene_name = uniprot_client.get_gene_name(ctrl_id)
                    if hgnc_client.is_kinase(ctrl_gene_name):
                        ctrl_is_kinase = True
                # Map human gene names to UniProt IDs
                elif ctrl_ns == 'HGNC':
                    gene_name = hgnc_client.get_hgnc_name(ctrl_id)
                    if hgnc_client.is_kinase(gene_name):
                        ctrl_is_kinase = True
                    up_id = hgnc_client.get_uniprot_id(ctrl_id)
                    if up_id:
                        ctrl_ns = 'UP'
                        ctrl_gene_name = gene_name
                        ctrl_id = up_id
                elif ctrl_ns == 'FPLX':
                    children = bio_ontology.get_children('FPLX', ctrl_id,
                                                         ns_filter={'HGNC'})
                    for _, hgnc_id in children:
                        gene_name = hgnc_client.get_hgnc_name(hgnc_id)
                        if hgnc_client.is_kinase(gene_name):
                            ctrl_is_kinase = True
                            break

                # We can now make a full key that contains the controller
                # as well as the target and final site
                final_annot_key = tuple([ctrl_ns, ctrl_id, ctrl_gene_name,
                                         ctrl_is_kinase] + final_site)
                # We use this full key to store evidences and mapping details
                if final_annot_key not in site_info:
                    site_info[final_annot_key] = idx
                    idx += 1
                # Note: we do get multiple pieces of evidence, e.g.,
                # from biopax
                for ev in stmt.evidence:
                    site_evidence[final_annot_key].append([ev, source, ms])

    # Now make the actual export tables
    def sanitize_ev_text(txt):
        if txt is None:
            return ''
        else:
            txt = txt.replace('\n', ' ')
            return txt

    export_rows = [export_header]
    evidence_rows = [evidence_header]
    for key, idx in site_info.items():
        (ctrl_ns, ctrl_id, ctrl_gene_name, ctrl_is_kinase, target_up_id,
            target_gene_name, target_res, target_pos) = key
        export_row = [str(idx),
                      ctrl_ns, ctrl_id, ctrl_gene_name, ctrl_is_kinase,
                      target_up_id, target_gene_name, target_res, target_pos]
        # Now get evidences
        evs = site_evidence[key]
        sources = sorted(list({s for e, s, m in evs}))
        export_row.append(','.join(sources))
        export_rows.append(export_row)
        for evidence, source, ms in evs:
            if source == 'bel':
                source_id = evidence.source_id[:16]
            else:
                source_id = evidence.source_id
            row = [str(idx), source, evidence.pmid,
                   source_id, sanitize_ev_text(evidence.text),
                   ms.description,
                   ms.up_id, ms.orig_res, ms.orig_pos,
                   ms.mapped_id, ms.mapped_res, ms.mapped_pos]
            evidence_rows.append(row)

    with open(export_file, 'wt') as fh:
        csvwriter = csv.writer(fh)
        csvwriter.writerows(export_rows)
    with open(evs_file, 'wt') as fh:
        csvwriter = csv.writer(fh)
        csvwriter.writerows(evidence_rows)


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

    # Print data for Figure 1B on MAPK1 sites
    stats = defaultdict(list)
    sites = [('T', 182), ('T', 183), ('T', 185), ('Y', 184), ('Y', 185),
             ('Y', 187)]
    for _, row in site_df[site_df['GENE_NAME'] == 'MAPK1'].iterrows():
        pos = int(row['ORIG_POS'])
        stats[(row['ORIG_RES'], pos)].append((row['SOURCE'], row['FREQ']))
    print('Table for Figure 1B\n------------')
    print('Site', 'Occ.', 'Sources')
    for res, pos in sites:
        src = len({s for s, c in stats[(res, pos)]})
        occ = sum([c for s, c in stats[(res, pos)]])
        print('%s%s' % (res, pos), occ, src)

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
    site_df = pd.read_csv(csv_file)
    # Drop rows with error_code (invalid gene names in BEL)
    site_df = site_df[site_df.ERROR_CODE.isna()]
    #sources = site_df.SOURCE.unique()
    sources = ['psp', 'signor', 'hprd', 'pid', 'reactome', 'bel',
               'reach', 'sparser', 'rlimsp']
    source_labels = ['PSP', 'SIGNOR', 'HPRD', 'NCI-PID',
                     'Reactome', 'BEL', 'Reach', 'Sparser', 'RLIMS-P']

    # First, generate plot of absolute valid/invalid sites
    # (Figure 1B and Figure 2E)
    fig_valid = plt.figure(figsize=(1.8, 2.2), dpi=150)
    fig_mapped = plt.figure(figsize=(1.8, 2.2), dpi=150)
    for ix, source in enumerate(sources):
        source_anns = site_df[site_df.SOURCE == source]
        valid = source_anns[source_anns.VALID == True]
        invalid = source_anns[source_anns.VALID != True]
        unmapped = invalid[invalid.MAPPED_RES.isna() |
                           invalid.MAPPED_POS.isna()]
        mapped = source_anns[~source_anns.MAPPED_RES.isna() &
                             ~source_anns.MAPPED_POS.isna()]
        print("%s, %d, %d, %d, %d" % (source, len(valid), len(invalid),
                                      len(unmapped), len(mapped)))
        # Plot for valid stats
        fig_valid.gca().bar(ix, height=len(valid), color='blue',
                            label='In Ref Seq')
        fig_valid.gca().bar(ix, height=len(invalid), bottom=len(valid),
                            color='red', label='Not in Ref Seq')
        # Plot for mapped stats
        fig_mapped.gca().bar(ix, height=len(valid), color='blue',
                             label='In Ref Seq')
        fig_mapped.gca().bar(ix, height=len(mapped), bottom=len(valid),
                             color='green')
        fig_mapped.gca().bar(ix, height=len(unmapped),
                             bottom=(len(valid) + len(mapped)), color='red')
    for fig in [fig_valid, fig_mapped]:
        ax = fig.gca()
        pf.format_axis(ax)
        ax.xaxis.set_ticks(range(len(sources)))
        ax.set_xticklabels(labels=source_labels, rotation='vertical')
        ax.set_ylabel('Unique Site Annotations')
        fig.subplots_adjust(left=0.31, bottom=0.25, top=0.95)
        #fig.legend(loc='upper right', fontsize=pf.fontsize)
    fig_valid.savefig('%s_valid_counts.pdf' % output_base)
    fig_mapped.savefig('%s_mapped_counts.pdf' % output_base)

    # Now generate the plot of invalid site percentages
    # (Figure 1C and Figure 2F)
    results = []
    fig_valid = plt.figure(figsize=(1.8, 2.4), dpi=150)
    fig_mapped = plt.figure(figsize=(1.8, 2.4), dpi=150)
    for ix, source in enumerate(sources):
        source_anns = site_df[site_df.SOURCE == source]
        valid = source_anns[source_anns.VALID == True]
        invalid = source_anns[source_anns.VALID != True]
        unmapped = invalid[invalid.MAPPED_RES.isna() |
                           invalid.MAPPED_POS.isna()]
        mapped = source_anns[~source_anns.MAPPED_RES.isna() &
                             ~source_anns.MAPPED_POS.isna()]
        pct_mapped = 100 * len(mapped) / len(source_anns)
        pct_unmapped = 100 * len(unmapped) / len(source_anns)
        pct_valid = 100 * len(valid) / len(source_anns)
        pct_invalid = 100 * len(invalid) / len(source_anns)
        results.append([
            source, len(source_anns), len(valid),
            round(100 * len(valid) / len(source_anns), 1),
            len(invalid),
            round(100 * len(invalid) / len(source_anns), 1),
            len(mapped),
            round(100 * len(mapped) / len(invalid), 1)])
        #print("%s, %d, %d, %d, %d" % (source, len(valid), len(invalid),
        #                              len(unmapped), len(mapped)))
        fig_valid.gca().bar(ix, height=pct_valid, color='blue')
        fig_valid.gca().bar(ix, height=pct_invalid, bottom=pct_valid,
                            color='red')
        fig_mapped.gca().bar(ix, height=pct_valid, color='blue')
        fig_mapped.gca().bar(ix, height=pct_mapped,
                             bottom=pct_valid, color='green')
        fig_mapped.gca().bar(ix, height=pct_unmapped,
                             bottom=(pct_valid + pct_mapped), color='red')
    for fig in [fig_valid, fig_mapped]:
        ax = fig.gca()
        pf.format_axis(ax)
        ax.xaxis.set_ticks(range(len(sources)))
        ax.set_xticklabels(labels=source_labels, rotation='vertical')
        ax.set_ylabel('Pct. Unique Site Annotations')
        fig.subplots_adjust(left=0.31, bottom=0.25, top=0.90)
    fig_valid.savefig('%s_valid_pcts.pdf' % output_base)
    fig_mapped.savefig('%s_mapped_pcts.pdf' % output_base)

    # This is where Table 1 comes from
    result_df = pd.DataFrame.from_records(results, columns=
            ['SOURCE', 'TOTAL', 'VALID', 'VALID_PCT', 'INVALID', 'INVALID_PCT',
             'MAPPED', 'MAPPED_PCT_INVALID'])
    result_df.to_csv('%s_result_df.csv' % output_base)
    return result_df


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
    elif sys.argv[1] == 'plot_annot_stats':
        input_file = sys.argv[2]
        output_base = sys.argv[3]
        result_df = plot_annot_stats(input_file, output_base)
    elif sys.argv[1] == 'export':
        site_pkl_file = sys.argv[2]
        mapping_results_file = sys.argv[3]
        export_file = sys.argv[4]
        evs_file = sys.argv[5]
        with open(site_pkl_file, 'rb') as fh:
            site_stmts = pickle.load(fh)
        with open(mapping_results_file, 'rb') as f:
            mapping_results = pickle.load(f)
        create_export(site_stmts, mapping_results, export_file, evs_file)
