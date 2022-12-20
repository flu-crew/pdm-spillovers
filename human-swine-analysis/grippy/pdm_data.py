# -*- coding: utf-8 -*-
import datetime
import re
from dendropy import Tree, Node
from Bio import SeqIO

from grippy import helpers, IAVSequence, pdm09
from grippy.helpers import get_strain_name_list


def rename_vdl2_seq(iav: IAVSequence, redo=False):
    state = iav.state
    if not redo and (iav.usda_id and state):
        id = iav.usda_id
        iav.orig_name = iav.name
    else:
        if not state:
            state = 'USA'
        if redo:
            id = iav.orig_name[8:].lstrip('0')
        else:
            id = iav.name[8:].lstrip('0')
            iav.orig_name = iav.name
    name = 'A/swine/%s/%s/%d' % (state, id, iav.date.year)
    iav.name = name
    if not redo and (iav.usda_id and iav.state):
        iav.full_name = '%s|A|H1|%s|Swine|1A.3.3.2|%s' % (name, state, iav.date.strftime('%Y-%m-%d'))
    else:
        iav.full_name = 'ISUVDL|%s|H1|%s|%s' % (name, state, iav.date.strftime('%Y-%m-%d'))


def extract_active_barcode(id: str) -> str:
    if id.startswith('A0'):
        return id
    else:
        id = id.lstrip('0')
        return re.search(r'([0-9]+)([a-zA-Z]+[0-9]*)?', id).group(1)


def filter_IRD_data(orig_path: str, pdm_list: str, new_path: str, min_seq_len=1690, start_date='2009-01-1',
                    end_date='2021-10-31'):
    # Process swine sequences:
    iav_sequences = helpers.read_iav_sequences(orig_path, host='swine')
    ird_pdm_list = get_strain_name_list(pdm_list)
    # print(ird_pdm_list)
    filtered_seqs = []

    start_date = datetime.datetime.strptime(start_date, '%Y-%m-%d')
    end_date = datetime.datetime.strptime(end_date, '%Y-%m-%d')
    for iav in iav_sequences:
        if iav.full_name.replace('|', '_') in ird_pdm_list and start_date <= iav.date <= end_date and\
                len(iav.seq) >= min_seq_len:
            filtered_seqs.append(iav)
            print(f'Saved {iav.full_name}')

    print(f'Total saved {len(filtered_seqs)}')
    helpers.write_iav_sequences(new_path, filtered_seqs)


def filter_gisaid_data(orig_path: str, pdm_list: str, new_path: str, min_seq_len=1690, start_date='2009-01-1',
                    end_date='2021-05-31'):
    # Process swine sequences:
    iav_sequences = helpers.read_iav_sequences(orig_path, host='human')
    gisaid_pdm_list = get_strain_name_list(pdm_list)
    # print(gisaid_pdm_list)
    filtered_seqs = []

    start_date = datetime.datetime.strptime(start_date, '%Y-%m-%d')
    end_date = datetime.datetime.strptime(end_date, '%Y-%m-%d')
    for iav in iav_sequences:
        if iav.full_name.replace('|', '_') in gisaid_pdm_list and start_date <= iav.date <= end_date and\
                len(iav.seq) >= min_seq_len:
            filtered_seqs.append(iav)
            print(f'Saved {iav.full_name}')

    print(f'Total saved {len(filtered_seqs)}')
    helpers.write_iav_sequences(new_path, filtered_seqs)


def merge_swine_data(ird_path, vdl1_path, vdl2_path, active_path, vdl2_pdm_list, active_pdm_list,
                     output_path):
    ird_iavs = helpers.read_iav_sequences(ird_path, protein='HA', host='swine')
    vdl1_iavs = helpers.read_iav_sequences(vdl1_path, protein='HA', host='swine')
    vdl2_iavs = helpers.read_iav_sequences(vdl2_path, protein='HA', host='swine')
    active_iavs = helpers.read_iav_sequences(active_path, protein='HA', host='swine')
    vdl2_pdm_names = get_strain_name_list(vdl2_pdm_list)
    active_pdm_names = get_strain_name_list(active_pdm_list)
    vdl2_pdms = [iav for iav in vdl2_iavs if iav.full_name in vdl2_pdm_names and len(iav.seq) >= 1690]
    active_pdms = [iav for iav in active_iavs if iav.full_name in active_pdm_names and len(iav.seq) >= 1690]

    all_pdms = ird_iavs + vdl1_iavs

    initial_ids = [iav.get_swine_name_id() for iav in ird_iavs + vdl1_iavs]
    all_ids = set(initial_ids)
    original_vdl_pdm_ids = set()
    for vdl2_iav in vdl2_pdms:
        original_vdl_pdm_ids.add(vdl2_iav.get_swine_name_id())
        rename_vdl2_seq(vdl2_iav)
        if vdl2_iav.get_swine_name_id() not in all_ids:
            rename_vdl2_seq(vdl2_iav, redo=True)
            all_ids.add(vdl2_iav.get_swine_name_id())
            all_pdms.append(vdl2_iav)
            print(vdl2_iav.full_name)

    for active_iav in active_pdms:
        barcode = extract_active_barcode(active_iav.get_swine_name_id())
        if barcode not in all_ids:
            if barcode in original_vdl_pdm_ids:
                print('%s seen among vdl2' % active_iav.full_name)
                continue
            all_pdms.append(active_iav)
            all_ids.add(barcode)
            print(active_iav.full_name)

    helpers.write_iav_sequences(output_path, all_pdms)


def merge_vdl_with_existing():
    existing_sequences = helpers.read_iav_sequences('../../genes/Swine2010-2021_pdm09-HA_US_all.fasta', protein='HA',
                                               host='swine')
    existing_name_ids = [iav.get_swine_name_id() for iav in existing_sequences]
    ird_all_iavs = helpers.read_iav_sequences('../../genes/Swine2010-21_IRD_Sep07-2021_filtered.fasta')
    vdl_recent = helpers.read_iav_sequences('../../genes/ISUVDL_HA_May2020-June2021.fasta')
    pdm_list = get_strain_name_list('../../genes/ISUVDL_recent_HApdm.txt')
    vdl_pdm = [iav for iav in vdl_recent if iav.full_name in pdm_list and len(iav.seq) >= 1650]
    active_vdl = helpers.read_iav_sequences('../../genes/Active_H1_all.fasta')
    active_pdm_list = get_strain_name_list('../../genes/Active_HA_pdm09_2020-2021.txt')
    active_pdm = [iav for iav in active_vdl if iav.full_name in active_pdm_list and len(iav.seq) >= 1650]

    all_sequences = existing_sequences
    all_name_ids = existing_name_ids

    # Merge in all missing IRD pdm sequences (except for TOSU -- as it was previously filtered).
    for ird_iav in ird_all_iavs:
        if not ird_iav.in_header('TOSU') and ird_iav.get_swine_name_id() not in all_name_ids:
            all_sequences.append(ird_iav)
            all_name_ids.append(ird_iav.get_swine_name_id())
            print(ird_iav.full_name)

    # Merge in all VDL sequences as of June 2021.
    original_vdl_pdm_ids = set()
    for vdl2_iav in vdl_pdm:
        original_vdl_pdm_ids.add(vdl2_iav.get_swine_name_id())
        rename_vdl2_seq(vdl2_iav)
        if vdl2_iav.name == 'A/swine/Iowa/71212/2020':  # Skip this unusual pdm.
            continue
        if vdl2_iav.get_swine_name_id() not in existing_name_ids:
            all_sequences.append(vdl2_iav)
            all_name_ids.append(vdl2_iav.get_swine_name_id())
            print(vdl2_iav.full_name)

    # Merge in all DARPA sequences as of June 2021.
    for active_vdl in active_pdm:
        barcode = extract_active_barcode(active_vdl.get_swine_name_id())
        if barcode not in all_name_ids:
            if barcode in original_vdl_pdm_ids:
                print('%s seen among vdl2' % active_vdl.full_name)
                continue
            all_sequences.append(active_vdl)
            all_name_ids.append(barcode)
            print(active_vdl.full_name)

    helpers.write_iav_sequences('../../pdm09-data/swine/Swine2010-2021_pdm09-HA_US_June2021-merged2_all_lengths.fasta', all_sequences)

    # all_sequences = sorted(all_sequences, key=lambda x: (str(x.seq), x.date))
    # # merged_seqs = []
    # same_seqs = [all_sequences[0]]
    # for i in range(1, len(all_sequences) + 1):
    #     if i < len(all_sequences) and str(all_sequences[i].seq) == str(same_seqs[0].seq)\
    #             and all_sequences[i].date == same_seqs[0].date:
    #         same_seqs.append(all_sequences[i])
    #     else:
    #         if len(same_seqs) > 1:
    #             seqs_by_type = {}
    #             for seq in same_seqs:
    #                 print(seq.full_name)
    #                 if seq.full_name.count('ISU/') > 0:
    #                     if seq.usda_id:
    #                         seqs_by_type['vdl2_a0'] = seq
    #                     else:
    #                         seqs_by_type['vdl2'] = seq
    #                 elif seq.contains_token('VDL'):
    #                     seqs_by_type['vdl'] = seq
    #                 else:
    #                     seqs_by_type['fludb'] = seq
    #             if seqs_by_type.get('fludb'):
    #                 merged_seqs.append(seqs_by_type['fludb'])
    #             elif seqs_by_type.get('vdl2_a0'):
    #                 merged_seqs.append(seqs_by_type['vdl2_a0'])
    #             elif seqs_by_type.get('vdl'):
    #                 merged_seqs.append(seqs_by_type['vdl'])
    #             else:
    #                 merged_seqs.append(seqs_by_type['vdl2'])
    #         else:
    #             merged_seqs.append(same_seqs[0])
    #         if i < len(all_sequences):
    #             same_seqs = [all_sequences[i]]
    # print('----------')
    # all_usda_ids = [iav.usda_id for iav in merged_seqs if iav.usda_id]
    # for usda_id in set(all_usda_ids):
    #     if all_usda_ids.count(usda_id) > 1:
    #         for iav in [iav for iav in merged_seqs if iav.usda_id == usda_id]:
    #             print(iav.full_name)


def merge_NA_pdms():
    zeb_seqs = helpers.read_iav_sequences('../../genes/DARPA/Zeb-NA-USA-pdms.fasta')
    gisaid_seqs = helpers.read_iav_sequences('../../genes/Human2019-2021-Aug2021-NA-pdm-US-gisaid-filtered.fasta')
    zeb_names = {seq.name for seq in zeb_seqs}
    all_seqs = zeb_seqs
    for gisaid_seq in gisaid_seqs:
        if gisaid_seq.name not in zeb_names:
            # print(gisaid_seq.full_name)
            all_seqs.append(gisaid_seq)

    vdl_recent = helpers.read_iav_sequences('../../genes/ISUVDL_NA_May2020-June2021.fasta')
    pdm_list = get_strain_name_list('../../genes/ISUVDL_June2021_pdm_NAs.txt')
    vdl_pdm = [iav for iav in vdl_recent if iav.full_name in pdm_list]
    active_vdl = helpers.read_iav_sequences('../../genes/DARPA/Active_N1_all.fasta')
    active_pdm_list = get_strain_name_list('../../genes/Active_N1_pdms.txt')
    active_pdm = [iav for iav in active_vdl if iav.full_name in active_pdm_list]

    all_seq_ids = [seq.get_swine_name_id() for seq in all_seqs if seq.full_name.lower().count('swine') > 0]
    for vdl2_iav in vdl_pdm:
        rename_vdl2_seq(vdl2_iav)
        if vdl2_iav.get_swine_name_id() not in all_seq_ids:
            all_seqs.append(vdl2_iav)
            all_seq_ids.append(vdl2_iav.get_swine_name_id())
            print(vdl2_iav.full_name)

    for active_vdl in active_pdm:
        if extract_active_barcode(active_vdl.get_swine_name_id()) not in all_seq_ids:
            all_seqs.append(active_vdl)
            all_seq_ids.append(active_vdl.get_swine_name_id())
            print(active_vdl.full_name)
            active_vdl.full_name = 'ISUVDL|%s|USA|%s' % (active_vdl.name, active_vdl.date.strftime('%Y-%m-%d'))

    helpers.write_iav_sequences('../../genes/human-swine-NA-Aug2021-pdms.fasta', all_seqs)


def finalize_NAs():
    na_cds_seqs = helpers.read_iav_sequences('../../genes/human-swine-NA-Aug2021-pdms-cds.aln')
    added_names = set()
    filtered = []
    for seq in na_cds_seqs:
        if seq.name not in added_names:
            filtered.append(seq)
            added_names.add(seq.name)
            start = ''
            start_tokens = ['ISUVDL', 'wgs']
            for token in start_tokens:
                if seq.contains_token(token):
                    start = token + '|'
                    break
            seq.full_name = start + '%s|%s|USA|%s' % (seq.name, seq.subtype if seq.subtype != 'unknown' else '',
                                                      seq.date.strftime('%Y-%m-%d'))
    helpers.write_iav_sequences('../../genes/human-swine-NA-Aug2021-pdms-cds-filtered.aln', filtered)


def finalize_HAs():
    ha_cds_seqs = helpers.read_iav_sequences('../../genes/human-swine-Aug2021-cds.aln')
    added_names = set()
    filtered = []
    for seq in ha_cds_seqs:
        if seq.name not in added_names:
            filtered.append(seq)
            added_names.add(seq.name)
            start = ''
            start_tokens = ['ISUVDL', 'wgs']
            for token in start_tokens:
                if seq.contains_token(token):
                    start = token + '|'
                    break
            seq.full_name = start + '%s|%s|USA|%s' % (seq.name, seq.subtype if seq.subtype != 'unknown' else '',
                                                      seq.date.strftime('%Y-%m-%d'))
    helpers.write_iav_sequences('../../genes/human-swine-HA-Aug2021-pdms-cds.aln', filtered)


def standardize_name(seq: IAVSequence):
    if seq.host == 'human' and seq.subtype == 'unknown':
        seq.subtype = 'H1N1'
    start = ''
    start_tokens = ['ISUVDL', 'DARPA', 'wgs', 'VDL']
    for token in start_tokens:
        if seq.contains_token(token):
            if token == 'VDL' or token == 'DARPA':
                token = 'ISUVDL'
            start = token + '|'
            break
    seq.full_name = start + '%s|%s|USA|%s' % (seq.name, seq.subtype if seq.subtype != 'unknown' else '',
                                              seq.date.strftime('%Y-%m-%d'))


def standardize_swine_HA_fasta():
    swine_has = helpers.read_iav_sequences('../../pdm09-data/swine/Swine2010-2021_pdm09-HA_US_June2021-merged.fasta')
    added_names = set()
    filtered = []
    for seq in swine_has:
        if str(seq.seq).count('n') > 0:
            print('Unknown characters', seq.full_name, seq.seq.count('n'))
        if seq.name not in added_names:
            filtered.append(seq)
            added_names.add(seq.name)
            standardize_name(seq)
        else:
            print('Duplicate!', seq.full_name)
    # helpers.write_iav_sequences('../../pdm09-data/swine/Swine2010-2021_pdm09-HA_US_June2021-merged-std.fasta',
    #                             filtered)


def merge_gisaid_datasets():
    original_gisaid_filtered = helpers.read_iav_sequences('../../genes/Human2010-2021_pdm09-HA_US_gisaid-filtered.fasta')
    new_gisaid_filtred = helpers.read_iav_sequences('../../genes/Human2019-2021-Aug2021-US-gisaid-filtered.fasta')
    original_names = {iav.name for iav in original_gisaid_filtered}
    for new_iav in new_gisaid_filtred:
        if new_iav.name not in original_names and new_iav.date.month == 1 and new_iav.date.day == 1:
            print(new_iav.full_name)


def get_final_VDL_list(fasta_path='../pdm09-data/swine/Swine2009_2021-10_pdm_merged.fasta'):
    iavs = helpers.read_iav_sequences(fasta_path)
    vdl_iavs = [iav for iav in iavs if iav.full_name.count('VDL') > 0 or iav.full_name.count('DARPA') > 0]
    with open('../pdm09-data/swine/VDL_list.csv', 'w') as vdl_list:
        vdl_list.write('barcode, strain_name, note\n')
        vdl_seqs = []
        for iav in vdl_iavs:
            is_darpa = iav.full_name.count('DARPA') > 0
            if not is_darpa:
                barcode = f'{str(iav.date.year)}{int(iav.get_swine_name_id()):06}'
            else:
                barcode = f'{iav.date.year}{int(extract_active_barcode(iav.get_swine_name_id())):06}'
            print(iav.name, barcode, is_darpa)
            vdl_list.write(f'{barcode}, {iav.name}, {"DARPA" if is_darpa else ""}\n')
            record = iav.to_bio()
            record.id = barcode
            record.description = ''
            vdl_seqs.append(record)
        SeqIO.write(vdl_seqs, '../pdm09-data/swine/VDL_list.fasta', format='fasta')


def color_tree_w_hosts(tree_path='../trees/pdm09_US_all.timetree.hosts.tre',
                       timetree_path='../trees/pdm09_US_all_timetree/timetree.nexus',
                       #tree_path='../data/human-swine-pdms-cds-Sept2021.iqtree.pruned.hosts.tre',
                       clock_outliers_path='../trees/clock-outliers_Feb22_iqtree.txt', root_date=2008.72,
                       ha1_path='../trees/pdm09_US_all.ha1.aln',
                       ancestral_ha1_path='../trees/pdm09_US_all_aasub_ancestral/ancestral_sequences.fasta',
                       spillover_stats='../trees/spillover_stats.csv',
                       variants='../trees/variants.csv',
                       prune_outliers=True):
    timetree = Tree.get(path=timetree_path, schema='nexus', preserve_underscores=True)
    root_date = float(timetree.seed_node.annotations.get_value('date'))
    print('Inferred root date:', root_date)

    # to_prune = ["Utah/46106", "USA/67620", "Ohio/21141"]
    to_prune = []
    if prune_outliers:
        # clock_outliers = get_strain_name_list('../data/full_treetime/clock-outlier-names.txt')
        clock_outliers = get_strain_name_list(clock_outliers_path)
        to_prune = to_prune + clock_outliers
    tree = Tree.get(path=tree_path, schema='nexus', preserve_underscores=True)
    taxa_to_prune = [taxon.label for taxon in tree.taxon_namespace if
                     any(taxon.label.count(prune_token) > 0 for prune_token in to_prune)]
    ha1s = list(SeqIO.parse(ha1_path, format='fasta'))
    ancestral_ha1s = list(SeqIO.parse(ancestral_ha1_path, format='fasta'))
    print(f'Pruning {len(taxa_to_prune)} outliers')
    if taxa_to_prune:
        tree.prune_taxa_with_labels(taxa_to_prune)

    hu_to_sw, sw_to_hu, avg_dist, total_variants = 0, 0, 0, 0
    with open(variants, 'w') as variants_log:
        variants_log.write('hu-strain,spillover_id,spillover_season,time_drifted,ha1_drift\n')
        for node in tree.nodes():
            assert isinstance(node, Node)
            if node.annotations.get_value('host') == 'swine':
                node.annotations.add_new('!color', '#FFA500')

            # 2007.85
            if node.parent_node and node.parent_node.distance_from_root() + root_date >= 2009:
            # if node.parent_node:
                if node.annotations.get_value('host') == 'swine' and node.parent_node.annotations.get_value('host') == 'human':
                    hu_to_sw += 1
                if node.annotations.get_value('host') == 'human' and node.parent_node.annotations.get_value('host') == 'swine':
                    hu_ancestor = node.parent_node
                    time_lapsed = node.edge_length
                    while hu_ancestor.annotations.get_value('host') == 'swine':
                        hu_ancestor = hu_ancestor.parent_node
                        time_lapsed += hu_ancestor.edge_length
                    spillover_date = node.parent_node.distance_from_root() + root_date
                    spillover_year = int(spillover_date)
                    spillover_month = int((spillover_date - spillover_year) * 365 / 31) + 1  # Approximate month estimation.
                    season = pdm09.date_to_season(datetime.datetime(year=spillover_year, month=spillover_month, day=1),
                                                  swine=True)
                    print('Sw-to-hu spillover season:', season, spillover_year, spillover_month)
                    if len(node.leaf_nodes()) <= 20:
                        print(hu_ancestor.label)
                        ancestral_ha1 = [ha1 for ha1 in ancestral_ha1s if ha1.id == hu_ancestor.label][0]
                        for leaf in node.leaf_nodes():
                            leaf_ha1 = [ha1 for ha1 in ha1s if ha1.id == leaf.taxon.label][0]
                            dist = helpers.aligned_dist(ancestral_ha1.seq, leaf_ha1.seq, normalized=False, ignore_tails=False)
                            print(leaf.taxon.label, dist)
                            if leaf.annotations.get_value('host') == 'human' and season != '9-10':
                                avg_dist += dist
                                total_variants += 1
                            if leaf.annotations.get_value('host') == 'human':
                                variants_log.write(f'{leaf.taxon.label},{sw_to_hu + 1},{season},{time_lapsed},{dist}\n')
                    else:
                        print(len(node.leaf_nodes()))
                    print('-----------')
                    sw_to_hu += 1

    print('Human to swine:', hu_to_sw)
    print('Swine to human:', sw_to_hu)
    print('Avg HA1 dist: ', avg_dist / total_variants)
    with open(spillover_stats, 'w') as spillovers:
        spillovers.write('root-date,hu-to-sw,sw-to-hu\n')
        spillovers.write(f'{root_date},{hu_to_sw},{sw_to_hu}\n')

    tree.write_to_path(helpers.append_token_to_file_name(tree_path, 'colored'), schema='nexus')


def sample_zoonotic_subtrees(tree_path='../data/human-swine-pdms-cds-Sept2021.iqtree.pruned.hosts.tre'):
    tree = Tree.get(path=tree_path, schema='nexus', preserve_underscores=True)
    taxa_tokens = (['A/Iowa/01/2021', 'A02245873'],
                   ['A/Iowa/22/2020', 'A/swine/Iowa/A02479344/2020'])

    # Standardize names.
    for taxon in tree.taxon_namespace:
        parsed_seq = IAVSequence(taxon.label, 'atg', host='swine' if taxon.label.count('swine') > 0 else 'human')
        standardize_name(parsed_seq)
        taxon.label = parsed_seq.full_name

    # Extract subtrees.
    for i, lca_tokens in enumerate(taxa_tokens):
        lca_taxa = [taxon for taxon in tree.taxon_namespace if
                    any(taxon.label.count(lca_token) > 0 for lca_token in lca_tokens)]
        assert isinstance(tree, Tree)
        subtree_root = tree.mrca(taxa=lca_taxa)
        assert isinstance(subtree_root, Node)
        subtree = tree.extract_tree_with_taxa([leaf.taxon for leaf in subtree_root.leaf_nodes()])
        subtree.write_to_path('../data/zoonotic_subtree%d.tre' % (i + 1), schema='newick')


if __name__ == '__main__':
    # merge_vdl_with_existing()
    # merge_NA_pdms()
    # finalize_NAs()
    # finalize_HAs()
    # merge_gisaid_datasets()
    # standardize_swine_HA_fasta()
    # color_tree_w_hosts()
    # color_tree_w_hosts(tree_path='../data/treetime_hosts/treetime.hosts.tre')
    sample_zoonotic_subtrees()
    # pass
