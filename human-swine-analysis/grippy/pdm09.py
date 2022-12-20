# -*- coding: utf-8 -*-
from datetime import datetime
import matplotlib
from Bio import SeqIO
from dendropy import Node, Tree
from matplotlib import pyplot as plt
from sklearn import svm
import numpy as np
from typing import List
from random import shuffle

import grippy.helpers
from .similarity_matrix import HuSwSimilarityMatrix
from .sequence import IAVSequence
from . import helpers, pdm_data
from .helpers import read_iav_sequences, write_iav_sequences


def season_to_date_range(season='19-20', swine=False, midyear=False):
    year1, year2 = [int('20' + year_str) for year_str in season.split('-')]
    if midyear:
        return datetime(year=year1, month=7, day=1), datetime(year=year2, month=6, day=30)
    else:
        if swine:
            return datetime(year=year1, month=11, day=1), datetime(year=year2, month=10, day=31)
        else:
            return datetime(year=year1, month=9, day=1), datetime(year=year2, month=4, day=30)


def date_to_season(date: datetime, swine=False, midyear=False):
    y = date.year % 100
    test_seasons = ['%d-%d' % (y, y + 1), '%d-%d' % (y - 1, y)]
    for season in test_seasons:
        range = season_to_date_range(season, swine=swine, midyear=midyear)
        if range[0] <= date <= range[1]:
            return season
    return None


def next_season(season='18-19', step=1):
    year1, year2 = [int(year_str) for year_str in season.split('-')]
    return '%d-%d' % (year1 + step, year2 + step)


def sim_average(sim_list):
    return sum([x[0] for x in sim_list]) / len(sim_list)


def trim_sequences_by_reference(similarity_matrix: HuSwSimilarityMatrix,
                                sequences: List[IAVSequence], reference_token='A/Texas/11913/2019'):
    ref_seq = [hu_seq for hu_seq in similarity_matrix.human_sequences if hu_seq.contains_token(reference_token)][0]
    start_trim, end_trim = ref_seq.seq_start(), ref_seq.seq_end()
    trimmed_seqs = helpers.trim_alignment(sequences,
                                          start_pos=start_trim, end_pos=end_trim)
    return trimmed_seqs


# Returns a tuple: (all_sims, top p percent similarities, top n sim)
def hu_season_stats(hu_season: str, sw_seq_ind: int,
                    similarity_matrix: HuSwSimilarityMatrix, p=10, n=10):
    hu_season_range = season_to_date_range(hu_season, swine=False)
    sorted_sims = similarity_matrix.find_n_closest_cols(
        sw_seq_ind, n=100000,
        col_filter=(lambda i, seq: hu_season_range[0] <= seq.date <= hu_season_range[1])
    )
    top_p = sorted_sims[:int(round(len(sorted_sims) * p / 100.0))]
    top_n = sorted_sims[:n]
    return (sorted_sims, top_p, top_n)


def hu_strains_per_season(similarity_matrix: HuSwSimilarityMatrix):
    for hu_season in ['14-15', '15-16', '16-17', '17-18', '18-19', '19-20', '20-21']:
        all, _, _ = hu_season_stats(hu_season, 0, similarity_matrix)
        print('%s: %d' % (hu_season, len(all)))


def get_seqs_for_tree(swine_seq, hu_seq, start_date=datetime(year=2018, month=1, day=1),
                      end_date=datetime(year=2021, month=6, day=30), trim=True,
                      host_csv_path=None):
    filtered_swine = [seq for seq in swine_seq if start_date <= seq.date <= end_date]
    filtered_hu = [seq for seq in hu_seq if start_date <= seq.date <= end_date]
    merged_sequences = None
    if trim:
        ref_seq = [hu_seq for hu_seq in hu_seq if hu_seq.contains_token('A/Texas/11913/2019')][0]
        start_trim, end_trim = ref_seq.seq_start(), ref_seq.seq_end()
        trimmed_seqs = helpers.trim_alignment(filtered_swine + filtered_hu,
                                              start_pos=start_trim, end_pos=end_trim)
        merged_sequences = trimmed_seqs
        # write_iav_sequences('hu-sw-pdm.fasta', trimmed_seqs)
    else:
        merged_sequences = filtered_swine + filtered_hu
    if host_csv_path:
        csv_output = open(host_csv_path, 'w+')
        csv_output.write('name, host\n')
        for iav in merged_sequences:
            csv_output.write('%s, %s\n' % (iav.full_name, iav.host))
        csv_output.close()
    return merged_sequences


main_tree = None
saved_season_spillovers = {}


def find_spillovers_phylogenetic(season, similarity_matrix: HuSwSimilarityMatrix,
                                 tree_path='../trees/pdm09_US_all.timetree.hosts.tre', swine_season=None):
    """
    returns a set of swine sequences' names that are spillovers in 'season'.
    :param swine_season: The season for swine sequences to be considered.
                         if None, the swine season is assumed to be the same as 'season'.
    """

    def hu_season_mrcas(tree: Tree, node: Node):
        """
        Recursively finds all nodes that are maximal common ancestors to only those human sequences
        that are from current season or later.
        :return: A pair (bool, int). 1: True if all human sequences below 'node' are from 'season' or later.
                 2: Number of human leaves below the node
        """
        if node.is_leaf():
            iav_seq = node.iav_seq
            hu_season_range = season_to_date_range(season, swine=False)
            if iav_seq.host == 'human':
                return (iav_seq.date >= hu_season_range[0], 1)
            else:
                return (True, 0)
        else:
            children_status = []
            mixed = False
            humans_below = 0
            for child in node.child_nodes():
                current_only, ch_humans_below = hu_season_mrcas(tree, child)
                children_status.append((child, current_only, ch_humans_below))
                mixed = mixed or (not current_only)
                humans_below += ch_humans_below
            if mixed and humans_below > 0:
                for child, current_only, ch_humans_below in children_status:
                    if current_only and ch_humans_below > 0:
                        season_mrcas.append((child, humans_below))
            return (not mixed, humans_below)

    # global saved_season_spillovers
    # if season in saved_season_spillovers.keys() and not swine_season:
    #     return saved_season_spillovers[season]
    start_date = season_to_date_range(next_season(season, step=-1), swine=False)[0]
    end_date = season_to_date_range(season, swine=True)[1]
    # alignment_path = 'data/tree-%s.fasta' % season
    # if not os.path.exists(alignment_path):
    #     tree_sequences = get_seqs_for_tree(similarity_matrix.swine_sequences, similarity_matrix.human_sequences,
    #                                        start_date, end_date)
    #     write_iav_sequences(alignment_path, tree_sequences)
    # else:
    #     tree_sequences = get_seqs_for_tree(similarity_matrix.swine_sequences, similarity_matrix.human_sequences,
    #                                        start_date, end_date, trim=False)
    # global main_tree
    # if not main_tree:
    #     # tree_path = 'data/pdm09_all_June2021.rooted.hosts.nexus'
    #     # tree_path = 'data/human-swine-pdms-cds-Sept2021.iqtree.pruned.hosts.tre'
    #     # tree_path = 'data/treetime_hosts/treetime.hosts.tre'
    #     tree_path = '../trees/pdm09_US_all.timetree.hosts.tre'
    #     tree_sequences = similarity_matrix.swine_sequences + similarity_matrix.human_sequences
    #     tree = helpers.get_annotated_tree(tree_path, tree_sequences, schema='nexus')
    #     main_tree = tree
    # else:
    #     tree = main_tree
    tree_sequences = similarity_matrix.swine_sequences + similarity_matrix.human_sequences
    tree = helpers.get_annotated_tree(tree_path, tree_sequences, schema='nexus')

    # Identify most recent common ancestors for this season's human strains.
    season_mrcas = []
    root_current_only, _ = hu_season_mrcas(tree, tree.seed_node)
    if root_current_only:
        season_mrcas.append((tree.seed_node, len(similarity_matrix.human_sequences)))

    # Label all swine sequences below nodes in 'season_mrcas' as same-season spillovers.
    same_season_spills = set()
    if swine_season:
        sw_season_range = season_to_date_range(swine_season, swine=True)
    else:
        sw_season_range = season_to_date_range(season, swine=True)
    # print('Overall leaves: %d' % len(tree.leaf_nodes()))
    leaves_traversed = 0
    for hu_season_mrca, humans_below in season_mrcas:
        percent_human = humans_below / len(hu_season_mrca.leaf_nodes())
        # Skip a small clade with few human seqs.
        # if humans_below <= 10 and percent_human < 0.7:
        #     continue
        # if len(hu_season_mrca.leaf_nodes()) > 0:
        #     print('MRCA leaves: %d' % len(hu_season_mrca.leaf_nodes()))
        if hu_season_mrca.annotations['host'].value == 'swine' or\
           humans_below < 0.7 * len(hu_season_mrca.leaf_nodes()):
            continue
        if hu_season_mrca.annotations['host'].value == 'swine' and\
           humans_below >= 0.7 * len(hu_season_mrca.leaf_nodes()):
            print('>70% human, but annotated as swine.')
        earliest_iav = None
        for leaf in hu_season_mrca.leaf_iter():
            if not earliest_iav or leaf.iav_seq.date < earliest_iav.date:
                earliest_iav = leaf.iav_seq
        # if earliest_iav.host is not 'human':
        #     continue
        # if humans_below < 0.7 * len(hu_season_mrca.leaf_nodes()):
        #     # We only consider those clades, where # of humans is at least 70%,
        #     # or the earliest sampled leaf is a human IAV.
        #     earliest_iav = None
        #     for leaf in hu_season_mrca.leaf_iter():
        #         if not earliest_iav or leaf.iav_seq.date < earliest_iav.date:
        #             earliest_iav = leaf.iav_seq
        #     if earliest_iav.host is not 'human':
        #         continue

        for leaf in hu_season_mrca.leaf_iter():
            leaves_traversed += 1
            if leaf.iav_seq.host == 'swine' and sw_season_range[0] <= leaf.iav_seq.date <= sw_season_range[1]:
                same_season_spills.add(leaf.iav_seq.full_name)
                if season == '20-21':
                    print(leaf.iav_seq.full_name, humans_below, percent_human)
    # print('Leaves Traversed: %d' % leaves_traversed)
    labeled_seqs = [(seq, 1 if seq.full_name in same_season_spills else 0) for seq in similarity_matrix.swine_sequences
                    if sw_season_range[0] <= seq.date <= sw_season_range[1]]
    zeros = len([0 for seq, label in labeled_seqs if label == 0])
    print('Season %s. Zeroes: %d, Ones: %d' % (season, zeros, len(labeled_seqs) - zeros))
    if not swine_season:
        saved_season_spillovers[season] = same_season_spills
    return same_season_spills


def identify_closest(similarity_matrix: HuSwSimilarityMatrix):
    # seasons = ['16-17', '17-18', '18-19']
    season = '11-19'
    sw_date_range = season_to_date_range(season, swine=True)
    hu_date_range = season_to_date_range(season, swine=False)
    above_nn, total, in_season, prior = 0, 0, 0, 0
    for ind, sw_seq in enumerate(similarity_matrix.swine_sequences):
        if not sw_date_range[0] <= sw_seq.date <= sw_date_range[1]:
            continue
        sim, human_seq = similarity_matrix.find_n_closest_cols(ind, n=1)[0]
        total += 1
        if sim >= 0.99:
            above_nn += 1
            # print(human_seq.date.year)
            # print(sw_seq.date, sim, human_seq.date)
        if hu_date_range[0] <= human_seq.date <= hu_date_range[1]:
            in_season += 1
        if human_seq.date <= sw_seq.date:
            prior += 1
    print("Above 99%%: %d/%d; Human before: %d; In-season: %d" % (above_nn, total, prior, in_season))


def get_training_sequences(similarity_matrix: HuSwSimilarityMatrix, window=100, all=False):
    if all:
        seasons = ['13-14', '14-15', '15-16', '16-17', '17-18', '18-19', '19-20']
        prev_spills = []
        same_season_spills = []
        for season in seasons:
            same_season_spill_names = find_spillovers_phylogenetic(season, similarity_matrix)
            sw_season_range = season_to_date_range(season, swine=True)
            for swine_seq in similarity_matrix.swine_sequences:
                if sw_season_range[0] <= swine_seq.date <= sw_season_range[1]:
                    if swine_seq.full_name in same_season_spill_names:
                        same_season_spills.append(swine_seq)
                    else:
                        prev_spills.append(swine_seq)
        helpers.write_iav_sequences('data/spillover-training-all-recent.fasta', same_season_spills)
        helpers.write_iav_sequences('data/spillover-training-all-prev.fasta', prev_spills)
    else:
        seasons = ['18-19']
        season_spills = []
        continued_spills = []

        for season in seasons:
            hu_season_range = season_to_date_range(season, swine=False)
            season_spill_names = find_spillovers_phylogenetic(season, similarity_matrix)
            for spill_name in season_spill_names:
                # We find a closest human sequence and make sure its within the last
                # 'window' days of 'spill_name'.
                season_spill = similarity_matrix.swine_sequences[
                    similarity_matrix.get_row_id_by_name(spill_name)]
                closest_human = similarity_matrix.find_n_closest_cols(
                    similarity_matrix.get_row_id_by_name(season_spill), n=1,
                    col_filter=(lambda i, seq: hu_season_range[0] <= seq.date <= season_spill.date),
                    latest_first=True
                )[0][1]
                # if (season_spill.date - closest_human.date).days <= window:
                if season_spill.date <= datetime(year=2019, month=4, day=30)\
                        and (season_spill.date - closest_human.date).days <= window:
                    season_spills.append(season_spill)
            subsequent_spill_names = find_spillovers_phylogenetic(season, similarity_matrix,
                                                                  swine_season=next_season(season, step=1))
            next_spill_names = find_spillovers_phylogenetic(next_season(season, step=1), similarity_matrix)
            continued_circ_names = list(set(subsequent_spill_names) - set(next_spill_names))
            for spill_name in continued_circ_names:
                spill = similarity_matrix.swine_sequences[
                    similarity_matrix.get_row_id_by_name(spill_name)]
                continued_spills.append(spill)
        print('Recent spills: %d; Continued circulation: %d' % (len(season_spills), len(continued_spills)))
        trimmed_sequences = trim_sequences_by_reference(similarity_matrix, season_spills + continued_spills)
        trimmed_same_season = trimmed_sequences[:len(season_spills)]
        trimmed_continued = trimmed_sequences[len(season_spills):]
        # Randomly subsample continued circulation spills to match the number of recent spills.
        shuffle(trimmed_continued)
        trimmed_continued = trimmed_continued[:len(trimmed_same_season)]

        helpers.write_iav_sequences('data/spillover-training-recent.fasta', trimmed_same_season)
        helpers.write_iav_sequences('data/spillover-training-continued.fasta', trimmed_continued)
        training_output = open('data/spillover-training-2.csv', 'w+')
        for season_spill in season_spills:
            training_output.write('%s, %d\n' % (season_spill.full_name, 1))
        for prev_season_spill in continued_spills:
            training_output.write('%s, %d\n' % (prev_season_spill, 0))
        training_output.close()


def process_gisaid(path: str, output_path: str):
    # 1. Read in all strains.
    iav_sequences = read_iav_sequences(path)
    # 2. Filter swine-passaged ones out.
    filtered = [seq for seq in iav_sequences if not seq.in_header('_swl')]
    final = []
    unique_names = {seq.name for seq in filtered}
    print(len(unique_names))
    # 3. Keep at most one sequence per unique name -- ideally original (non-passaged) and non 20xx-01-01.
    for name in unique_names:
        seqs = [seq for seq in iav_sequences if seq.name == name]
        if len(seqs) > 1:
            original = {seq for seq in seqs if seq.contains_token('Original')}
            non_jan_first = {seq for seq in seqs if not (seq.date.month == 1 and seq.date.day == 1)}
            has_original = len(original) > 0
            if len(original.intersection(non_jan_first)) > 0:
                final.append(original.intersection(non_jan_first).pop())
            elif len(original) > 0:
                final.append(original.pop())
            elif len(non_jan_first) > 0:
                final.append(non_jan_first.pop())
            else:
                final.append(seqs[0])
        else:
            final.append(seqs[0])
    final_bio = [seq.to_bio() for seq in final]
    SeqIO.write(final_bio, output_path, 'fasta')


def save_seasonal_sequences(similarity_matrix: HuSwSimilarityMatrix, seasons=('17-18', '18-19', '19-20')):
    for season in seasons:
        hu_season_range = season_to_date_range(season, swine=False)
        sw_season_range = season_to_date_range(season, swine=True)
        hu_iavs = [hu_iav for hu_iav in similarity_matrix.human_sequences
                   if hu_season_range[0] <= hu_iav.date <= hu_season_range[1]]
        sw_iavs = [sw_iav for sw_iav in similarity_matrix.swine_sequences
                   if sw_season_range[0] <= sw_iav.date <= sw_season_range[1]]
        helpers.write_iav_sequences('data/hu_%s.fasta' % season, hu_iavs)
        helpers.write_iav_sequences('data/sw_%s.fasta' % season, sw_iavs)


def save_all_20_21_sequences_trimmed(similarity_matrix: HuSwSimilarityMatrix):
    sw_season_range = season_to_date_range('20-21', swine=True)
    sw_seqs = [sw_seq for sw_seq in similarity_matrix.swine_sequences if
               sw_season_range[0] <= sw_seq.date <= sw_season_range[1]]
    ref_seq = [hu_seq for hu_seq in similarity_matrix.human_sequences if hu_seq.contains_token('A/Texas/11913/2019')][0]
    start_trim, end_trim = ref_seq.seq_start(), ref_seq.seq_end()
    trimmed_sw_seqs = helpers.trim_alignment(sw_seqs, start_pos=start_trim, end_pos=end_trim)
    helpers.write_iav_sequences('data/20-21-all.fasta', trimmed_sw_seqs)


def predict_20_21(similarity_matrix: HuSwSimilarityMatrix):
    human_season = '19-20'
    hu_season_range = season_to_date_range(human_season, swine=False)
    hu_seqs = [hu_seq for hu_seq in similarity_matrix.human_sequences if
               hu_season_range[0] <= hu_seq.date <= hu_season_range[1]]
    ref_seq = [hu_seq for hu_seq in hu_seqs if hu_seq.contains_token('A/Texas/11913/2019')][0]
    start_trim, end_trim = ref_seq.seq_start(), ref_seq.seq_end()
    trimmed_hu_seqs = helpers.trim_alignment(hu_seqs, start_pos=start_trim, end_pos=end_trim)
    # write_iav_sequences('19-20-human.fasta', trimmed_hu_seqs)

    svm_reg = 75

    tt_pred_labels, tt_prev_seasons, tt_seqs, _ = \
        infer_spillover_seasons(similarity_matrix, ['13-14', '14-15', '15-16', '16-17', '17-18', '18-19', '19-20'], predict_season='20-21',
                                svm_reg=svm_reg, plot=True)
    for i, pred_label in enumerate(tt_pred_labels):
        if pred_label:
            closest_human = similarity_matrix.find_n_closest_cols(
                similarity_matrix.get_row_id(tt_seqs[i]), n=1)[0][1]
            print(tt_seqs[i].full_name, tt_prev_seasons[i, 0], tt_prev_seasons[i, 1], closest_human.full_name)
    predicted_swine = [seq for i, seq in enumerate(tt_seqs) if tt_pred_labels[i]]
    trimmed_predicted = helpers.trim_alignment(predicted_swine, start_pos=start_trim, end_pos=end_trim)
    # write_iav_sequences('20-21-predicted.fasta', helpers.sort_by_date(trimmed_predicted))


def is_same_season_spillover(sw_seq: IAVSequence, seq_ind: int, season: str,
                             similarity_matrix: HuSwSimilarityMatrix, phylogenetic=True) -> bool:
    if phylogenetic:
        return sw_seq.full_name in find_spillovers_phylogenetic(season, similarity_matrix)
    else:
        cur_top_all, cur_top_p, cur_top_n = hu_season_stats(season, seq_ind, similarity_matrix, p=5, n=10)
        prev_top_all, prev_top_p, prev_top_n = hu_season_stats(next_season(season, step=-1),
                                                               seq_ind, similarity_matrix, p=5, n=10)
        prev2_top_all, prev2_top_p, prev2_top_n = hu_season_stats(next_season(season, step=-2),
                                                                  seq_ind, similarity_matrix, p=5, n=10)
        if cur_top_all[0][0] >= 0.985 and cur_top_p[-1][0] > prev_top_p[-1][0] and \
                cur_top_p[-1][0] > prev2_top_p[-1][0]:
            return True
        return False


def find_same_season_spillovers(similarity_matrix: HuSwSimilarityMatrix, phylogenetic=False):
    seasons = ['11-12', '12-13', '13-14', '14-15', '15-16', '16-17', '17-18', '18-19', '19-20']
    total_spill, percent_spill = [], []
    for season in seasons:
        sw_season_range = season_to_date_range(season, swine=True)
        spills, total_seq = 0, 0
        for ind, swine_seq in enumerate(similarity_matrix.swine_sequences):
            if sw_season_range[0] <= swine_seq.date <= sw_season_range[1]:
                total_seq += 1
                if is_same_season_spillover(swine_seq, ind, season, similarity_matrix, phylogenetic=phylogenetic):
                    spills += 1
        print('Season %s: spills %d (%.3f)' % (season, spills, spills / total_seq))
        total_spill.append(spills)
        percent_spill.append(spills / total_seq)
    print(', '.join(['%.3f' % p for p in percent_spill]))


def print_all_ggt(similarity_matrix: HuSwSimilarityMatrix):
    signature = {969: 'g', 981: 'g', 1405: 'c'}
    # signature = {981: 'a'}
    hu_seasons = {'14-15': 0, '15-16': 0, '16-17': 0, '17-18': 0, '18-19': 0, '19-20': 0, '20-21': 0}
    for iav in similarity_matrix.human_sequences + similarity_matrix.swine_sequences:
        has_signature = True
        for pos in signature.keys():
            if iav.seq[pos - 1].lower() != signature[pos]:
                has_signature = False
                break
        if has_signature:
            season = date_to_season(iav.date)
            if iav.host == 'human' and season and season in hu_seasons.keys():
                hu_seasons[season] += 1
            # print(iav.host, date_to_season(iav.date), iav.full_name)
    print(hu_seasons)


def validate_predictions(similarity_matrix: HuSwSimilarityMatrix, phylogenetic=True):
    # validation_seasons = {'13-14', '17-18', '15-16', '18-19', '19-20'}
    validation_seasons = {'13-14', '14-15', '15-16', '16-17', '17-18', '18-19', '19-20'}
    svm_reg = 75

    # Do leave one season out analysis.
    avg_error, avg_fp = 0, 0
    avg_precision, avg_recall, avg_f1 = 0, 0, 0
    for test_season in validation_seasons:
        train_seasons = validation_seasons - {test_season}
        predicted, _, _, true_labels =\
            infer_spillover_seasons(similarity_matrix, train_seasons, predict_season=test_season,
                                    svm_reg=svm_reg, plot=False, phylogenetic=phylogenetic,
                                    positive_weight=2)
        fn, fp, tn, tp = 0, 0, 0, 0
        for i, label in enumerate(true_labels):
            p_label = predicted[i]
            if p_label == label == 1:
                tp += 1
            if p_label == label == 0:
                tn += 1
            if p_label == 1 and label == 0:
                fp += 1
            if p_label == 0 and label == 1:
                fn += 1
        print('Season %s. FP: %d, TP: %d, FN: %d, TN: %d' % (test_season, fp, tp, fn, tn))
        if len(true_labels) == 0 or tp == 0:
            continue
        precision = tp / (fp + tp)
        recall = tp / (tp + fn) if (tp + fn) > 0 else -1
        f1 = tp / (tp + 0.5 * (fp + fn))
        error_rate = (fn + fp) / len(true_labels)
        fp_rate = fp / len(true_labels)
        print('Error: %.2f; FP rate: %.2f; Precision: %.2f; Recall: %.2f; F1 : %.3f; Ratio: %.3f' %
              (error_rate, fp_rate, precision, recall, f1, (tp + fn) / len(true_labels)))
        avg_error += error_rate
        avg_fp += fp_rate
        avg_precision += precision
        avg_recall += recall
        avg_f1 += f1
    avg_error /= len(validation_seasons)
    avg_fp /= len(validation_seasons)
    avg_precision /= len(validation_seasons)
    avg_recall /= len(validation_seasons)
    avg_f1 /= len(validation_seasons)
    print('Avg Error: %.2f; Avg FP: %.2f; Avg Precision: %.3f; Avg Recall: %.3f; Avg F1: %.3f' %
          (avg_error, avg_fp, avg_precision, avg_recall, avg_f1))

    tt_pred_labels, tt_prev_seasons, tt_seqs, _ =\
        infer_spillover_seasons(similarity_matrix, validation_seasons, predict_season='20-21',
                                svm_reg=svm_reg, positive_weight=3, plot=True, phylogenetic=phylogenetic)
    for i, pred_label in enumerate(tt_pred_labels):
        if pred_label:
            closest_human = similarity_matrix.find_n_closest_cols(
                similarity_matrix.get_row_id(tt_seqs[i]), n=1)[0][1]
            print(tt_seqs[i].full_name, tt_prev_seasons[i, 0], tt_prev_seasons[i, 1], closest_human.full_name)
    # write_iav_sequences('20-21-predicted.fasta', [seq for i, seq in enumerate(tt_seqs) if tt_pred_labels[i]])


saved_prev_swine_similarities = {}


def infer_spillover_seasons(similarity_matrix: HuSwSimilarityMatrix, tree_path='../trees/pdm09_US_all.timetree.hosts.tre',
                            training_seasons=('17-18', '18-19', '19-20'), predict_season='20-21',
                            svm_reg=50, positive_weight=2, phylogenetic=False, plot_predictors=True,
                            plot=True, last_spill_only=False, print_predicted=False, correct_p=None):
    # seasons = ['17-18', '18-19', '19-20']
    # trainin_seasons = ['17-18']
    # sw_date_range = (datetime(year=2011, month=9, day=1), datetime(year=2020, month=3, day=30))
    prev_season_sims = []
    prev2_sims = []
    prev_sims = []
    max_prev_sims = []
    prev_sw_sims = []
    prev_normalized_sims = []
    cur_spill_labels = []
    in_season_best = 0
    train_swine = 0
    predict_prev2_sims = []
    predict_prev_sims = []
    predict_max_prev_sims = []
    predict_prev_sw_sims = []
    predict_prev_normalized = []
    predict_seqs = []
    train_seqs = []
    true_labels = []
    all_seasons = (list(training_seasons) + [predict_season]) if predict_season else list(training_seasons)
    for i, season in enumerate(all_seasons):
        if phylogenetic:
            same_season_spills = find_spillovers_phylogenetic(season, similarity_matrix, tree_path=tree_path)
        sw_season_range = season_to_date_range(season, swine=True)
        prev_sw_season_range = season_to_date_range(next_season(season, step=-1), swine=True)
        hu_season_range = season_to_date_range(season, swine=False)
        print('Human sequences in season %s:' % season, len([seq for seq in similarity_matrix.human_sequences
                                                             if hu_season_range[0] <= seq.date <= hu_season_range[1]]))
        if season == '20-21':
            for seq in [seq for seq in similarity_matrix.human_sequences if hu_season_range[0] <= seq.date <= hu_season_range[1]]:
                print(f'\t{seq.full_name}')
        for ind, swine_seq in enumerate(similarity_matrix.swine_sequences):
            if sw_season_range[0] <= swine_seq.date <= sw_season_range[1]:
                p = correct_p if correct_p else 5
                cur_top_all, cur_top_p, cur_top_n = hu_season_stats(season, ind, similarity_matrix, p=p, n=10)
                prev_top_all, prev_top_p, prev_top_n = hu_season_stats(next_season(season, step=-1),
                                                                       ind, similarity_matrix, p=p, n=10)
                prev2_top_all, prev2_top_p, prev2_top_n = hu_season_stats(next_season(season, step=-2),
                                                                          ind, similarity_matrix, p=p, n=10)
                prev_season_swine = [sw_seq for sw_seq in similarity_matrix.swine_sequences
                                     if prev_sw_season_range[0] <= sw_seq.date <= prev_sw_season_range[1]]
                prev_normalized = (1-prev_top_all[0][0]) / (swine_seq.date - prev_top_all[0][1].date).days
                global saved_prev_swine_similarities
                if saved_prev_swine_similarities.get(swine_seq.full_name):
                    prev_swine_sim = saved_prev_swine_similarities[swine_seq.full_name]
                else:
                    prev_swine_sim, prev_swine_closest = helpers.find_closest_sequence(swine_seq, prev_season_swine,
                                                                                       ignore_tails=False)
                    prev_swine_sim = float('%.6f' % prev_swine_sim)
                    saved_prev_swine_similarities[swine_seq.full_name] = prev_swine_sim
                # if season_of_interest == season:
                #     print(prev_top_all[0][1], prev_top_all[0][0], prev_top_all[1][1], prev_top_all[1][0])
                if predict_season and i == len(all_seasons) - 1:
                    predict_seqs.append(swine_seq)
                    predict_prev2_sims.append(sim_average(prev2_top_all[:1]))
                    if correct_p:
                        predict_prev_sims.append(sim_average(prev_top_p[-1:]))  # Use (100-correct_p)th percentile value.
                    else:
                        predict_prev_sims.append(sim_average(prev_top_all[:1]))  # Otherwise use max similarity
                    predict_prev_sw_sims.append(prev_swine_sim)
                    predict_prev_normalized.append(prev_normalized)
                    if prev_swine_sim <= prev_top_all[0][0]:
                        print(swine_seq.full_name)
                else:
                    if last_spill_only:
                        if not sim_average(prev_top_p[:1]) >= sim_average(prev2_top_p[:1]):
                            continue

                    train_swine += 1
                    # has_next_season = False
                    # if season != '19-20' and season != '20-21':
                    #     has_next_season = True
                    #     next_top_all, next_top_p, next_top_n = hu_season_stats(next_season(season, step=1),
                    #                                                            ind, similarity_matrix, p=10, n=10)
                    prev2_sims.append(sim_average(prev2_top_all[:1]))
                    if correct_p:
                        prev_sims.append(sim_average(prev_top_p[-1:]))  # Use (100-correct_p)th percentile value.
                    else:
                        prev_sims.append(sim_average(prev_top_all[:1]))  # Otherwise use max similarity
                    prev_sw_sims.append(prev_swine_sim)
                    prev_normalized_sims.append(prev_normalized)
                    train_seqs.append(swine_seq)

                if len(cur_top_all) >= 20:
                    cur_spillover = 0
                    if phylogenetic:
                        cur_spillover = swine_seq.full_name in same_season_spills
                        if cur_spillover:
                            if prev_top_all[0][0] <= prev_swine_sim:
                                if prev_top_all[0][0] == prev_swine_sim:
                                    print('Same as swine:', swine_seq.full_name, prev_top_all[0][0], prev_swine_sim)
                                else:
                                    print('Closer to swine:', swine_seq.full_name, prev_top_all[0][0], prev_swine_sim)
                            prev_season_sims.append(())
                    else:
                        # if avg_sim > avg_prev_sim and avg_sim > next_avg_sim:
                        if cur_top_all[0][0] >= 0.985 and sim_average(cur_top_p[-1:]) > sim_average(prev_top_p[-1:]) and\
                                sim_average(cur_top_p[-1:]) > sim_average(prev2_top_p[-1:]):
                            # We assume this means that the spillover happened this season
                            cur_spillover = 1
                            prev_season_sims.append(())
                    if predict_season and i == len(all_seasons) - 1:
                        true_labels.append(cur_spillover)
                    else:
                        cur_spill_labels.append(cur_spillover)

    print('%d/%d' % (len(prev_season_sims), train_swine))
    # plt.hist(x=[y[0] for y in prev_season_sims], bins='auto', color='#0504aa',
    #          alpha=0.5, rwidth=0.85)
    # plt.hist(x=[y[0] for y in prev_season_sims], bins='auto', color='green',
    #          alpha=0.5, rwidth=0.85)
    # plt.hist(x=[y[1] for y in prev_season_sims], bins='auto', color='red',
    #          alpha=0.5, rwidth=0.85)
    # train_set = np.array([prev2_sims, prev_sims])
    # train_set = np.array([prev2_sims, prev_normalized_sims])
    # train_set = np.array([prev2_sims, prev_sims, prev_sw_sims])
    train_set = np.array([prev_sims, prev_sw_sims])
    train_set = train_set.transpose()
    # predict_set = np.array([predict_prev2_sims, predict_prev_sims])
    # predict_set = np.array([predict_prev2_sims, predict_prev_normalized])
    # predict_set = np.array([predict_prev2_sims, predict_prev_sims, predict_prev_sw_sims])
    predict_set = np.array([predict_prev_sims, predict_prev_sw_sims])
    predict_set = predict_set.transpose()
    train_labels = np.array(cur_spill_labels)
    if plot_predictors:
        plot_train_points(train_set, train_labels, predictors=predict_set)
    predicted_labels = train_and_predict(train_set, predict_set, train_labels, train_seqs, svm_reg=svm_reg,
                                         positive_weight=positive_weight, plot=plot)
    if print_predicted:
        for i, label in enumerate(predicted_labels):
            if label == 1:
                print(predict_seqs[i], predict_set[i, 0], predict_set[i, 1])
    return predicted_labels, predict_set, predict_seqs, true_labels


def train_and_predict(train_set, predict_set, train_labels, train_seqs, svm_reg=50, positive_weight=2,
                      plot=True):
    colors = ['red', 'blue']
    # plt.scatter(x=X[:, 0], y=X[:, 1], c=y,
    #             cmap=matplotlib.colors.ListedColormap(colors))
    # plt.show()
    # return 1
    svc = svm.SVC(kernel='rbf', C=svm_reg, class_weight={0: positive_weight, 1: 1}, probability=False).fit(train_set, train_labels)
    print('SVM fitted. Support vectors: %d' % len(svc.support_vectors_))

    predicted = svc.predict(train_set)
    fn, fp, tn, tp = 0, 0, 0, 0
    for i, label in enumerate(train_labels):
        p_label = predicted[i]
        if p_label == label == 1:
            tp += 1
        if p_label == label == 0:
            tn += 1
        if p_label == 1 and label == 0:
            # print(train_seqs[i])
            fp += 1
        if p_label == 0 and label == 1:
            fn += 1
    print('FP: %d, TP: %d, FN: %d, TN: %d' % (fp, tp, fn, tn))
    if predict_set.all():
        predicted_labels = svc.predict(predict_set)
        print('predicted spillovers: %d/%d' % (np.count_nonzero(predicted_labels), predict_set.shape[0]))

    if plot:
        h=0.001
        print(max(train_set[:, 1]))
        xx, yy = np.meshgrid(np.arange(0.92, 1.01, h),
                             np.arange(0.92, 1.01, h))
        Z = svc.predict(np.c_[xx.ravel(), yy.ravel()])
        # # Put the result into a color plot
        Z = Z.reshape(xx.shape)
        plt.contourf(xx, yy, Z, cmap=plt.cm.coolwarm, alpha=0.6)
        # Plot also the training points
        plt.scatter(train_set[:, 0], train_set[:, 1], c=train_labels, cmap=plt.cm.coolwarm, s=0.5)
        plt.plot([0, 1], [0, 1], 'k-', color='black')
        # plt.ylim([0., max(train_set[:, 1])])
        plt.ylim([0.975, 1.0])
        plt.xlim([0.975, 1.0])
        # plt.show()
        plt.scatter(predict_set[:, 0], predict_set[:, 1], color='yellow', s=0.5)
        plt.show()

    return predicted_labels


def plot_train_points(train_set, train_labels, predictors=None):
    colors = ['blue', 'red']
    matplotlib.rcParams.update({'font.size': 12})
    fig1 = plt.figure(1, figsize=(6.5, 6))
    fig1.tight_layout()
    train_0 = train_set[train_labels == 0, :]
    train_1 = train_set[train_labels == 1, :]
    # train_0 = np.array([row for i, row in enumerate(train_set) if train_labels[i] == 0])
    # train_1 = np.array([row for i, row in enumerate(train_set) if train_labels[i] == 1])
    # plt.scatter(train_set[:, 0], train_set[:, 1], c=train_labels, cmap=matplotlib.colors.ListedColormap(colors), s=3)
    plt.scatter(train_0[:, 0], train_0[:, 1], c='blue', s=12, marker='o')
    plt.scatter(train_1[:, 0], train_1[:, 1], s=12, marker='o', facecolors='none', edgecolors='red')
    plt.plot([0, 1], [0, 1], 'k-', color='black', linewidth=0.5)
    # plt.ylim([0., max(train_set[:, 1])])
    plt.ylim([0.975, 1.0])
    plt.xlim([0.975, 1.0])
    plt.xlabel('Similarity to prior human season')
    plt.ylabel('Similarity to prior swine season')
    if predictors is not None:
        fig2 = plt.figure(2, figsize=(6.5, 6))
        fig2.tight_layout()
        plt.scatter(predictors[:, 0], predictors[:, 1], color='dimgray', s=12)
        plt.plot([0, 1], [0, 1], 'k-', color='black', linewidth=0.5)
        plt.ylim([0.975, 1.0])
        plt.xlim([0.975, 1.0])
        plt.xlabel('Similarity to prior human season')
        plt.ylabel('Similarity to prior swine season')

    plt.show()


def get_dates(dates_path='../trees/ancestral_dates.csv'):
    node_dates = {}
    with open(dates_path, 'r') as dates_file:
        header_skipped = False
        for line in dates_file:
            if not header_skipped:
                header_skipped = True
                continue
            node, date_str = [token.strip() for token in line.split(',')[:2]]
            if date_str.count('--') > 0:
                continue
            date = datetime.strptime(date_str, '%Y-%m-%d')
            node_dates[node] = date
    return node_dates


def assess_long_spillovers(tree_path='data/treetime_hosts/treetime.hosts.tre',
                           timetree_path='../trees/pdm09_US_all_timetree/timetree.nexus',
                           log_path='../trees/spillover_lengths.csv',
                           # dates_path='../trees/ancestral_dates.csv', root_age=2007.85,
                           threshold=1,
                           seasons=('9-10', '10-11', '11-12', '12-13', '13-14', '14-15', '15-16', '16-17', '17-18', '18-19', '19-20'),
                           outlier_path='../trees/clock-outliers_Feb22_iqtree.txt', prune_outliers=True):
    # Find root age
    timetree = Tree.get(path=timetree_path, schema='nexus', preserve_underscores=True)
    root_age = float(timetree.seed_node.annotations.get_value('date'))

    # Open log and write the header.
    spillover_log = open(log_path, 'w')
    for (i, season) in enumerate(seasons):
        spillover_log.write(f'{season}_num,{season}_plong,{season}_pslong,{season}_max,{season}_avg')
        if i < len(seasons) - 1:
            spillover_log.write(',')
    spillover_log.write('\n')

    # Read a time-scaled phylogeny and prune outliers if needed.
    tree = Tree.get(path=tree_path, schema='nexus', preserve_underscores=True)
    if prune_outliers:
        clock_outliers = grippy.helpers.get_strain_name_list(outlier_path)
        taxa_to_prune = [taxon.label for taxon in tree.taxon_namespace if
                         any(taxon.label.count(prune_token) > 0 for prune_token in clock_outliers)]
        print(f'Pruning {len(taxa_to_prune)} outliers')
        if taxa_to_prune:
            tree.prune_taxa_with_labels(taxa_to_prune)

    # node_dates = get_dates(dates_path)
    spillovers_by_season = {}
    zoonotic_spillovers = {season: [] for season in (list(seasons) + ['20-21'])}
    long_spillovers_by_season = {}
    long_sustained_spills_by_season = {}
    total_spill_length = {}
    longest_spillovers = {}
    for node in tree.nodes():
        parent = node.parent_node
        if parent and node.annotations.get_value('host') == 'swine' and parent.annotations.get_value('host') == 'human':
            # Calculate the height of the subtree + the edge length to the parent
            assert isinstance(node, Node)
            assert isinstance(tree, Tree)
            height = node.distance_from_tip()
            shortest_height = sorted([leaf.distance_from_root() - node.distance_from_root() for leaf in node.leaf_nodes()])[0]
            # spillover_date1 = node_dates[parent.label]
            spillover_length = height + node.edge_length
            initial_length = shortest_height + node.edge_length
            spillover_date = root_age + parent.distance_from_root()
            spillover_year = int(spillover_date)
            spillover_month = int((spillover_date - spillover_year) * 365 / 31) + 1  # Approximate month estimation.
            spillover_date2 = datetime(year=spillover_year, month=spillover_month, day=1)
            # if spillover_date2.month != spillover_date1.month:
                # print(parent.label, spillover_date1, spillover_date2)
            spillover_season = date_to_season(datetime(year=spillover_year, month=spillover_month, day=1), midyear=True)
            if spillover_length >= threshold:
                long_spillovers_by_season[spillover_season] = long_spillovers_by_season.get(spillover_season, 0) + 1
                if initial_length < threshold:
                    long_sustained_spills_by_season[spillover_season] = long_sustained_spills_by_season.get(spillover_season, 0) + 1
            if spillover_length > longest_spillovers.get(spillover_season, 0):
                longest_spillovers[spillover_season] = spillover_length
            spillovers_by_season[spillover_season] = spillovers_by_season.get(spillover_season, 0) + 1
            total_spill_length[spillover_season] = total_spill_length.get(spillover_season, 0) + spillover_length
        elif parent and node.annotations.get_value('host') == 'human' and parent.annotations.get_value('host') == 'swine':
            # Swine-to-human spillover.
            spillover_date = root_age + node.distance_from_root()
            spillover_year = int(spillover_date)
            spillover_month = int((spillover_date - spillover_year) * 365 / 31) + 1  # Approximate month estimation.
            spillover_season = date_to_season(datetime(year=spillover_year, month=spillover_month, day=1), swine=True)
            zoonotic_spillovers[spillover_season].append(node.leaf_nodes())
            print(f'Zoonotic spill season: {spillover_season}')
            for leaf in node.leaf_nodes():
                print(f'\t{leaf.taxon.label}')

    for i, season in enumerate(seasons):
        total_spills = spillovers_by_season.get(season, 0)
        long_spills = long_spillovers_by_season.get(season, 0)
        long_sustained_spills = long_sustained_spills_by_season.get(season, 0)
        p_long = (long_spills / total_spills) * 100 if total_spills > 0 else 0
        p_sust_long = (long_sustained_spills / total_spills) * 100 if total_spills > 0 else 0
        avg_spill = total_spill_length[season] / total_spills if total_spills > 0 else 0
        longest_spill = longest_spillovers[season] if total_spills > 0 else 0
        spillover_log.write(f'{total_spills},{p_long},{p_sust_long},{longest_spill},{avg_spill}')
        if i < len(seasons) - 1:
            spillover_log.write(',')

        if total_spills > 0:
            p_long = (long_spills / total_spills) * 100
            print('Season %s. Long spillovers: %.2f (%d/%d). Long sustained: %.2f. Average length: %.3f. Longest spillover: %.3f. Zoonotic spills: %d' %
                  (season, p_long, long_spills, total_spills, (long_sustained_spills / total_spills) * 100,
                   total_spill_length[season] / total_spills,
                   longest_spillovers[season], len(zoonotic_spillovers[season])))

    spillover_log.write('\n')
    spillover_log.close()

    with open('data/spillovers.csv', 'w') as spillovers_log:
        spillovers_log.write('season,reverse-zoonotic,zoonotic\n')
        for season in seasons:
            spillovers_log.write(f'{season},{spillovers_by_season.get(season, 0)},{len(zoonotic_spillovers[season])}\n')


def spillovers_by_season(sim_matrix: HuSwSimilarityMatrix, tree_path='../trees/pdm09_US_all.timetree.hosts.tre',
                         log_path='../trees/spillover_contributions.csv', query_seasons=None):
    seasons = ['10-11', '11-12', '12-13', '13-14', '14-15', '15-16', '16-17', '17-18', '18-19', '19-20', '20-21']
    if not query_seasons:
        query_seasons = seasons
    # Open log and write the header.
    contribs_log = open(log_path, 'w')
    for (i, season) in enumerate(seasons):
        if season not in query_seasons:
            continue
        for contrib_season in (['original'] + seasons[:i+1]):
            contribs_log.write(f'{season}_{contrib_season}')
            if not (season == query_seasons[-1] and contrib_season == query_seasons[-1]):
                contribs_log.write(',')
    contribs_log.write('\n')

    season_totals = {}
    season_breakdowns = {}
    for i, season in enumerate(seasons):
        if season not in query_seasons:
            continue
        sw_season_range = season_to_date_range(season, swine=True)
        sw_seqs = [sw_seq for sw_seq in sim_matrix.swine_sequences if
                   sw_season_range[0] <= sw_seq.date <= sw_season_range[1]]
        total_seqs = len(sw_seqs)
        season_totals[season] = total_seqs

        by_season_spills = {}
        future_spills = 0
        for spill_season in reversed(seasons[:i+1]):
            spill_season_seqs = len(find_spillovers_phylogenetic(spill_season, sim_matrix, tree_path=tree_path, swine_season=season))
            by_season_spills[spill_season] = spill_season_seqs - future_spills
            future_spills = spill_season_seqs
        by_season_spills['original'] = total_seqs - future_spills
        # Write the log for later aggregation:
        for contrib_season in ['original'] + seasons[:i+1]:
            contribs_log.write(str(by_season_spills[contrib_season]))
            if not (season == query_seasons[-1] and contrib_season == query_seasons[-1]):
                contribs_log.write(',')

        season_breakdowns[season] = by_season_spills

    contribs_log.write('\n')
    contribs_log.close()

    with open('data/season_spills.csv', 'w+') as output:
        for season in ['total', 'original'] + seasons:
            output.write(',' + season)
        output.write('\n')

        for i, season in enumerate(seasons):
            if season not in query_seasons:
                continue
            output.write(season)
            output.write(',%d' % season_totals[season])
            for spill_season in ['original'] + seasons[:i+1]:
                output.write(',%d' % season_breakdowns[season][spill_season])
            output.write('\n')
