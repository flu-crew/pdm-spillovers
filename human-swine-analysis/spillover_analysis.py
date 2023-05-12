#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Author: Alexey Markin
from grippy import HuSwSimilarityMatrix, pdm_data, pdm09, hi_selections
from grippy.pdm09 import infer_spillover_seasons, assess_long_spillovers
from helpers import read_saved_sequences
from analysis_info import analysis_name, matrix_file

if __name__ == '__main__':
    # We assume that 'run-analyses-multi.sh' was executed 20 times and the results are stored in "../trees/results/".

    # Read the sequences and the similarity matrix:
    swine_sequences = read_saved_sequences(species='swine', name=analysis_name, aligned=True)
    human_sequences = read_saved_sequences(species='human', name=analysis_name, aligned=True)
    matrix = HuSwSimilarityMatrix('HA', 'data/' + matrix_file, swine_sequences, human_sequences,
                                  force_compute=False, aligned=True, ignore_tails=False)

    # We perform phylogenetic analysis across all 20 replicates and save
    # spillover statistics (and other statistics) to appropriate log files.
    # Later these log files are processed and analyzed using ../trees/log_combine.R script.
    for i in range(1, 21):
        print('==========================')
        path_prefix = f'../trees/results/{analysis_name}-{i}/'

        # Color a time-scaled tree by hosts (swine/human) and save a log with overall hu-to-sw spillover statistics
        # Additionally, save detailed information on every detected variant case (sw-to-hu transmission):
        pdm_data.color_tree_w_hosts(tree_path=path_prefix + 'pdm09_US_all.timetree.hosts.tre',
                                    timetree_path=path_prefix + 'pdm09_US_all_timetree/timetree.nexus',
                                    ha1_path=path_prefix + 'pdm09_US_all.ha1.aln',
                                    ancestral_ha1_path=path_prefix + 'pdm09_US_all_aasub_ancestral/ancestral_sequences.fasta',
                                    spillover_stats=path_prefix + 'spillover_stats.csv',
                                    variants=path_prefix + 'variants.csv',
                                    hosts_confidence_path=path_prefix + 'pdm09_US_all_mugration/confidence.csv'
                                    )

        # Analyze the geography of all human-to-swine spillovers.
        # For each spillover, determine what was the first US state of detection
        # + the total number of US states reached by the spillovers:
        pdm_data.spillover_geo_analysis(tree_path=path_prefix + 'pdm09_US_all.timetree.hosts.tre',
                                        log_path=path_prefix + 'spillovers_geo.csv')

        # Analyze all hu-to-sw spillovers and find how long they persisted in swine, save stats
        # (including the percent of spillovers that were maintained for at least a year):
        pdm09.assess_long_spillovers(tree_path=path_prefix + 'pdm09_US_all.timetree.hosts.tre',
                                     timetree_path=path_prefix + 'pdm09_US_all_timetree/timetree.nexus',
                                     log_path=path_prefix + 'long_spillovers_t1.csv',
                                     threshold=1)

        # Identify the source of swine pdm09 viruses from season 2020-21:
        pdm09.spillovers_by_season(matrix, tree_path=path_prefix + 'pdm09_US_all.timetree.hosts.tre',
                                   log_path=path_prefix + 'spillover_contributions.csv', query_seasons=['20-21'])

        # Rule mining analysis: establish conclusively that 2020-21 swine pdm09 viruses are a result of onward
        # circulation of pdm09 from prior seasons (rather than novel human-to-swine spillovers).
        # This will plot parts of Fig 4, which interrupts the further analysis until plots are closed.
        # To turn it off: set plot_predictors=False.
        infer_spillover_seasons(matrix, tree_path=path_prefix + 'pdm09_US_all.timetree.hosts.tre',
                                training_seasons=('13-14', '14-15', '15-16', '16-17', '17-18', '18-19', '19-20'), phylogenetic=True,
                                predict_season='20-21', svm_reg=75, positive_weight=2, plot_predictors=True, plot=False,
                                print_predicted=True, last_spill_only=False, correct_p=None)

    # Select (contemporary) representative swine pdm09 strains for each season of spillover:
    hi_selections.strain_selection_by_season(matrix)
