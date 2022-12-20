#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from datetime import datetime

from grippy import HuSwSimilarityMatrix, pdm09
from helpers import process_original_sequences, align_and_split, read_saved_sequences
from analysis_info import analysis_name, matrix_file


def build_similarity_matrix() -> HuSwSimilarityMatrix:
    # Aligns human+swine sequences and builds a similarity matrix between them.
    process_original_sequences(swine_path='../pdm09-data/Swine2009_2021-10_pdm_merged.fasta',
                               human_path='../pdm09-data/gisaid/Feb2022/Human2009-2022_H1N1_USA_pdms-filtered.fasta',
                               protein=analysis_name, sample_size=15000, rnd_shuffle=False)
    align_and_split('data/%s-%s.fasta' % ('swine', analysis_name), 'data/%s-%s.fasta' % ('human', analysis_name),
                    'data/%s-%s.aln' % ('swine', analysis_name), 'data/%s-%s.aln' % ('human', analysis_name),
                    reference_token='A/Texas/11913/2019')
    swine_sequences = read_saved_sequences(species='swine', name=analysis_name, aligned=True)
    human_sequences = read_saved_sequences(species='human', name=analysis_name, aligned=True)
    husw_matrix = HuSwSimilarityMatrix('HA', 'data/' + matrix_file, swine_sequences, human_sequences,
                                       force_compute=True, aligned=True, ignore_tails=False)
    return husw_matrix


def make_hosts_csv(husw_matrix: HuSwSimilarityMatrix):
    csv_output = open(f'data/{analysis_name}.hosts.csv', 'w+')
    csv_output.write('name, host\n')
    for iav in husw_matrix.swine_sequences + husw_matrix.human_sequences:
        csv_output.write('%s, %s\n' % (iav.full_name, iav.host))
    csv_output.close()


if __name__ == '__main__':
    husw_matrix = build_similarity_matrix()
    make_hosts_csv(husw_matrix)
