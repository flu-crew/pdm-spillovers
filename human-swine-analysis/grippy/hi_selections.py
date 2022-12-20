# -*- coding: utf-8 -*-
import datetime
import subprocess

from grippy import HuSwSimilarityMatrix, pdm09, helpers

seasons = ['9-10', '10-11', '11-12', '12-13', '13-14', '14-15', '15-16', '16-17', '17-18', '18-19', '19-20']


def get_strains_by_season(matrix: HuSwSimilarityMatrix):
    spills_by_season = {}
    for season in seasons:
        spills_by_season[season] = pdm09.find_spillovers_phylogenetic(season, matrix, swine_season='19-21')

    swine_has_by_season = {season: [] for season in seasons}
    for swine_ha in matrix.swine_sequences:
        if swine_ha.date >= datetime.datetime(year=2020, month=7, day=1):
            spillover_season = 'original'
            for season in reversed(seasons):
                spills = spills_by_season[season]
                if swine_ha.full_name in spills:
                    spillover_season = season
                    break
            swine_has_by_season[spillover_season].append(swine_ha)
    for season in seasons:
        print(f'New swine HAs from {season}: {len(swine_has_by_season[season])}')
        if len(swine_has_by_season[season]) <= 20:
            for ha in swine_has_by_season[season]:
                print('\t', ha.full_name)
        # helpers.write_iav_sequences('data/contemporary_%s.fasta' % season, swine_has_by_season[season])
    return swine_has_by_season


def strain_selection_by_season(matrix: HuSwSimilarityMatrix):
    swine_has_by_season = get_strains_by_season(matrix)
    for season in seasons:
        if len(swine_has_by_season[season]) > 1:
            print(f'------------\nSeason {season}')
            season_path = 'data/contemporary_%s.fasta' % season
            ha1_path = 'data/contemporary_%s_ha1.fasta' % season
            consensus_pre_path = 'data/contemporary_%s_consensus_pre.fasta' % season
            consensus_path = 'data/contemporary_%s_consensus.fasta' % season
            helpers.write_iav_sequences(season_path, swine_has_by_season[season])
            with open(ha1_path, 'w') as ha1_file:
                subprocess.call(['flutile', 'trim', 'ha1', '--subtype', 'H1', '--conversion', 'dna2aa', season_path],
                                stdout=ha1_file, stderr=subprocess.DEVNULL)
            with open(consensus_pre_path, 'w') as consensus_pre_file:
                subprocess.call(['smof', 'consensus', ha1_path], stdout=consensus_pre_file)
            with open(consensus_path, 'w') as consensus_file:
                subprocess.call(['sed', r's/Consensus/A\/swine\/Consensus|2022/g', consensus_pre_path],
                                stdout=consensus_file)

            ha1_seqs = helpers.read_iav_sequences(ha1_path)
            consensus = helpers.read_iav_sequences(consensus_path)[0]
            consensus_sims = []
            for ha1_seq in ha1_seqs:
                if ha1_seq.name.count('A0') > 0:
                    sim = 1 - helpers.aligned_dist(ha1_seq.seq, consensus.seq, normalized=True, ignore_tails=False)
                    consensus_sims.append((sim, ha1_seq))
            consensus_sims.sort(key=lambda x: x[0], reverse=True)
            for sim, ha1_seq in consensus_sims[:20]:
                print(ha1_seq.full_name, sim)
