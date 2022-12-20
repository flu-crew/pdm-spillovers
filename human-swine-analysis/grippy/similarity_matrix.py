# -*- coding: utf-8 -*-
import os
from typing import Tuple, List
from Levenshtein import ratio

from . import helpers
from .sequence import IAVSequence


class SequenceSimilarityMatrix:

    def __init__(self, protein: str, path: str, row_sequences, col_sequences, force_compute=False,
                 aligned=False, ignore_tails=True):
        """
        :param protein: Name of the protein for which the similarity matrix is constructed.
        :param path: Path to the matrix (to load from or to save to)
        :param force_compute: Recompute the matrix even if a matrix was already saved
        :param aligned: The given sequences are aligned
        :param ignore_tails: Only applicable if aligned is True. Ignores missing nt/aa at the tails.
        """
        assert len(row_sequences) > 0 and len(col_sequences) > 0
        self.row_sequences = list(row_sequences)
        self.col_sequences = list(col_sequences)
        self.protein = protein

        if force_compute or not os.path.exists(path):
            processed = 0
            self.matrix = []
            for swine_seq in row_sequences:
                row = []
                for human_seq in col_sequences:
                    # Compute the similarity between sequences and append to row.
                    processed += 1
                    if processed % 10000 == 0:
                        print(processed)
                    if aligned:
                        # Aligned distance
                        sim = 1 - helpers.aligned_dist(str(human_seq.seq), str(swine_seq.seq), normalized=True,
                                                   ignore_tails=ignore_tails)
                    else:
                        # Levenshtein ratio
                        sim = ratio(str(human_seq.seq), str(swine_seq.seq))
                    row.append(sim)
                self.matrix.append(row)
            self._save_matrix(path)
        else:
            self._load_matrix(path)

    def find_n_closest_cols(self, row_ind: int, n=1, col_filter=None, latest_first=True)\
            -> [Tuple[float, IAVSequence]]:
        """
        Returns a list on n closest column sequences in the [(similarity, sequence), ...] format.
        :type col_filter: Callable[[int, IAVSeqence], bool]
        :param col_filter: Takes column index and the respective sequence and returns True
                           if the column needs to be considered
        :param latest_first: If multiple sequences have the same similarity to 'row_ind',
                             the method will sort them from earliest to latest by default.
                             Specify True to reverse that order.
        """
        # For each column that passes col_filter get (similarity, sequence) tuple:
        col_sims = [(sim, self.col_sequences[col_ind]) for col_ind, sim in enumerate(self.matrix[row_ind])
                    if (not col_filter or col_filter(col_ind, self.col_sequences[col_ind]))]
        # Sort tuples by similarity in the descending order:
        if latest_first:
            col_sims = sorted(col_sims, key=lambda x: (x[0], x[1].date), reverse=True)
        else:
            col_sims = sorted(col_sims, key=lambda x: (-x[0], x[1].date), reverse=False)
        return col_sims[:n]

    def get_row_id(self, row_seq: IAVSequence) -> int:
        return self.row_sequences.index(row_seq)

    def get_col_id(self, col_seq: IAVSequence) -> int:
        return self.col_sequences.index(col_seq)

    def get_row_id_by_name(self, name: str) -> int:
        for i, seq in enumerate(self.row_sequences):
            if seq.full_name == name:
                return i
        return -1

    def get_col_id_by_name(self, name: str) -> int:
        for i, seq in enumerate(self.col_sequences):
            if seq.full_name == name:
                return i
        return -1

    def _save_matrix(self, path='clustering/HA-matrix'):
        matrix_file = open(path, 'w+')
        for row in self.matrix:
            matrix_file.write(' '.join(['%.6f' % sim for sim in row]))
            matrix_file.write('\n')
        matrix_file.close()

    def _load_matrix(self, path='clustering/HA-matrix'):
        matrix_file = open(path, 'r')
        self.matrix = []
        for line in matrix_file.readlines():
            if len(line) > 0:
                row = [float(sim) for sim in line.split(' ')]
                self.matrix.append(row)
        matrix_file.close()


class HuSwSimilarityMatrix(SequenceSimilarityMatrix):
    def __init__(self, protein: str, path: str, row_sequences, col_sequences, force_compute=False,
                 aligned=False, ignore_tails=True):
        super().__init__(protein, path, row_sequences, col_sequences, force_compute=force_compute,
                         aligned=aligned, ignore_tails=ignore_tails)

    def get_swine_sequences(self) -> List[IAVSequence]:
        return self.row_sequences

    def get_human_sequences(self) -> List[IAVSequence]:
        return self.col_sequences

    swine_sequences = property(get_swine_sequences)
    human_sequences = property(get_human_sequences)
