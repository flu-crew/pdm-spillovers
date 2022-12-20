# -*- coding: utf-8 -*-
import subprocess
from typing import Tuple, List, Optional
from dendropy import Tree
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment, AlignInfo

import os

from grippy import IAVSequence
from .treetime_root import root_tree

FASTTREE_PATH = 'FastTree'


def read_iav_sequences(path: str, protein='HA', host='swine') -> [IAVSequence]:
    return [IAVSequence(seq.name, seq.seq, protein=protein, host=host)
            for seq in SeqIO.parse(path, 'fasta')]


def write_iav_sequences(path: str, sequences: [IAVSequence]):
    bio_seqs = [seq.to_bio() for seq in sequences]
    SeqIO.write(bio_seqs, path, 'fasta')


def sort_by_date(sequences: [IAVSequence]) -> [IAVSequence]:
    return sorted(sequences, key=lambda seq: seq.date)


def get_consensus_seq(sequences: List[IAVSequence], name='Consensus') -> IAVSequence:
    max_date = max([iav.date for iav in sequences])
    min_date = min([iav.date for iav in sequences])
    bio_alignment = MultipleSeqAlignment([seq.to_bio() for seq in sequences])
    alignment_summary = AlignInfo.SummaryInfo(bio_alignment)
    consensus_seq = alignment_summary.gap_consensus(threshold=0.5, ambiguous='X')
    defline = '%s/%s/%s|%s' % (name, min_date.strftime('%Y-%m'), max_date.strftime('%Y-%m'),
                                      max_date.strftime('%Y-%m-%d'))
    iav_consensus = IAVSequence(defline, consensus_seq)
    return iav_consensus


def join_genes_by_name(sequences_by_gene: List[List[IAVSequence]]):
    """
    Takes two or more IAVSequence lists and looks at the intersection of those lists based on sequences names.
    Returns the updated lists that contain only common sequence names.
    """
    intersected_names = set()
    for gene_seqs in sequences_by_gene:
        unique_names = {seq.name for seq in gene_seqs}
        if intersected_names:
            intersected_names = intersected_names.intersection(unique_names)
        else:
            intersected_names = unique_names

    joint_sequence_lists = []
    for gene_seqs in sequences_by_gene:
        common_seqs = [seq for seq in gene_seqs if seq.name in intersected_names]
        joint_sequence_lists.append(common_seqs)
    return joint_sequence_lists


def remove_name_duplicates(sequences: List[IAVSequence]) -> List[IAVSequence]:
    filtered_sequences = []
    seen_names = set()
    for seq in sequences:
        if seq.name not in seen_names:
            filtered_sequences.append(seq)
        seen_names.add(seq.name)
    return filtered_sequences


def find_seq_with_name(sequences: List[IAVSequence], name: str) -> Optional[IAVSequence]:
    seqs_w_name = [seq for seq in sequences if seq.name == name]
    if not seqs_w_name:
        return None
    else:
        return seqs_w_name[0]


def build_rooted_tree(sequences: List[IAVSequence], alignment_path: str, force_compute=False,
                      schema='newick') -> Tree:
    extension = '.tre' if schema == 'newick' else '.nexus'
    unrooted_tree = '.'.join(alignment_path.split('.')[:-1]) + extension
    rooted_path = '.'.join(alignment_path.split('.')[:-1]) + '.rooted.%s' % extension
    if not os.path.exists(rooted_path) or force_compute:
        subprocess.call([FASTTREE_PATH, '-nt', '-gtr', '-out', unrooted_tree, alignment_path],
                        stderr=subprocess.STDOUT)
        rooted_path = root_tree(unrooted_tree, alignment_path, rooted_path=rooted_path)
    rooted_tree = get_annotated_tree(rooted_path, sequences, schema=schema)
    return rooted_tree


def get_annotated_tree(tree_path: str, leaf_sequences: List[IAVSequence], schema='newick'):
    rooted_tree = Tree.get(path=tree_path, schema=schema, preserve_underscores=True)
    # Augment leaves' attributes with 'iav_seq':
    sequence_map = {seq.full_name: seq for seq in leaf_sequences}
    for leaf in rooted_tree.leaf_node_iter():
        leaf.iav_seq = sequence_map[leaf.taxon.label.replace(' ', '_')]
    return rooted_tree


def find_closest_sequence(query_seq: IAVSequence, reference_sequences: [IAVSequence],
                          filter=None, ignore_tails=True) -> [Tuple[int, IAVSequence]]:
    if filter:
        reference_sequences = [seq for seq in reference_sequences if filter(seq)]
    similarities = [(1 - aligned_dist(query_seq.seq, ref_seq.seq, ignore_tails=ignore_tails, normalized=True), ref_seq)
                    for ref_seq in reference_sequences]
    sim, closest_ref = max(similarities, key=lambda x: x[0])
    return sim, closest_ref


def trim_alignment(iav_sequences: List[IAVSequence],
                   reference_seq=None, start_pos=None, end_pos=None) -> List[IAVSequence]:
    """
    Given a list of aligned sequences, trims them down to a desired interval.
    :param iav_sequences: Aligned sequences to be trimmed.
    :param reference_seq: If given, the start position (non '-') and the end position of that
                          sequence will be used as an interval
    :param start_pos: First index for the desired interval.
    :param end_pos: Last index for the desired interval.
    :return:
    """
    assert reference_seq or (start_pos and end_pos)
    if reference_seq:
        assert isinstance(reference_seq, IAVSequence)
        start_pos = reference_seq.seq_start()
        end_pos = reference_seq.seq_end()

    for iav_seq in iav_sequences:
        iav_seq.seq = iav_seq.seq[start_pos:end_pos + 1]
    return iav_sequences


def aligned_dist(seq1: str, seq2: str, normalized=True, ignore_tails=True):
    """
    Computes the hamming distance between two aligned sequences
    :param normalized: if True, the distance will be in the [0,1] range.
                       Otherwise, the number of substitutions will be returned.
    :param ignore_tails: if True, the flanking '---' parts will be ignored on both sequences.
    """
    # Assert that 1. sequences have the same lengths and
    # 2. both sequences contain at least one informative (not '-') symbol.
    assert len(seq1) == len(seq2)
    assert seq1.count('-') != len(seq1)
    assert seq2.count('-') != len(seq2)
    length = len(seq1)

    # Compute the number of differing sites.
    difference = 0
    for site, symbol1 in enumerate(seq1):
        symbol2 = seq2[site]
        if symbol1 != symbol2:
            difference += 1

    # Correct the count based on missing tails.
    if ignore_tails:
        seq1_start, seq2_start = 0, 0
        seq1_end, seq2_end = length - 1, length - 1
        while seq1[seq1_start] == '-':
            seq1_start += 1
        while seq2[seq2_start] == '-':
            seq2_start += 1
        while seq1[seq1_end] == '-':
            seq1_end -= 1
        while seq2[seq2_end] == '-':
            seq2_end -= 1
        for i in range(max(seq1_start, seq2_start)):
            if seq1[i] != seq2[i]:
                difference -= 1
        for i in range(length - 1, min(seq1_end, seq2_end), -1):
            if seq1[i] != seq2[i]:
                difference -= 1
        length = min(seq1_end, seq2_end) - max(seq1_start, seq2_start) + 1
    if normalized:
        return difference / length
    else:
        return difference


def append_token_to_file_name(file_name: str, token: str) -> str:
    if file_name.count('.') > 0:
        augmented_name = '.'.join(file_name.split('.')[:-1]) + '.' + token + '.' + file_name.split('.')[-1]
    else:
        augmented_name = file_name + '.' + token
    return augmented_name


def get_strain_name_list(file_name: str) -> List[str]:
    list_handle = open(file_name)
    strain_list = [line[:-1] for line in list_handle.readlines()]  # Remove '\n' from lines.
    list_handle.close()
    return strain_list


def filter_sequences(fasta_path: str, filter_path: str, out_path: str, keep=True, full_token=False):
    iavs = read_iav_sequences(fasta_path)
    filtered_iavs = filter_iavs(iavs, filter_path, keep, full_token)
    write_iav_sequences(out_path, filtered_iavs)


def filter_iavs(iavs: List[IAVSequence], filter_path: str, keep=True, full_token=False) -> List[IAVSequence]:
    filtered_iavs = []
    filter_tokens = get_strain_name_list(filter_path)
    for iav in iavs:
        if full_token:
            in_filter = any([iav.contains_token(filter_token) for filter_token in filter_tokens])
        else:
            in_filter = any([iav.full_name.count(filter_token) > 0 for filter_token in filter_tokens])
        if (keep and in_filter) or (not keep and not in_filter):
            filtered_iavs.append(iav)
        else:
            print(f"Removing {iav.full_name}")
    return filtered_iavs


def filter_tree(tree_path: str, schema: str, filter_path: str, out_path: str, keep=True, full=False):
    filter_tokens = get_strain_name_list(filter_path)
    tree = Tree.get(path=tree_path, schema=schema, preserve_underscores=True)
    assert isinstance(tree, Tree)
    ending = '|' if full else ''
    filter_taxa = {taxon.label for taxon in tree.taxon_namespace if
                   any(taxon.label.count(prune_token + ending) > 0 for prune_token in filter_tokens)}
    if keep:
        prune_taxa = {taxon.label for taxon in tree.taxon_namespace}
        prune_taxa = prune_taxa.symmetric_difference(filter_taxa)
    else:
        prune_taxa = filter_taxa
    tree.prune_taxa_with_labels(prune_taxa)
    tree.write_to_path(out_path, schema=schema)
