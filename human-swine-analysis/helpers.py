# -*- coding: utf-8 -*-
import subprocess
from random import shuffle

from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment

from grippy import IAVSequence
from grippy.helpers import trim_alignment, write_iav_sequences


def align_and_split(path1: str, path2: str, aln_path1: str, aln_path2: str, reference_token=None):
    path1_records_num = len(list(SeqIO.parse(path1, 'fasta')))
    cat_out = subprocess.Popen(['cat', path1, path2], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    with open(aln_path1, 'w+') as aln_file1:
        mafft = subprocess.Popen(['mafft', '/dev/stdin'], stdin=cat_out.stdout, stdout=aln_file1)
        mafft.communicate()
    records = list(AlignIO.read(aln_path1, 'fasta'))
    if reference_token:
        # Prune all sequences as guided by the reference sequence.
        iavs = [IAVSequence(record.id, record.seq) for record in records]
        ref_iav = ref_seq = [iav for iav in iavs if iav.contains_token(reference_token)][0]
        trim_alignment(iavs, ref_iav)
        write_iav_sequences(aln_path1, iavs[:path1_records_num])
        write_iav_sequences(aln_path2, iavs[path1_records_num:])
    else:
        AlignIO.write(MultipleSeqAlignment(records[:path1_records_num]), aln_path1, 'fasta')
        AlignIO.write(MultipleSeqAlignment(records[path1_records_num:]), aln_path2, 'fasta')


def save_processed(sequences: list, species='human', name='HA', aligned=False):
    if aligned:
        seq_file_name = 'data/%s-%s.aln' % (species, name)
    else:
        seq_file_name = 'data/%s-%s.fasta' % (species, name)
    SeqIO.write(sequences, seq_file_name, 'fasta')


def read_saved_sequences(species='human', name='HA', aligned=True):
    if aligned:
        path = 'data/%s-%s.aln' % (species, name)
    else:
        path = 'data/%s-%s.fasta' % (species, name)
    return [IAVSequence(seq.name, seq.seq, host=species) for seq in SeqIO.parse(path, 'fasta')]


def subsample(sequences: list, sample_size=1000, save=False, species='human', protein='HA',
              rnd_shuffle=True) -> list:
    assert species in ['human', 'swine']

    reshuffled_seqs = sequences.copy()
    if rnd_shuffle:
        shuffle(reshuffled_seqs)
    sample = reshuffled_seqs[:sample_size]
    if save:
        save_processed(sample, species=species, name=protein)
    return sample


def process_original_sequences(swine_path='../genes/Swine2010-2020_HA_all-states.fasta',
                               human_path='../genes/Human2010-2020_HA_all-states.fasta',
                               protein='HA', sample_size=1000, rnd_shuffle=True):
    swine_seq = list(SeqIO.parse(swine_path, 'fasta'))
    human_seq = list(SeqIO.parse(human_path, 'fasta'))
    subsample(swine_seq, sample_size=sample_size, save=True, species='swine', protein=protein, rnd_shuffle=rnd_shuffle)
    subsample(human_seq, sample_size=sample_size, save=True, species='human', protein=protein, rnd_shuffle=rnd_shuffle)
