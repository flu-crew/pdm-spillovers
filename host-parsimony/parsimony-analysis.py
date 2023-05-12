#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Alexey Markin
"""
from Bio import SeqIO
import subprocess
from dendropy import Tree, StandardCharacterMatrix, new_standard_state_alphabet
from dendropy.model import parsimony
import sys
import random as rnd


def save_hosts_as_fasta(hosts_path: str):
    host_map = {}
    with open(hosts_path) as hosts_file:
        hosts_file.readline()  # Skip the first line.
        for line in hosts_file:
            line = line.strip('\n').strip()
            if line:
                taxon, host = line.split(', ')
                host = host.lower()[:1]
                host_map[taxon] = host
    
    fasta_path = hosts_path + '.fasta'
    with open(fasta_path, 'w') as fasta_file:
        for taxon, host in host_map.items():
            fasta_file.write(f'>{taxon}\n{host}\n')
    return fasta_path


def compute_host_parsimony(tree_path: str, hosts_path: str):
    tree: Tree = Tree.get(path=tree_path, schema='newick', preserve_underscores=True)
    fasta_path = save_hosts_as_fasta(hosts_path)
    taxon_characters: StandardCharacterMatrix = StandardCharacterMatrix.get_from_path(
        fasta_path, schema='fasta',
        taxon_namespace=tree.taxon_namespace,
        default_state_alphabet=new_standard_state_alphabet('hs')
    )
    taxon_states = taxon_characters.taxon_state_sets_map(gaps_as_missing=False)
    p_score = parsimony.fitch_down_pass(tree.postorder_node_iter(), taxon_state_sets_map=taxon_states)
    print('Minimum # of host transitions:', p_score)

    hu_to_sw = 0
    node: Node
    for node in tree.preorder_node_iter():
        parent: Node = node.parent_node
        site = 0
        colored = False
        if parent:
            parent_state = parent.state_sets[site]
            # print(parent_state)
            if parent_state in node.state_sets[site]:
                node.state_sets[site] = parent_state
                if node.state_sets[site] == 1:
                    node.annotations.add_new('!color', '#FFA500')
                continue
            elif parent_state == 0:
                hu_to_sw += 1
            else:
                for leaf in node.leaf_nodes():
                    leaf.taxon.annotations.add_new('!color', '#0000FF')
                # print(len(node.leaf_nodes()))
        # choose a random state and assign
        state_sets = list(node.state_sets[site])
        if len(state_sets) > 1:
            print('Randomly resolving between hu/sw')
        rnd_state = rnd.choice(state_sets)
        node.state_sets[site] = rnd_state
        if node.state_sets[site] == 1:
            node.annotations.add_new('!color', '#FFA500')
    print('Minimum number of human-to-swine transitions:', hu_to_sw)
    tree.write(path=tree_path + '.parsimony.color.tre', schema='nexus')


if __name__ == '__main__':
    for i in range(1, 21):
        print(f'------------Replicate {i}--------------')
        compute_host_parsimony(f'../trees/results/PDM-HA-USA-09-21-Mar2022-merged-cds-{i}/iqtree.treefile',
                               '../trees/results/PDM-HA-USA-09-21-Mar2022-merged-cds.hosts.csv')
