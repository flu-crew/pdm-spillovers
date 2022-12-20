# -*- coding: utf-8 -*-

"""
Read in 'data/spillover-training-recent.fasta' and 'data/spillover-training-continued.fasta'
files and determine key differences between the sequences using a random forest classifier.
"""

import numpy as np
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import precision_score, recall_score
from sklearn import tree
from sklearn.model_selection import cross_val_score, cross_validate
from matplotlib import pyplot
from . import pdm09, helpers, HuSwSimilarityMatrix


def preprocess_and_predict(similarity_matrix: HuSwSimilarityMatrix, exclude_third=True):
    # recent_spills = helpers.read_iav_sequences('data/spillover-training-recent.fasta')
    # continued_spills = helpers.read_iav_sequences('data/spillover-training-continued.fasta')
    recent_spills = helpers.read_iav_sequences('data/spillover-training-all-recent.fasta')
    continued_spills = helpers.read_iav_sequences('data/spillover-training-all-prev.fasta')
    predict_spills = helpers.read_iav_sequences('data/20-21-all.fasta')
    all_spills = [(iav, 1) for iav in recent_spills] + [(iav, 0) for iav in continued_spills] +\
                 [(iav, None) for iav in predict_spills]
    alignment_length = len(recent_spills[0].seq)
    features = []
    feature_names = []
    labels = [label for _, label in all_spills if label is not None]
    no_insertion_ind = 0
    for i in range(alignment_length):
        base_frequency = {'a': 0, 't': 0, 'g': 0, 'c': 0, '-': 0}
        for iav, _ in all_spills:
            nt = iav.seq[i].lower()
            if nt in base_frequency.keys():
                base_frequency[nt] += 1
        most_frequent_base = max(base_frequency.items(), key=lambda x: x[1])[0]
        if most_frequent_base == '-':
            continue  # Skip this site.
        is_consensus = base_frequency[most_frequent_base] > len(all_spills) / 2  # >50% support.

        if exclude_third and (no_insertion_ind + 1) % 3 == 0:
            no_insertion_ind += 1
            continue
        variations = 0  # Number of sequences with a base different from 'most_frequent_base'.
        nts = []
        for iav, label in all_spills:
            nt = iav.seq[i].lower()
            if nt not in {'a', 't', 'g', 'c'}:
                nt = most_frequent_base
                if not is_consensus:
                    print('Site: %d. Unusual nt: %s. IAV: %s. No consensus' % (i, nt, iav.full_name))
            if nt != most_frequent_base:
                variations += 1
            nts.append(nt)
        if variations >= 2:
            # We then treat site i as informative and add it to the feature set.
            base_to_index = {}
            for base in ['a', 't', 'g', 'c']:
                if base_frequency[base] > 0:
                    feature = []
                    feature_name = '%d%s%d' % (i + 1, base, no_insertion_ind + 1)
                    feature_names.append(feature_name)
                    for nt in nts:
                        feature.append(1 if nt == base else 0)
                    features.append(feature)

        no_insertion_ind += 1
    print('Features #: %d' % len(features))
    X = np.array(features)
    X = X.transpose()
    train_X = X[:len(recent_spills) + len(continued_spills)]
    predict_X = X[len(recent_spills) + len(continued_spills):]
    dtc = DecisionTreeClassifier(max_depth=3, criterion='gini', class_weight={0: 1, 1: 1}, random_state=938870)
    # dtc = RandomForestClassifier(n_estimators=2000)
    dtc.fit(train_X, labels)
    for feature_ind, importance in enumerate(dtc.feature_importances_):
        if importance > 0.005:
            print(feature_ind, feature_names[feature_ind], importance)
    print('Precision: ', precision_score(labels, dtc.predict(train_X)))
    print('Recall: ', recall_score(labels, dtc.predict(train_X)))

    pos_predicted = 0
    for i, label in enumerate(dtc.predict(predict_X)):
        if label == 1:
            print(predict_spills[i].full_name)
            pos_predicted += 1
    print('Predicted recent: %d/%d' % (pos_predicted, len(predict_spills)))

    if isinstance(dtc, DecisionTreeClassifier):
        pyplot.figure(figsize=(12, 6))
        tree.plot_tree(dtc, feature_names=feature_names, fontsize=5)
    print(cross_val_score(dtc, train_X, labels, scoring='precision').mean())
    print(cross_val_score(dtc, train_X, labels, scoring='recall').mean())

    pyplot.show()
