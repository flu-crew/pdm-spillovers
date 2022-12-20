#!/usr/bin/env python
import sys
import os
import subprocess
from Bio import SeqIO
from datetime import datetime
import re


def extract_dates(path: str, format='%Y-%m-%d') -> str:
    records = SeqIO.parse(path, 'fasta')
    file_name = path + '.dates.csv'
    dates_file = open(file_name, 'w+')
    dates_file.write('name, date\n')
    for record in records:
        name = record.name
        # date_str = name.split('|')[-1]
        for token in name.split('|'):
            if re.fullmatch(r'[\d\-/]{4,}', token):
                if token.count('/') == 2:
                    date = datetime.strptime(token, '%m/%d/%Y')
                elif token.count('/') == 1:
                    date = datetime.strptime(token, '%m/%Y')
                elif token.count('-') == 2:
                    date = datetime.strptime(token, '%Y-%m-%d')
                elif token.count('-') == 1:
                    date = datetime.strptime(token, '%Y-%m')
                else:
                    date = datetime.strptime(token, '%Y')
        dec_date = date.year + ((date.month-1)*30 + date.day)/365.0
        dates_file.write('%s, %.2f\n' % (name, dec_date))
    dates_file.close()
    return file_name


def root_tree(tree_path: str, alignment_path: str, rooted_path=None) -> str:
    dates_file = extract_dates(alignment_path)
    if not rooted_path:
        rooted_path = alignment_path + '.rooted.tre'
    treetime_dir = alignment_path + '.treetime'
    print(' '.join(['treetime', 'clock', '--tree', tree_path,
                     '--aln', alignment_path, '--dates', dates_file,
                     '--outdir', treetime_dir]))
    subprocess.call(['treetime', 'clock', '--tree', tree_path,
                     '--aln', alignment_path, '--dates', dates_file,
                     '--outdir', treetime_dir], stderr=subprocess.STDOUT)
    os.replace(treetime_dir + '/rerooted.newick', rooted_path)
    return rooted_path


if __name__ == '__main__':
    args = sys.argv[1:]
    root_tree(args[0], args[1])
