# "gisaid_H1N1_human_2019-2021-accessions.txt" contains the list of EPI GISAID accessions.
# To obtain the sequences, you can download all H1N1-associated genes (human host) from GISAID
# collected between 2019 and 2021. Please make sure that strain names are as follows:
# Isolate ID | DNA Accession no. | Isolate name | Segment | Collection date
$ Assuming that the downloaded sequences are stored in gisaid-all.fasta:
# run "smof grep -f gisaid_H1N1_human_2019-2021-accessions.txt gisaid-all.fasta > gisaid.fasta"
# to keep only those sequences that match the provided list of accessions.
# Next, run the following commands:


cat gisaid.fasta IRD_H1_swine_2019-2021.fasta > merged_2019-21.fasta
makeblastdb -in merged_2019-21.fasta -title "hu-sw-H1-IAV" -dbtype nucl

blastn -db merged_2019-21.fasta -query ../variants/variant-orfs.fasta -max_target_seqs 10 -outfmt 6 > blast-output_10.txt

# Find results for individual genes as follows:
cat blast-output_10.txt | grep A/North_Carolina/01/2021 | grep PB2
