# Supplemental Table S2.
# AA diff table between H1N1pdm vaccine strains and swine pdm representative strains

flutile trim ha1 --subtype H1 swine-reps-pdm09.fasta > swine-reps-pdm09.ha1.fasta
flutile trim ha1 --subtype H1 h1pdm_vaccines_HA.fasta > h1pdm_vaccines_HA.ha1.fasta

cat h1pdm_vaccines_HA.ha1.fasta swine-reps-pdm09.ha1.fasta > combined.fasta

flutile aadiff --subtype H1 combined.fasta > aadiff_combined.tab

# Import .tab file to .xlsx
