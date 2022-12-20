# Download all HA (H1N1) sequences collected between 2009 and 2021 from humans in the US from the GISAID database.
# The strain names should be:  Isolate ID | DNA Accession no. | Isolate name | Segment | Collection date
# Assuming that the downloaded file is called gisaid.fasta, then run

smof grep -f Human2009-2021_H1N1_USA_pdms_epi.txt gisaid.fasta > Human2009-2022_H1N1_USA_pdms-filtered.fasta
