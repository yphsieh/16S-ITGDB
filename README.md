# 16S-ItgDB
This is an integrated database for improving taxonomic classification of 16S Ribosomal RNA sequences. Taxonomic classification requires 16S rRNA classifiers and reference databases. Several reference databases have been proposed, e.g., the commonly used SILVA, RDP, and Greengenes. Each database has its own pros and cons, and therefore, 16S-ItgDB is integrated with an aim to combine the advantages of different databases.

We proposed two types of integrations, including sequence-based and taxonomy-based ITGDB. Both of them are currently integrated from RDP (version NO18), SILVA (version 138), and Greengenes (version 13_8) OTU clustering databases with 99% of sequence similarity. 

## Sequence-based ITGDB
With the boom of third-generation sequencing (TGS), seuqence lengths have been extended, redering higher opportunity to classify deeper taxonomic levels (e.g., at species level). Seq_ItgDB is generated based on this idea by comparing the lengths of sequences with the same taxonomy and keeping the longer sequence.

## Taxonomy-based ITGDB
There are some shortcomings of the databases. Greengenes has fewer sequences with species information since it has not been updated since 2015. Additionally, both of Greengenes and SILVA contain one fifth anomalies, while RDP has fewer anomalies. To resolve these issues, the sequences with no specific species information are firstly cleansed out from the databases. Then, Taxa_ItgDB is integrated from the cleansed databases.

Seq_ItgDB and Taxa_ItgDB can both be regenerated when either SILVA, RDP, or Greengenes is updated.

The sequence-based and taxonomy-based ItgDBs are provided in ```data/```.
There are five datasets for validation provided in ```validation data/```.

