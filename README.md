# 16S-ITGDB
16S-ITGDB is an integrated database for improving taxonomic classification of 16S ribosomal RNA (rRNA) sequences. Taxonomic classification requires 16S rRNA classifiers and reference databases. Reference databases such as SILVA, RDP, and Greengenes were used for classification tasks. Each database has its own pros and cons, 16S-ITGDB takes the advantages of these databases and were designed for better taxonomy assignment performance.

We proposed two types of integrations, including sequence-based and taxonomy-based ITGDB. Both of them are currently integrated from RDP (version NO18), SILVA (version 138), and Greengenes (version 13_8) OTU clustering databases with 99% of sequence similarity. 

Seq_ItgDB and Taxa_ItgDB can both be regenerated when either RDP, SILVA, or Greengenes is updated.

## Sequence-based ITGDB
With the rapid advancement of third-generation sequencing (TGS), 16S full-length provides opportunities for us to access deeper taxonomic levels (e.g., species level). Sequence-based ITGDB was developed by integrating all the unique sequences of RDP, SILVA, and Greengenes.

## Taxonomy-based ITGDB
There are some shortcomings of the databases. Greengenes has fewer sequences with species information since it has not been updated since 2013. Additionally, both of Greengenes and SILVA contain one fifth anomalies, while RDP has fewer anomalies. To resolve these issues, the sequences with no specific species information are firstly cleansed out from the databases. Then, Taxonomy-based ITGDB is integrated from the cleansed databases.

## Usage
The sequence-based and taxonomy-based ITGDBs are provided in ```data/```.
There are five datasets for validation provided in ```validation data/```.

Two files are provided for each 16S-ITGDB and validation dataset, including a sequence file and a taxonomy file. Sequence files are provided in FASTA formats and named with ```*_seq.txt```, and taxonomy files are provided in TXT formats and named with ```*_taxa.txt```.

In each sequence file, each sequence follows its seuqence ID. For example, in the following picture, ```AJ000684``` and ```EF599163``` are the seuqence IDs and followed by their corresponding sequences.

<img width="1057" alt="Screen Shot 2021-09-03 at 9 46 28 AM" src="https://user-images.githubusercontent.com/47639979/131938695-b98a4a4b-2040-4985-a507-8fcb8ae825df.png">

In each taxonomy file, each ID and its corresponding taxonomy are separated with a tabular key. We can indicate the taxonomies of ```AJ000684``` and ```EF599163``` from the following picture.

<img width="988" alt="Screen Shot 2021-09-03 at 9 47 13 AM" src="https://user-images.githubusercontent.com/47639979/131938704-5d3b3e9f-b637-4358-bc6f-79f2c199c048.png">

### QIIME2, SPINGO, and mothur
These files are all compatilble with QIIME2, SPINGO and mothur. Therefore, we can directly serve them as the input to QIIME2, SPINGO, and mothur.

### SINTAX
SINTAX algorithms requires a single FASTA file containing both taxonomy and sequence information. By simply combining the sequence and taxonomy files into a sintax compatible file, we can input the FASTA file into SINTAX. The required file is as below, where each ID is followed by ```;tax=``` indicating the corresponding taxonomy.

<img width="1037" alt="Screen Shot 2021-09-03 at 10 02 52 AM" src="https://user-images.githubusercontent.com/47639979/131939495-8c3a85ec-4c04-473c-82b1-f817c42f97c9.png">
