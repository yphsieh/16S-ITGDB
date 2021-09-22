# 16S-ITGDB
16S-ITGDB is an integrated database for improving taxonomic classification of 16S ribosomal RNA (rRNA) sequences. Taxonomic classification requires 16S rRNA classifiers and reference databases. Reference databases such as SILVA, RDP, and Greengenes were used for classification tasks. Each database has its own pros and cons, 16S-ITGDB takes the advantages of these databases and were designed for better taxonomy assignment performance.

We proposed two types of integrations, including sequence-based and taxonomy-based ITGDB. Both of them are currently integrated from RDP (version NO18), SILVA (version 138), and Greengenes (version 13_8) OTU clustering databases with 99% of sequence similarity. 

Seq_ItgDB and Taxa_ItgDB can both be regenerated when either RDP, SILVA, or Greengenes is updated.

## Sequence-based ITGDB
Sequence-based ITGDB integrates all the unique sequences from RDP, SILVA, and Greengenes without limiting the collected sequences need to have exact species names (no taxonomic depth limitation).

## Taxonomy-based ITGDB
Taxonomy-based ITGDB integrates the sequences based on their taxonomies. All the taxonomies of RDP, SILVA, and Greengenes were taken into consideration and all the collected sequences had exact species names. Taxonomy-based ITGDB was suggested for 16S full-length classification. 

## File format explanation
The sequence-based and taxonomy-based ITGDBs are in ```data/``` directory.
The validation datasets are in ```validation data/``` directory.

Two files are provided for each 16S-ITGDB and validation dataset, including a sequence file and a taxonomy file. Sequence files are in FASTA formats (named as ```*_seq.fasta```), and the extension of taxonomy files are "txt" (named as ```*_taxa.txt```).

In sequence files, each sequence has its seuqence ID and accompanied with its corresponding sequence. For example, in the following picture, ```AJ000684``` and ```EF599163``` are the seuqence IDs and followed by their corresponding sequences.

<img width="1057" alt="Screen Shot 2021-09-03 at 9 46 28 AM" src="https://user-images.githubusercontent.com/47639979/131938695-b98a4a4b-2040-4985-a507-8fcb8ae825df.png">

In taxonomy files, each ID and its corresponding taxonomy is separated by a tab key. The corresponding taxonomies of ```AJ000684``` and ```EF599163``` are shown the following picture. The IDs links the sequence files and taxonomy files together.

<img width="988" alt="Screen Shot 2021-09-03 at 9 47 13 AM" src="https://user-images.githubusercontent.com/47639979/131938704-5d3b3e9f-b637-4358-bc6f-79f2c199c048.png">

## Usage
Taxonomy-based ITGDB is suggested for 16S full-length classification. The following content shows how to use ITGDBs in SINTAX, SPINGO, Mothur, and QIIME2 classifiers.

### SINTAX
To serve as the reference database in SINTAX algorithms, the taxonomy and sequence files should be combined and converted into a ```UDB``` file in the following format:<br/>
<img width="1037" alt="Screen Shot 2021-09-03 at 10 02 52 AM" src="https://user-images.githubusercontent.com/47639979/131939495-8c3a85ec-4c04-473c-82b1-f817c42f97c9.png"><br/>
The converted taxonomy-based integrated database is provided as ```taxa_itgdb.udb``` in ```data/``` directory.<br/>

Then, we used the following command to assign taxonomies:<br/>
```./usearch11.0.667_i86linux32 -sintax <input file> -db <reference database> -tabbedout <output file> -strand <plus/both> -sintax_cutoff <bootstrap cutoff>```<br/>
For example, to assign the taxonomies of the ```Intersection``` dataset, we used:<br/>
```./usearch11.0.667_i86linux32 -sintax ./intersect.fasta -db taxa_itgdb.udb -tabbedout taxa_itgdb_taxonomy_intersect.sintax -strand both -sintax_cutoff 0.8```<br/>
Detailed tutorials can be found in: https://www.drive5.com/usearch/manual/cmd_sintax.html.<br/>


### SPINGO
SPINGO requires a species specific database. We provided ```taxa_itgdb_spingo.fa``` in ```data/``` directory as the species specific database. The reference database is converted internally to an efficient indexed structure. To reuse the index, ```--write-index``` option is used. For instance, in this study, we used
```./spingo --write-index -p 48 -d database/taxa_itgdb_spingo.fa -i intersect_seq.fasta > spingo_intersect_itgdb_results.tsv```
to assign the taxonomies of the ```Intersection``` dataset.<br/>

### Mothur classifier
[[[Mothur tutorial]]]<br/>
[[[Upload the organized Mothur files]]]<br/>
[[[Give command pattern and usage example]]]

### QIIME2 classifier
The trained QIIME2 artifacts are in ```data``` directory (seq_itgdb.qza and taxa_itgdb.qza). These artifacts files are trained by QIIME2 version 2020.8, which means these ITGDB artifacts are compatible with QIIME2 version higher than 2020.8. The usage is shown below.[[[tutorial]]]<br/>

[[[Upload *.qza files]]]<br/>
[[[Give command pattern and usage example]]]

### Integrate the newly released RDP, SILVA, and Greengenes
Here shows how to use the source code to integrate the newly released RDP, SILVA and Greengenes. [[[tutorial]]]

