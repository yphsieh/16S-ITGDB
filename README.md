# 16S-ITGDB
16S-ITGDB is an integrated database for improving taxonomic classification of 16S ribosomal RNA (rRNA) sequences. Taxonomic classification requires 16S rRNA classifiers and reference databases. Reference databases such as SILVA, RDP, and Greengenes were used for classification tasks. Each database has its own pros and cons, 16S-ITGDB takes the advantages of these databases and were designed for better taxonomy assignment performance.

We proposed two types of integrations, including sequence-based and taxonomy-based ITGDB. Both of them are currently integrated from RDP (version NO18), SILVA (version 138), and Greengenes (version 13_8) OTU clustering databases with 99% of sequence similarity. 

## Sequence-based ITGDB
Sequence-based ITGDB integrates all the unique sequences from RDP, SILVA, and Greengenes without limiting the collected sequences need to have exact species names (no taxonomic depth limitation).

## Taxonomy-based ITGDB
Taxonomy-based ITGDB integrates the sequences based on their taxonomies. All the taxonomies of RDP, SILVA, and Greengenes were taken into consideration and all the collected sequences had exact species names. Taxonomy-based ITGDB was suggested for 16S full-length classification. 

## File format explanation
The sequence-based and taxonomy-based ITGDBs are in ```data/``` directory.
The validation datasets are in ```test cases/``` directory.

Two files are provided for each 16S-ITGDB and validation dataset, including a sequence file and a taxonomy file. Sequence files are in FASTA formats (named as ```*_seq.fasta```), and the extension of taxonomy files are "txt" (named as ```*_taxa.txt```).

In sequence files, each sequence has its seuqence ID and accompanied with its corresponding sequence. For example, in the following picture, ```AJ000684``` and ```EF599163``` are the seuqence IDs and followed by their corresponding sequences.

<img width="1057" alt="Screen Shot 2021-09-03 at 9 46 28 AM" src="https://user-images.githubusercontent.com/47639979/131938695-b98a4a4b-2040-4985-a507-8fcb8ae825df.png">

In taxonomy files, each ID and its corresponding taxonomy is separated by a tab key. The corresponding taxonomies of ```AJ000684``` and ```EF599163``` are shown the following picture. The IDs links the sequence files and taxonomy files together.

<img width="988" alt="Screen Shot 2021-09-03 at 9 47 13 AM" src="https://user-images.githubusercontent.com/47639979/131938704-5d3b3e9f-b637-4358-bc6f-79f2c199c048.png">

## Usage
Taxonomy-based ITGDB is suggested for 16S full-length classification. The following content shows how to use ITGDBs in SINTAX, SPINGO, Mothur, and QIIME2 classifiers. All the classifiers listed here provide cut-off value settings to prevent over-classification issue when usgin NGS short sequences (length = 200 ~ 500 bp) for taxonomy assignment tasks. However, third generation sequencing (TGS) of 16S rRNA reads include V1-V9 regions (full-length is about 1500 bp). Applying default setting of cut-off values will result in very conservative taxonomy assignment results (many sequences could not be assigned to the species level). Therefore, deactivate the cut-off settings is suggested for 16S full-length assignment (but not for NGS short sequence data). Users could refer our 16S-ITGDB paper for detail discussion of taxonomic depth and accuracy of deactivating cut-off settings for 16S full-length sequence assignment.

### SINTAX
To serve as the reference database in SINTAX algorithms, the taxonomy and sequence files should be combined and converted into a ```UDB``` file in the following format:<br/>
<img width="1037" alt="Screen Shot 2021-09-03 at 10 02 52 AM" src="https://user-images.githubusercontent.com/47639979/131939495-8c3a85ec-4c04-473c-82b1-f817c42f97c9.png"><br/>
The converted sequence-based and taxonomy-based integrated database are provided as ```seq_itgdb.udb``` and ```taxa_itgdb.udb``` in ```data/``` directory.<br/>

We assume your SINTAX execution(binary) file is named as "usearch" and use 8 processors(threads) for parallel computation, the following command is used to assign taxonomies:<br/>
```
./usearch -sintax <input file> -db <reference database> -tabbedout <output file> -strand <plus/both> -sintax_cutoff <bootstrap cutoff> -threads <threads number>
```
As previously suggested, the cut-off is deactivated by setting bootstrap cutoff as 0. For example, to assign the taxonomies of the ```Intersection``` dataset, we used:<br/>
```
./usearch -sintax intersect_seq.fasta -db taxa_itgdb.udb -tabbedout sintax_intersect_itgdb_results.tsv -strand both -sintax_cutoff 0 -threads 8
```
Detailed tutorials can be found in: https://www.drive5.com/usearch/manual/cmd_sintax.html.<br/>

### SPINGO
SPINGO requires a species specific database in the following format: <br/>
<img width="1217" alt="Screen Shot 2021-09-22 at 11 01 45 AM" src="https://user-images.githubusercontent.com/47639979/134276818-8b3ef7b7-f1f0-4eaf-95cd-d289b5e4d9ab.png">

We provided ```taxa_itgdb_spingo.fa``` in ```data/``` directory as the species specific database. The reference database is converted internally to an efficient indexed structure. To reuse the index, ```--write-index``` option is used. SPINGO provides parallel computation. Assume your SPINGO execution file is named as "spingo", the utilized command is as follows:
```
./spingo --write-index -p <number of processors> -d <reference database> -i <input file> > <output file>
```
For example, we used the command:
```
./spingo --write-index -p 8 -d database/taxa_itgdb_spingo.fa -i intersect_seq.fasta > spingo_intersect_itgdb_results.tsv
```
to assign the taxonomies of the sequences in the ```Intersection``` dataset with 8 processors for the program. <br/>
Detailed descriptions are in: https://github.com/GuyAllard/SPINGO.<br/>

### Mothur classifier
The classify.seqs command uses reference files to assign the taxonomies of the sequences in your fasta file with designated number of processors. The taxonomy file compatible with mothur classifier is provided as `taxa_itgdb_taxa_mothur.txt` in `data/`.
```
mothur > classify.seqs(fasta=<input file>, processors=<number of processors>, reference=<sequence file of the reference database>, taxonomy=<taxonomy file of the reference database>, methods=<wang>, cutoff=<bootstrap cutoff>)
```
mothur will output two files from the classify.seqs command: a ```*.taxonomy``` file which contains a taxonomy string for each sequence, and a ```*.tax.summary``` file which contains a taxonomic outline indicating the number of sequences that were found for your collection at each level.<br/>
We used the following command to assgin taxonomy for ```Intersection``` dataset using 8 processors and deactivated cut-off setting (cutoff=0):
```
mothur > classify.seqs(fasta=intersect_seq.fasta, processors=8, reference=taxa_itgdb_seq.fasta, taxonomy=taxa_itgdb_taxa_mothur.txt, methods=wang, cutoff=0)
```
Detailed usage can be found in: https://mothur.org/wiki/classify.seqs/.

### QIIME2 classifier
The trained QIIME2 classifiers are in ```data/``` directory. The file names are "taxa_itgdb_q2_2020_08_clf.qza", "taxa_itgdb_q2_2021_08_clf.qza", and  "taxa_itgdb_q2_2022_11_clf.qza". In these files, "2020_08" corresponds to QIIME2 version 2020.08, "2021_08" as QIIME2 version 2021_.08, and "2022_11" indicates QIIME2 version 2022.11. Please download the version that compatible with your QIIME2 pipeline. The usage of taxonomic assignment with desired number of processors is shown below.
```
# import the input file as a QIIME2 artifact
qiime tools import --type 'FeatureData[Sequence]' --input-path <input file> --output-path <input QIIME2 artifact>
# taxonomic assignment
qiime feature-classifier classify-sklearn --i-classifier <trained QIIME2 artifact> --i-reads <input QIIME2 artifact> --o-classification <output QIIME2 artifact> --p-n-jobs <number of processors> --p-confidence <disable or bootstrap cutoff>
```
The output file is a ```*.qza``` file, which can be exported to ```taxonomy.tsv``` by:
```
qiime tools export --input-path <output file> --output-path ./
```
As mentioned before, the cut-off value should be disabled for 16S full-length sequences. For instance, the commands used to classify the sequences in ```Intersection``` dataset with 8 processors are:
```
qiime tools import --type 'FeatureData[Sequence]' --input-path intersect_seq.fasta --output-path intersect_seq.qza

qiime feature-classifier classify-sklearn --i-classifier taxa_itgdb.qza --i-reads intersect_seq.qza --o-classification qiime2_intersect_itgdb_results.qza --p-n-jobs 8 --p-confidence disable

qiime tools export --input-path qiime2_intersect_itgdb_results.qza --output-path ./

mv taxonomy.tsv qiime2_intersect_itgdb_results.tsv
```
Detailed tutorials of QIIME2 usage are in : https://docs.qiime2.org/2020.8/tutorials/. <br/>



## License
Shield: [![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

This datasets of 16S-ITGDB are licensed under a
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa].

[![CC BY-NC-SA 4.0][cc-by-nc-sa-image]][cc-by-nc-sa]

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg

