'''
Progress: no add in!! only replace
Author: Hsieh, Yu-Peng
Last modified: Apr. 18 2020

'''

import pandas as pd
import argparse
import csv
import os
import time
import copy

def createCSV(txtFile):
	filename = txtFile.split('.')[0]+".csv"
	df = pd.read_csv(txtFile, sep="\t",header=None)
	df.to_csv(filename, header=False, index=False)
	print(f'finish converting to {filename}')

def Preprocess(fastaFile, txtFile):
	if os.path.isfile(txtFile[:-4]+'.csv'): print('CSV file already exist... passing ...')
	else: createCSV(txtFile)

	taxFile = txtFile[:-4]+'.csv'
	fastaReader = open(fastaFile)
	taxa_seq_dict = {}
	max_length, min_length, idx = 0, 100000, 1
	
	for line in fastaReader:
		if line.startswith('>'):
			if idx != 1:
				taxa_seq_dict[id_key] = seq
				max_length = max(max_length, len(seq))
				min_length = min(min_length, len(seq))
			id_key = line[1:-1]
			seq = ''

		else:
			tmp = line.strip('\n')
			seq += tmp.upper()
		idx += 1
	taxa_seq_dict[id_key] = seq
	fastaReader.close()
	print(f'The longest/shortest in {fastaFile}: {max_length} / {min_length}')

	cnt = 0 
	not_found = []
	with open(taxFile, newline='') as csvfile:
		taxa_reader = csv.reader(csvfile)
		for id_name in taxa_reader:
			if id_name[0] in taxa_seq_dict.keys():
				cnt += 1

				## cleanse gg taxa names
				tax_name = id_name[1].split(';')
				for i in range(len(tax_name)):
					if tax_name[i].endswith('_'):
						tax_name = tax_name[:i]
						break
					else:
						tax_name[i] = tax_name[i].split('__')[1].lower()

				## combine taxonomy into a string
				tax_name = ';'.join(tax_name)

				## replace id with name
				if tax_name in taxa_seq_dict.keys():
					taxa_seq_dict[tax_name].append([id_name[0], taxa_seq_dict[id_name[0]]])
					del taxa_seq_dict[id_name[0]]
				else:
					taxa_seq_dict[tax_name] = [[id_name[0], taxa_seq_dict[id_name[0]]]]
					del taxa_seq_dict[id_name[0]]

			else:
				not_found.append(id_name[0])

	cnt = 0
	for i in taxa_seq_dict.keys():
		if len(i.split(';')) == 7: cnt += 1
	print(f'{cnt} out of {len(taxa_seq_dict.keys())} names are complete in 7 levels.')
	print(f'{len(not_found)} ids are not found in fasta file.')
	return taxa_seq_dict

def integrateDB(candidate_dict, basis_dict):
	candidate_cnt = 0
	basis_cnt = 0
	for i in basis_dict.keys():
		basis_cnt += len(basis_dict[i])
	for i in candidate_dict.keys():
		candidate_cnt += len(candidate_dict[i]) 

	print(f'There are {candidate_cnt} in candidate, {basis_cnt} in basis, computing {candidate_cnt*basis_cnt}')

	Seq_ItgDB = copy.deepcopy(basis_dict)
	done, addin, replace = 0, 0, 0
	tStart = time.time()
	for candidate_key in candidate_dict.keys():

		for candidate_idx, candidate_val in enumerate(candidate_dict[candidate_key]):
			candidate_seq = candidate_val[1]

			done += 1
			found = 0
			
			for basis_key in basis_dict.keys():
				for idx, basis_val in enumerate(basis_dict[basis_key]):
					basis_seq = basis_val[1]

					if len(candidate_seq) > len(basis_seq):
						index = candidate_seq.find(basis_seq)
						if index > -1:
							print("found and replace: ", basis_key)
							found = 1
							Seq_ItgDB[basis_key][idx][1] = candidate_seq
							replace += 1
							break
					else:
						index = basis_seq.find(candidate_seq)
						if index > -1:
							found = 1

				if found and len(candidate_seq) > len(basis_seq): break
			
			if not found:
				addin += 1
				if candidate_key in Seq_ItgDB.keys():
					Seq_ItgDB[candidate_key].append([candidate_val[0], candidate_seq])
				else:
					Seq_ItgDB[candidate_key] = [[candidate_val[0], candidate_seq]]
			
			if done % 50 == 0 :
				print(f'Costing time: {(time.time()-tStart)/60} min, {done}/{candidate_cnt}, Estimate finish time: {((time.time()-tStart)/60) * candidate_cnt/done} min')

	return addin, replace, Seq_ItgDB

def saveFile(name, Seq_ItgDB):
	Seq_ItgDBFasta = open( name + '.fasta', 'w')
	Seq_ItgDBTxt = open( name + '.txt', 'w')

	for i in Seq_ItgDB.keys():
		for j in Seq_ItgDB[i]:
			ids = '>' + j[0]
			Seq_ItgDBFasta.write(ids+'\n')
			Seq_ItgDBFasta.write(j[1]+'\n')
			Seq_ItgDBTxt.write(ids[1:]+'\t'+i+'\n')

	Seq_ItgDBFasta.close()
	Seq_ItgDBTxt.close()

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('--name', 		type=str  , default='Seq_ItgDB')
	parser.add_argument('--ggFasta',	type=str  , default='../gg/99_otus.fasta')
	parser.add_argument('--ggTaxa',		type=str  , default='../gg/99_otu_taxonomy.txt')
	parser.add_argument('--slvFasta',	type=str  , default='../slv/silva-138-99-seqs.fasta')
	parser.add_argument('--slvTaxa',	type=str  , default='../slv/silva-138-99-tax.txt')
	parser.add_argument('--rdpFasta',	type=str  , default='../rdp/RefOTUs.fasta')
	parser.add_argument('--rdpTaxa',	type=str  , default='../rdp/Ref_taxonomy.txt')
	args = parser.parse_args()

	gg_dict = Preprocess(fastaFile=args.ggFasta, txtFile=args.ggTaxa)
	slv_dict = Preprocess(fastaFile=args.slvFasta, txtFile=args.slvTaxa)
	rdp_dict = Preprocess(fastaFile=args.rdpFasta, txtFile=args.rdpTaxa)

	# Change this order if needed
	addin_slv, replace_slv, RSIDB = integrateDB(slv_dict, rdp_dict)
	addin_gg, replace_gg, Seq_ItgDB = integrateDB(gg_dict, RSIDB)
	
	print(f'add {addin_slv} from SILVA, add {addin_gg} from Greengenes')
	print(f'replace {replace_slv} with SILVA, replace {replace_gg} with Greengenes')
	
	saveFile(args.name, Seq_ItgDB)
