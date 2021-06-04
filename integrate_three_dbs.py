'''
Progress: reduce "containing" sequences
Author: Hsieh, Yu-Peng
Last modified: Apr. 08 2020

'''

import pandas as pd
import argparse
import csv
import os
import time
import copy

def createCSV(db_txtFile):
	filename = db_txtFile.split('.')[0]+".csv"
	df = pd.read_csv(db_txtFile, sep="\t",header=None)
	df.to_csv(filename, header=False, index=False)
	print('finish converting')

def Preprocess(fastaFile, txtFile, isModify = 0):
	if os.path.isfile(txtFile[:-4]+'.csv'): print('CSV file already exist... passing ...')
	else: createCSV(txtFile)

	db_taxFile = txtFile[:-4]+'.csv'
	db = open(fastaFile)
	db_dict = {}
	max_length, min_length, idx = 0, 100000, 1
	
	for line in db:
		if line.startswith('>'):
			if idx != 1:
				db_dict[db_key] = seq
				max_length = max(max_length, len(seq))
				min_length = min(min_length, len(seq))
			db_key = line[1:-1]
			seq = ''

		else:
			tmp = line.strip('\n')
			seq += tmp.upper()
		idx += 1
	db.close()
	print(f'The longest/shortest in {fastaFile}: {max_length} / {min_length}')

	cnt = 0 
	not_found = []
	with open(db_taxFile, newline='') as csvfile:
		taxa_reader = csv.reader(csvfile)
		for id_name in taxa_reader:
			if id_name[0] in db_dict.keys():
				cnt += 1

				tax_name = id_name[1].split(';')

				## only preserve taxonomy with 7 levels for modify files
				if isModify and len(tax_name) < 7: continue

				## cleanse db taxa names
				for i in range(len(tax_name)):
					if tax_name[i].endswith('_'):
						tax_name = tax_name[:i]
						break
					else:
						tax_name[i] = tax_name[i].split('__')[1].lower()

				## combine taxonomy into a string
				tax_name = ';'.join(tax_name)

				## replace id with name
				if tax_name in db_dict.keys():
					db_dict[tax_name].append([id_name[0], db_dict[id_name[0]]])
					del db_dict[id_name[0]]
				else:
					db_dict[tax_name] = [[id_name[0], db_dict[id_name[0]]]]
					del db_dict[id_name[0]]

				if cnt % 10000 == 0 : print(tax_name, len(db_dict[tax_name]))
			else:
				not_found.append(id_name[0])

	cnt = 0
	for i in db_dict.keys():
		if len(i.split(';')) == 7: cnt += 1
	print(f'{cnt} out of {len(db_dict.keys())} names are complete in 7 levels.')
	print(f'{len(not_found)} ids are not found in fasta file.')
	return db_dict

def integrateDB(modify_dict, base_dict):
	modify_cnt = 0
	base_cnt = 0
	for i in base_dict.keys():
		base_cnt += len(base_dict[i])
	for i in modify_dict.keys():
		modify_cnt += len(modify_dict[i]) 

	print(f'There are {modify_cnt} in GG, {base_cnt} in Slv, computing {modify_cnt*base_cnt}')

	IDB = copy.deepcopy(base_dict)
	done = 0
	tStart = time.time()
	for modify_key in modify_dict.keys():

		for modify_idx, modify_val in enumerate(modify_dict[modify_key]):
			modify_seq = modify_val[1]

			done += 1
			found = 0
			
			for base_key in base_dict.keys():
				for idx, base_val in enumerate(base_dict[base_key]):
					base_seq = base_val[1]

					if len(modify_seq) > len(base_seq):
						index = modify_seq.find(base_seq)
						if index > -1:
							print("found and replace: ", base_key)
							found = 1
							IDB[base_key][idx][1] = modify_seq
							break
					else:
						index = base_seq.find(modify_seq)
						if index > -1:
							found = 1

				if found and len(modify_seq) > len(base_seq): break

			if not found:
				if modify_key in IDB.keys():
					IDB[modify_key].append([modify_val[0], modify_seq])
				else:
					IDB[modify_key] = [[modify_val[0], modify_seq]]
			
			if done % 50 == 0 :
				print(f'Costing time: {(time.time()-tStart)/60} min, {done}/{modify_cnt}, Estimate finish time: {((time.time()-tStart)/60) * modify_cnt/done} min')

	return IDB

def saveFile(name, IDB):
	IDBFasta = open( name + '.fasta', 'w')
	IDBTxt = open( name + '.txt', 'w')

	for i in IDB.keys():
		for j in IDB[i]:
			ids = '>' + j[0]
			IDBFasta.write(ids+'\n')
			IDBFasta.write(j[1]+'\n')
			IDBTxt.write(ids[1:]+'\t'+i+'\n')

	IDBFasta.close()
	IDBTxt.close()

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('--name', 		type=str  , default='RSGIDB')
	parser.add_argument('--ggFasta',	type=str  , default='../gg/99_otus.fasta')
	parser.add_argument('--ggTxt',		type=str  , default='../gg/99_otu_taxonomy.txt')
	parser.add_argument('--slvFasta',	type=str  , default='../slv/silva_132_99_16S.fna')
	parser.add_argument('--slvTxt',		type=str  , default='../slv/taxonomy_all_levels.txt')
	parser.add_argument('--rdpFasta',	type=str  , default='../rdp/RefOTUs.fasta')
	parser.add_argument('--rdpTxt',		type=str  , default='../rdp/Ref_taxonomy.txt')
	args = parser.parse_args()

	rdp_dict = Preprocess(fastaFile=args.rdpFasta, txtFile=args.rdpTxt, isModify = 0)
	gg_dict  = Preprocess(fastaFile=args.ggFasta,  txtFile=args.ggTxt,  isModify = 1)
	slv_dict = Preprocess(fastaFile=args.slvFasta, txtFile=args.slvTxt, isModify = 1)

	RSIDB = integrateDB(slv_dict, rdp_dict)
	RSGIDB = integrateDB(gg_dict, RSIDB)
	saveFile(args.name, RSGIDB)
