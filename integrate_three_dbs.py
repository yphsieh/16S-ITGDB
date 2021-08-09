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

def createCSV(gg_txtFile):
	filename = gg_txtFile.split('.')[0]+".csv"
	df = pd.read_csv(gg_txtFile, sep="\t",header=None)
	df.to_csv(filename, header=False, index=False)
	print('finish converting')

def Preprocess(fastaFile='gg_13_5.fasta', txtFile='gg_13_5_taxonomy.txt'):
	if os.path.isfile(txtFile[:-4]+'.csv'): print('CSV file already exist... passing ...')
	else: createCSV(txtFile)

	gg_taxFile = txtFile[:-4]+'.csv'
	gg = open(fastaFile)
	gg_dict = {}
	max_length, min_length, idx = 0, 100000, 1
	
	for line in gg:
		if line.startswith('>'):
			if idx != 1:
				gg_dict[gg_key] = seq
				max_length = max(max_length, len(seq))
				min_length = min(min_length, len(seq))
			gg_key = line[1:-1]
			seq = ''

		else:
			tmp = line.strip('\n')
			seq += tmp.upper()
		idx += 1
	gg_dict[gg_key] = seq
	gg.close()
	print(f'The longest/shortest in {fastaFile}: {max_length} / {min_length}')

	cnt = 0 
	not_found = []
	with open(gg_taxFile, newline='') as csvfile:
		taxa_reader = csv.reader(csvfile)
		for id_name in taxa_reader:
			if id_name[0] in gg_dict.keys():
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
				if tax_name in gg_dict.keys():
					gg_dict[tax_name].append([id_name[0], gg_dict[id_name[0]]])
					del gg_dict[id_name[0]]
				else:
					gg_dict[tax_name] = [[id_name[0], gg_dict[id_name[0]]]]
					del gg_dict[id_name[0]]

				#if cnt % 10000 == 0 : print(tax_name, len(gg_dict[tax_name]))
			else:
				not_found.append(id_name[0])

	cnt = 0
	for i in gg_dict.keys():
		if len(i.split(';')) == 7: cnt += 1
	print(f'{cnt} out of {len(gg_dict.keys())} names are complete in 7 levels.')
	print(f'{len(not_found)} ids are not found in fasta file.')
	return gg_dict

def integrateDB(gg_dict, slv_dict):
	gg_cnt = 0
	slv_cnt = 0
	for i in slv_dict.keys():
		slv_cnt += len(slv_dict[i])
	for i in gg_dict.keys():
		gg_cnt += len(gg_dict[i]) 

	print(f'There are {gg_cnt} in GG, {slv_cnt} in Slv, computing {gg_cnt*slv_cnt}')

	IDB = copy.deepcopy(slv_dict)
	done, addin, replace = 0, 0, 0
	tStart = time.time()
	for gg_key in gg_dict.keys():

		for gg_idx, gg_val in enumerate(gg_dict[gg_key]):
			gg_seq = gg_val[1]

			done += 1
			found = 0
			
			for slv_key in slv_dict.keys():
				for idx, slv_val in enumerate(slv_dict[slv_key]):
					slv_seq = slv_val[1]

					if len(gg_seq) > len(slv_seq):
						index = gg_seq.find(slv_seq)
						if index > -1:
							print("found and replace: ", slv_key)
							found = 1
							IDB[slv_key][idx][1] = gg_seq
							replace += 1
							break
					else:
						index = slv_seq.find(gg_seq)
						if index > -1:
							found = 1

				if found and len(gg_seq) > len(slv_seq): break
			
			if not found:
				addin += 1
				if gg_key in IDB.keys():
					IDB[gg_key].append([gg_val[0], gg_seq])
				else:
					IDB[gg_key] = [[gg_val[0], gg_seq]]
			
			if done % 50 == 0 :
				print(f'Costing time: {(time.time()-tStart)/60} min, {done}/{gg_cnt}, Estimate finish time: {((time.time()-tStart)/60) * gg_cnt/done} min')

	return addin, replace, IDB

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
	parser.add_argument('--name', 		type=str  , default='_RSGIDB')
	parser.add_argument('--ggFasta',	type=str  , default='../gg/99_otus.fasta')
	parser.add_argument('--ggTxt',		type=str  , default='../gg/99_otu_taxonomy.txt')
	parser.add_argument('--slvFasta',	type=str  , default='../slv/silva-138-99-seqs.fasta')
	parser.add_argument('--slvTxt',		type=str  , default='../slv/silva-138-99-tax.txt')
	parser.add_argument('--rdpFasta',	type=str  , default='../rdp/RefOTUs.fasta')
	parser.add_argument('--rdpTxt',		type=str  , default='../rdp/Ref_taxonomy.txt')
	args = parser.parse_args()

	gg_dict = Preprocess(fastaFile=args.ggFasta, txtFile=args.ggTxt)
	slv_dict = Preprocess(fastaFile=args.slvFasta, txtFile=args.slvTxt)
	rdp_dict = Preprocess(fastaFile=args.rdpFasta, txtFile=args.rdpTxt)

	addin_rdp, replace_rdp, RSIDB = integrateDB(slv_dict, rdp_dict)
	addin_gg, replace_gg, IDB = integrateDB(gg_dict, RSIDB)
	
	print(f'add {addin_rdp} from rdp, add {addin_gg} from gg')
	print(f'replace {replace_rdp} from rdp, replace {replace_gg} from gg')
	saveFile(args.name, IDB)
