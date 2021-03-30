'''
Progress: reduce "containing" sequences
Author: Hsieh, Yu-Peng
Last modified: Aug. 11 2020

--- Speed up ---
found then break

'''

import pandas as pd
import argparse
import csv
import os
import time

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
	for line in gg:
		if line.startswith('>'):
			gg_key = line[1:].replace('\n','')
		else:
			gg_dict[gg_key] = [line.replace('\n','')]
	gg.close()

	cnt = 0 
	not_found = []
	with open(gg_taxFile, newline='') as csvfile:
		taxa_reader = csv.reader(csvfile)
		for id_name in taxa_reader:
			if id_name[0] in gg_dict.keys():
				cnt += 1

				## cleanse gg taxa names
				tax_name = id_name[1].split(';')
				while tax_name[-1].endswith('_'): tax_name = tax_name[:-1]
				for i in range(len(tax_name)): tax_name[i] = tax_name[i].split('__')[1].lower()

				## include ids
				tax_name[0] = id_name[0] + ';'+ tax_name[0]

				## label names with complete 7 levels with '>' in the beginning
				if len(tax_name) >= 7: tax_name[0] = '>'+tax_name[0]
				# elif len(tax_name) > 7:
				# 	tax_name[0] = '!'+tax_name[0]
				# 	print(';'.join(tax_name))
				tax_name = ';'.join(tax_name)

				## replace id with name
				if tax_name in gg_dict.keys():
					gg_dict[tax_name].append(gg_dict[id_name[0]][0])
					del gg_dict[id_name[0]]
				else:
					gg_dict[tax_name] = gg_dict.pop(id_name[0])

				if cnt % 10000 == 0 : print(tax_name, len(gg_dict[tax_name]))
			else:
				not_found.append(id_name[0])

	cnt = 0
	for i in gg_dict.keys():
		if i.startswith('>'): cnt += 1
	print(f'{cnt} out of {len(gg_dict.keys())} names are complete in 7 levels.')
	print(f'{len(not_found)} ids are not found in fasta file.')
	return gg_dict

def compareSeqsCheck(longSeq, shortSeq):
	assert len(longSeq) >= len(shortSeq)
	cnt = -1
	index = 0
	while index < len(longSeq) and index > -1:
		index = longSeq.find(shortSeq, index)
		if index != -1: index += len(shortSeq)
	cnt += 1
	if index > -1 : print("found!")
	return index, cnt

def compareSeqs(longSeq, shortSeq):
	assert len(longSeq) >= len(shortSeq)
	index = 0
	index = longSeq.find(shortSeq, index)
	if index > -1 : print("found!")
	return index

def saveFile(IDB):
	IDBFasta = open('IDB_v5.fasta', 'w')
	IDBTxt = open('IDB_v5.txt', 'w')
	for i in IDB.keys():
		# if not i.startswith('>'): ids = '>' + i
		ids = i if i.startswith('>') else '>' + i
		ids = ids.split(';')[0]
		for j in IDB[i]:
			IDBFasta.write(ids+'\n')
			IDBFasta.write(j+'\n')
		IDBTxt.write(ids[1:]+'\t'+';'.join(i.split(';')[1:])+'\n')
	IDBFasta.close()
	IDBTxt.close()

def check(gg_dict):
	gg_array = []
	for i in gg_dict.keys():
		for j in range(len(gg_dict[i])):
			gg_array.append(gg_dict[i][j])
	for k in range(len(gg_array)):
		for l in range(k+1, len(gg_array)):
			if len(gg_array[k]) > len(gg_array[l]):
				index = gg_array[k].find(gg_array[l])
			else:
				index = gg_array[l].find(gg_array[k])
			if index > -1 :
				# print('oh no: ', gg_array[k], gg_array[l], index)

				for a in gg_dict.keys():
					if k > len(gg_dict[a]):
						k -= len(gg_dict[a])
						l -= len(gg_dict[a])
					else:
						print(a, k, l)
						break
		if k % 10000 == 0: print(f'passing {k} ...')
		# 172726 161964
'''
def integrateDB(gg_dict, slv_dict, IDB, reverse=0):
	for i in gg_dict.keys():
		for gg_seq in gg_dict[i]:
			found = 0
			print(len(gg_seq))
			for j in slv_dict.keys():
				for slv_seq in slv_dict[j]:
					if gg_seq == slv_seq and not reverse:
						found = 1
						if i.startswith('>'):
							if i in IDB.keys():
								IDB[i].append(slv_seq)
							else:
								IDB[i] = [slv_seq]
						else:
							if j in IDB.keys():
								IDB[j].append(slv_seq)
							else:
								IDB[j] = [slv_seq]
					elif len(gg_seq) > len(slv_seq): continue
					index, double = compareSeqs(slv_seq, gg_seq)
					if index > -1 and not double:
						found = 1
						if not reverse:
							if i.startswith('>'):
								if i in IDB.keys():
									IDB[i].append(slv_seq)
								else:
									IDB[i] = [slv_seq]
							else:
								if j in IDB.keys():
									IDB[j].append(slv_seq)
								else:
									IDB[j] = [slv_seq]
						else:
							if j.startswith('>'):
								if j in IDB.keys():
									IDB[j].append(slv_seq)
								else:
									IDB[j] = [slv_seq]
							else:
								if i in IDB.keys():
									IDB[i].append(slv_seq)
								else:
									IDB[i] = [slv_seq]
			if not found:
				if i in IDB.keys():
					IDB[i].append(gg_seq)
				else:
					IDB[i] = [gg_seq]
				if len(IDB) >= 3 : return IDB
			else: return IDB

def integrateDB(gg_dict, slv_dict):
	IDB = slv_dict
	for i in gg_dict.keys():
		for seq1 in gg_dict[i]:
			found = 0
			for j in slv_dict.keys():
				for seq2 in slv_dict[j]:
					if seq1 == seq2:
						continue
					elif len(seq1) > len(seq2):
						index, double = compareSeqs(seq1, seq2)
					else:
						index, double = compareSeqs(seq2, seq1)

					if index > -1 and not double:
						print(j, i)
						if j.startswith('>'):
							print("using gg ids")
							if j in IDB.keys():
								IDB[j].append(seq1 if len(seq1)>len(seq2) else seq2)
							else:
								IDB[j] = [seq1] if len(seq1)>len(seq2) else [seq2]
						else:
							print("using silva ids")
							if i in IDB.keys():
								IDB[i].append(seq1 if len(seq1)>len(seq2) else seq2)
							else:
								IDB[i] = [seq1] if len(seq1)>len(seq2) else [seq2]
					elif index > -1 and double:
						print("found twice in a sequence!")
					else:
						if i in IDB.keys():
							IDB[i].append(seq1)
						else:
							IDB[i] = [seq1]
						# if j in IDB.keys():
						# 	IDB[j].append(seq2)
						# else:
						# 	IDB[j] = [seq2]

			if not found:
				if i in IDB.keys():
					IDB[key]=[seq1]
				else:
					IDB[i] = [seq1]
'''
def rename(gg_tax, slv_tax):
	gg_tax = gg_tax.split(';')
	slv_tax = slv_tax.split(';')
	compare(gg_tax, slv_tax)

def integrateDB(gg_dict, slv_dict):
	IDB = slv_dict.copy()
	undone = len(gg_dict)
	tStart = time.time()
	for gg_key in gg_dict.keys():
		undone -= 1

		for gg_seq in gg_dict[gg_key]:
			found = 0
			for slv_key in slv_dict.keys():
				for idx, slv_seq in enumerate(slv_dict[slv_key]):
					if len(gg_seq) > len(slv_seq):
						index = gg_seq.find(slv_seq)
						if index > -1:
							print(f'Costing time: {(time.time()-tStart)/60}, {undone}/{len(gg_dict)}')
							print("found: ", gg_key, slv_key)
							found = 1
							IDB[slv_key][idx] = gg_seq
							# if found and gg sequence has replaced silva sequence, break
							break
					else:
						index = slv_seq.find(gg_seq)
						if index > -1:
							found = 1

				if found: break

			if not found:
				if gg_key in IDB.keys():
					IDB[gg_key].append(gg_seq)
				else:
					IDB[gg_key] = [gg_seq]

	return IDB
'''
def integrateDB(gg_dict, slv_dict):
	IDB = slv_dict.copy()
	undone = len(gg_dict)
	tStart = time.time()
	for gg_key in gg_dict.keys():
		undone -= 1
		tokens = gg_key.split(';')
		for gg_seq in gg_dict[gg_key]:
			found = 0
			for slv_key in slv_dict.keys():
				cnt = 0
				for token in tokens[1:]:
					if slv_key.find(token) != -1:
						cnt += 1
					else: break

				if cnt <= 3*len(tokens-1)/4: continue
				# else:
				# 	if print_cnt % 100000 == 0: print(f'similarity: {cnt}/{len(tokens)}\ngg: {gg_key}, slv: {slv_key}')			
		
				for idx, slv_seq in enumerate(slv_dict[slv_key]):
					if len(gg_seq) > len(slv_seq):
						index = gg_seq.find(slv_seq)
						if index > -1:
							print(f'Costing time: {(time.time()-tStart)/60}, {undone}/{len(gg_dict)}')
							print("found: ", gg_key, slv_key)
							found = 1
							IDB[slv_key][idx] = gg_seq
					else: # len(gg_seq) < len(slv_seq):
						index = slv_seq.find(gg_seq)
						if index > -1:
							found = 1

			if not found:
				if gg_key in IDB.keys():
					IDB[gg_key].append(gg_seq)
				else:
					IDB[gg_key] = [gg_seq]

	return IDB

def integrateDB(gg_dict, slv_dict):
	IDB = slv_dict.copy()
	for gg_key in gg_dict.keys():
		for gg_seq in gg_dict[gg_key]:
			found = 0
			for slv_key in slv_dict.keys():
				for slv_seq in slv_dict[slv_key]:
					if gg_seq == slv_seq:
						continue
					elif len(gg_seq) > len(slv_seq):
						index = compareSeqs(gg_seq, slv_seq)
					else:
						index = compareSeqs(slv_seq, gg_seq)
					
					if index > -1 #and not double:
						found = 1
						if len(gg_seq) > len(slv_seq):

						if gg_key.startswith('>'):
							print("using gg ids")
							if gg_key in IDB.keys():
								IDB[gg_key].append(gg_seq if len(gg_seq)>len(slv_seq) else slv_seq)
							else:
								IDB[gg_key] = [gg_seq] if len(gg_seq)>len(slv_seq) else [slv_seq]
						else:
							print("using silva ids")
							if i in IDB.keys():
								IDB[slv_key].append(gg_seq if len(gg_seq)>len(slv_seq) else slv_seq)
							else:
								IDB[slv_key] = [gg_seq] if len(gg_seq)>len(slv_seq) else [slv_seq]
					# elif index > -1 and double:
					# 	print("found twice in a sequence!")
					else:
						continue

			if not found:
				if gg_key in IDB.keys():
					IDB[gg_key].append(gg_seq)
				else:
					IDB[gg_key] = [gg_seq]
			# else:
			# 	return IDB

	return IDB
'''
if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('--ggFasta',	type=str  , default='99_otus.fasta')
	parser.add_argument('--ggTxt',		type=str  , default='99_otu_taxonomy.txt')
	parser.add_argument('--slvFasta',	type=str  , default='silva_132_99_16S.fna')
	parser.add_argument('--slvTxt',		type=str  , default='taxonomy_all_levels.txt')
	args = parser.parse_args()

	gg_dict = Preprocess(fastaFile=args.ggFasta, txtFile=args.ggTxt)
	slv_dict = Preprocess(fastaFile=args.slvFasta, txtFile=args.slvTxt)

	gg_cnt = 0
	slv_cnt = 0
	for i in slv_dict.keys():
		slv_cnt+=len(slv_dict[i])
	for i in gg_dict.keys():
		gg_cnt+=len(gg_dict[i])

	print(f'There are {gg_cnt} in GG, {slv_cnt} in Slv, computing {gg_cnt*slv_cnt}')
	# IDB = {}
	# check(gg_dict)
	IDB = integrateDB(gg_dict, slv_dict)
	# IDB = integrateDB(slv_dict, gg_dict, IDB, 1)
	saveFile(IDB)
