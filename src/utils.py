import os
import re
import csv
import pandas as pd

class Microbiome():
	def __init__(self):
		self.s_taxa_fullname = ""
		self.l_taxa = []
		self.s_seq = ""
		self.b_is_rdp_overlap = False
		self.b_is_slv_overlap = False
		self.b_is_gg_overlap = False
		self.b_is_union_collect = True

def Check_Microbiome_content(dict_mic):
	for key, val in dict_mic.items():
		print(key)
		print(val.s_taxa_fullname)
		print(val.l_taxa)
		print(val.s_seq)
		print("======================================")

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

def ParseFile(s_seq_file, s_taxa_file, s_db):
	dict_mic = {}
	
	if(s_db == "slv"):	#Parse SILVA database
		with open(s_seq_file, "r") as fp:
			
			for line in fp:
				sd_temp_mic = Microbiome()
				
				if(line[0] == '>'):
					s_ID = line[1:-1].strip()
					line = fp.readline()
					sd_temp_mic.s_seq = line.upper().strip()
					
					if(s_ID not in dict_mic):
						dict_mic[s_ID] = sd_temp_mic
					else:
						print("Repeated ID: " + s_ID)
		
		with open(s_taxa_file, "r") as fp:
			for line in fp:
				l_ele = line.strip().split("\t")
				s_ID = l_ele[0].strip()
				s_taxa_name = l_ele[1].strip()
				l_taxa = s_taxa_name.split(';')
				l_taxa = [re.sub(r"^[kpcofgsd]__", "", iter.strip()) for iter in l_taxa]
				
				dict_mic[s_ID].l_taxa = l_taxa
					
		return dict_mic
						
	elif(s_db == "rdp"):	#Parse RDP database
		
		with open(s_seq_file, "r") as fp:
			
			for line in fp:
				sd_temp_mic = Microbiome()
				
				if(line[0] == '>'):
					s_ID = line[1:-1].strip()
					line = fp.readline()
					sd_temp_mic.s_seq = line.upper().strip()
					
					if(s_ID not in dict_mic):
						dict_mic[s_ID] = sd_temp_mic
					else:
						print("Repeated ID: " + s_ID)
		
		with open(s_taxa_file, "r") as fp:
			for line in fp:
				l_ele = line.strip().split("\t")
				s_ID = l_ele[0].strip()
				s_taxa_name = l_ele[1].strip()
				l_taxa = s_taxa_name.split(';')
				l_taxa = [re.sub(r"^[kpcofgsd]__", "", iter.strip()) for iter in l_taxa]
				
				dict_mic[s_ID].l_taxa = l_taxa
					
		return dict_mic
	
	elif(s_db == "gg"):	#Parse Greengenes database
		with open(s_seq_file, "r") as fp:
			
			for line in fp:
				sd_temp_mic = Microbiome()
				
				if(line[0] == '>'):
					s_ID = line[1:-1].strip()
					line = fp.readline()
					sd_temp_mic.s_seq = line.upper().strip()
					
					if(s_ID not in dict_mic):
						dict_mic[s_ID] = sd_temp_mic
					else:
						print("Repeated ID: " + s_ID)
		
		set_sptaxa = set()
		
		with open(s_taxa_file, "r") as fp:
			for line in fp:
				l_ele = line.strip().split("\t")
				s_ID = l_ele[0].strip()
				s_taxa_name = l_ele[1].strip()
				l_taxa = s_taxa_name.split(';')
				l_taxa = [re.sub(r"^[kpcofgsd]__", "", iter.strip()) for iter in l_taxa]
				
				if(l_taxa[6] != ''):
					set_sptaxa.add(s_taxa_name)
				
				dict_mic[s_ID].l_taxa = l_taxa
				
		wp = open("sptaxa_gg.txt", 'w')

		for iter in set_sptaxa:
			wp.write(iter + "\n")
					
		return dict_mic
		
	elif(s_db == "itg"):	#Parse integrated database
		with open(s_seq_file, "r") as fp:
			
			for line in fp:
				sd_temp_mic = Microbiome()
				
				if(line[0] == '>'):
					s_ID = line[1:-1].strip()
					line = fp.readline()
					sd_temp_mic.s_seq = line.upper().strip()
					
					if(s_ID not in dict_mic):
						dict_mic[s_ID] = sd_temp_mic
					else:
						print("Repeated ID: " + s_ID)
		
		
		with open(s_taxa_file, "r") as fp:
			for line in fp:
				l_ele = line.strip().split("\t")
				s_ID = l_ele[0].strip()
				s_taxa_name = l_ele[1].strip()
				l_taxa = s_taxa_name.split(';')
				l_taxa = [re.sub(r"^[kpcofgsd]__", "", iter.strip()) for iter in l_taxa]
				
				dict_mic[s_ID].l_taxa = l_taxa
				
		return dict_mic

	else:
		print("Invalid database. The acceptable database are: SILVA, Greengenes, and RDP")
		sys.exit()
	