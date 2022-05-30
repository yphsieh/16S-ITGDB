import os
import re

from joblib import load

def ValidateSINTAXFullTaxa(s_correct_file, s_sintax_file):
	
	dict_taxaMap = load("/work1782/cil6758/Project/cpclf/val/taxa_map.joblib")

	fpc = open(s_correct_file, 'r')
	fps = open(s_sintax_file, 'r')
	
	dict_sintax_taxa = {}	#dict_sintax_taxa = {ID:l_taxa}
	
	for line in fps:
		line = re.sub(r"\(\d+\.\d+\)", "", line)
		line = re.sub(r"[dpcofgs]:", "", line)
		l_ele = line.strip().split('\t')
		
		s_ID = l_ele[0]
		l_taxa = l_ele[1].split(',')
		dict_sintax_taxa[s_ID] = l_taxa
	
	dict_answer_taxa = {}	#dict_correct_taxa = {ID:l_taxa}
	
	fps.close()
	
	if(s_correct_file.find("dairy") != -1):
		for line in fpc:
			l_ele = line.strip().split('\t')
			s_ID = l_ele[0]
			l_taxa = l_ele[8].strip().split(';')
			dict_answer_taxa[s_ID] = l_taxa
			
	elif(s_correct_file.find("udb") != -1):
		for line in fpc:
			l_ele = line.strip().split('\t')
			s_ID = l_ele[0]
			l_taxa = l_ele[1].strip().split('; ')
			l_taxa = [re.sub(r"[kdpcofgs]__", "", iter) for iter in l_taxa]
			dict_answer_taxa[s_ID] = l_taxa
	elif(s_correct_file.find("exclude") != -1 or s_correct_file.find("le2overlap") != -1 or s_correct_file.find("union") != -1):
		for line in fpc:
			l_ele = line.strip().split('\t')
			s_ID = l_ele[0]
			l_taxa = l_ele[1].strip().split('; ')
			l_taxa = [re.sub(r"[kdpcofgs]__", "", iter) for iter in l_taxa]
			dict_answer_taxa[s_ID] = l_taxa
	elif(s_correct_file.find("mock") != -1):
		for line in fpc:
			l_ele = line.strip().split('\t')
			s_ID = l_ele[0]
			l_taxa = l_ele[1].strip().split(';')
			dict_answer_taxa[s_ID] = l_taxa
	else:
		pass
	
	fpc.close()
	
	l_correct = [0,0,0,0,0,0,0]
	
	l_predict_class_count = [{},{},{},{},{},{},{}]
	l_answer_class_count = [{},{},{},{},{},{},{}]
	l_correct_assign_count = [{},{},{},{},{},{},{}]
	
	for skey,staxa in dict_sintax_taxa.items():
		for i in range(len(staxa)):
			if(staxa[i] not in l_predict_class_count[i].keys()):
				l_predict_class_count[i][staxa[i]] = 0
			else:
				l_predict_class_count[i][staxa[i]] += 1
			
	for akey,ataxa in l_answer_class_count.items():
		for i in range(len(ataxa)):
			if(ataxa[i] not in l_answer_class_count[i].keys()):
				l_answer_class_count[i][ataxa[i]] = 0
			else:
				l_answer_class_count[i][ataxa[i]] += 1
	
	for rkey,rtaxa in dict_sintax_taxa.items():
	
		b_isFullTaxaMatch = False
	
		if(rtaxa[5].lower() == dict_answer_taxa[rkey][5].lower()):
			l_sintax_sp = rtaxa[6].split('_')
			l_correct_sp = dict_answer_taxa[rkey][6].split('_')
			
			if(len(l_sintax_sp) > 1):
				if(len(l_correct_sp) > 1):
					if(l_sintax_sp[1].lower() == l_correct_sp[1].lower()):
						b_isFullTaxaMatch = True
				else:
					if(l_sintax_sp[1].lower() == l_correct_sp[0].lower()):
						b_isFullTaxaMatch = True
			else:
				if(len(l_correct_sp) > 1):
					if(l_sintax_sp[0].lower() == l_correct_sp[1].lower()):
						b_isFullTaxaMatch = True
				else:
					if(l_sintax_sp[0].lower() == l_correct_sp[0].lower()):
						b_isFullTaxaMatch = True
		
		if(b_isFullTaxaMatch == True):
			for i in range(7):
				l_correct[i] += 1
				
				if(rtaxa not in l_correct_assign_count[i].keys())
					l_correct_assign_count[i][rtaxa[i]] = 1
				else:
					l_correct_assign_count[i][rtaxa[i]] += 1
				
		else:
			for i in range(6):
				if(rtaxa[i].lower().find(dict_answer_taxa[rkey][i].lower()) != -1 or dict_answer_taxa[rkey][i].lower().find(rtaxa[i].lower()) != -1):
				#if(rtaxa[i].lower() == dict_answer_taxa[rkey][i].lower()):
					l_correct[i] += 1
					
					if(rtaxa not in l_correct_assign_count[i].keys())
						l_correct_assign_count[i][rtaxa[i]] = 1
					else:
						l_correct_assign_count[i][rtaxa[i]] += 1
					
				else:
					if(rtaxa[i] in dict_taxaMap.keys()):
						if(dict_answer_taxa[rkey][i] in dict_taxaMap[rtaxa[i]]):
							l_correct[i] += 1
							
							if(rtaxa not in l_correct_assign_count[i].keys())
								l_correct_assign_count[i][rtaxa[i]] = 1
							else:
								l_correct_assign_count[i][rtaxa[i]] += 1
							
						else:
							# if(i == 4):
								# print(rtaxa[i], dict_answer_taxa[rkey][i])
							break
					else:
						# if(i == 4):
							# print(rtaxa[i], dict_answer_taxa[rkey][i])
						break
	
	l_acc = [round(iter/len(dict_sintax_taxa), 4)  for iter in l_correct]
	
	l_precision = [{},{},{},{},{},{},{}]
	l_recall = [{},{},{},{},{},{},{}]
	
	for i in range(7):
		for ctaxa,cval in l_correct_assign_count[i].items():
			l_precision[i][ctaxa] = cval/l_predict_class_count[i][ctaxa]
			l_recall[i][ctaxa] = cval/l_answer_class_count[i][ctaxa]
			
	l_macro_precision = [0,0,0,0,0,0,0]
	l_macro_recall = [0,0,0,0,0,0,0]
	
	for i in range(7):
		f_sum = 0
		for ctaxa, cval in l_precision[i].items():
			f_sum += cval
			
		l_macro_precision[i] = f_sum / len(l_precision[i])
		l_macro_recall[i] = f_sum / len(l_recall[i])
	
	print("Accuracy: ", l_acc)
	print("Macro precision: ", l_macro_precision)
	print("Macro recall: ", l_macro_recall)
		
def ValidateQIIMEFullTaxa(s_correct_file, s_qiime_file):
	
	dict_taxaMap = load("/work1782/cil6758/Project/cpclf/val/taxa_map.joblib")

	fpc = open(s_correct_file, 'r')
	fps = open(s_qiime_file, 'r')
	
	dict_qiime_taxa = {}	#dict_sintax_taxa = {ID:l_taxa}
	
	fps.readline()
	
	for line in fps:
		line = re.sub(r"[kpcofgs]__", "", line)
		l_ele = line.strip().split('\t')
		
		s_ID = l_ele[0]
		l_taxa = l_ele[1].split('; ')
		l_taxa = [iter.strip() for iter in l_taxa]
		
		dict_qiime_taxa[s_ID] = l_taxa
	
	fps.close()
	
	dict_correct_taxa = {}	#dict_correct_taxa = {ID:l_taxa}
	
	if(s_correct_file.find("dairy") != -1):
		for line in fpc:
			l_ele = line.strip().split('\t')
			s_ID = l_ele[0]
			l_taxa = l_ele[8].strip().split(';')
			dict_correct_taxa[s_ID] = l_taxa
			
	elif(s_correct_file.find("udb") != -1):
		for line in fpc:
			l_ele = line.strip().split('\t')
			s_ID = l_ele[0]
			l_taxa = l_ele[1].strip().split('; ')
			l_taxa = [re.sub(r"[kdpcofgs]__", "", iter) for iter in l_taxa]
			dict_correct_taxa[s_ID] = l_taxa
	elif(s_correct_file.find("exclude") != -1 or s_correct_file.find("le2overlap") != -1 or s_correct_file.find("union") != -1):
		for line in fpc:
			l_ele = line.strip().split('\t')
			s_ID = l_ele[0]
			l_taxa = l_ele[1].strip().split('; ')
			l_taxa = [re.sub(r"[kdpcofgs]__", "", iter) for iter in l_taxa]
			dict_correct_taxa[s_ID] = l_taxa
	else:
		pass
	
	fpc.close()
	
	l_acc = [0,0,0,0,0,0,0]
	
	for rkey,rtaxa in dict_qiime_taxa.items():
	
		b_isFullTaxaMatch = False
	
		if(rtaxa[5].lower() == dict_correct_taxa[rkey][5].lower()):
			l_qiime_sp = rtaxa[6].split('_')
			l_correct_sp = dict_correct_taxa[rkey][6].split('_')
			
			if(len(l_qiime_sp) > 1):
				if(len(l_correct_sp) > 1):
					if(l_qiime_sp[1].lower() == l_correct_sp[1].lower()):
						b_isFullTaxaMatch = True
				else:
					if(l_qiime_sp[1].lower() == l_correct_sp[0].lower()):
						b_isFullTaxaMatch = True
			else:
				if(len(l_correct_sp) > 1):
					if(l_qiime_sp[0].lower() == l_correct_sp[1].lower()):
						b_isFullTaxaMatch = True
				else:
					if(l_qiime_sp[0].lower() == l_correct_sp[0].lower()):
						b_isFullTaxaMatch = True
		
		if(b_isFullTaxaMatch == True):
			for i in range(7):
				l_acc[i] += 1
		else:
			for i in range(6):
				if(rtaxa[i].lower().find(dict_correct_taxa[rkey][i].lower()) != -1 or dict_correct_taxa[rkey][i].lower().find(rtaxa[i].lower()) != -1):
				#if(rtaxa[i].lower() == dict_correct_taxa[rkey][i].lower()):
					l_acc[i] += 1
				else:
					if(rtaxa[i] in dict_taxaMap.keys()):
						if(dict_correct_taxa[rkey][i] in dict_taxaMap[rtaxa[i]]):
							l_acc[i] += 1
						else:
							# if(i == 4):
								# print(rtaxa[i], dict_correct_taxa[rkey][i])
							break
					else:
						# if(i == 4):
							# print(rtaxa[i], dict_correct_taxa[rkey][i])
						break
	
	l_acc = [round(iter/len(dict_qiime_taxa), 4)  for iter in l_acc]
	
	print(l_acc)
	


	
	
