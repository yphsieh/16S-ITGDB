import os
import sys
import csv
import joblib
import argparse
import numpy as np
import pandas as pd
from copy import deepcopy as cp
import matplotlib.pyplot as plt

base_dir = '/work1782/yphsieh/7IDB/'

taxa_map = joblib.load(os.path.join(base_dir, os.path.expanduser('taxa_map.joblib')))
tax_keys = list(taxa_map.keys())

for cnt, tax_key in enumerate(tax_keys):
	if tax_key == 'all' or tax_key == 'root': continue
	new_key = tax_key.lower().strip('[').strip(']').strip('\'') 
	taxa_map[new_key] = taxa_map[tax_key]
	for idx in range(len(taxa_map[new_key])):
		taxa_map[new_key][idx] = taxa_map[new_key][idx].lower().strip('[').strip(']').strip('\'')
	
	genus_key = new_key.split(' ')[0]
	if genus_key.lower() == 'unknown': continue
	if genus_key not in taxa_map.keys(): taxa_map[genus_key] = cp(taxa_map[new_key])
	else: taxa_map[genus_key].extend(taxa_map[new_key])

	for idx in range(1, len(taxa_map[genus_key])+1):
		taxa_map[genus_key][-idx] = taxa_map[genus_key][-idx].split(' ')[0] 
	
	taxa_map[genus_key] = list(set(taxa_map[genus_key]))
	
def removeNoInfo(taxlevel):
	# remove unclassified
	a = 1
	level = len(taxlevel)
	while (taxlevel[-a] == '?' or taxlevel[-a].startswith('unclassified') or taxlevel[-a].startswith('unidentified') or taxlevel[-a].startswith('unassigned') or taxlevel[-a].startswith('uncultured') or taxlevel[-a].find('incertae') != -1 or taxlevel[-a].find('sedis') != -1) and (a<level):
		a += 1
		level -= 1
	taxlevel = taxlevel[:level]
	
	# remove genus in species
	if level >= 7:
		if len(taxlevel[6].split(' ')) > 1 and taxlevel[6].startswith(taxlevel[5].split(' ')[0]):
			if taxlevel[6].find('symbiont of') != -1:
				taxlevel[6] = 'endosymbiont'
			else:
				taxlevel[6] = ' '.join(taxlevel[6].split(' ')[1:])
	
	taxlevel = list(filter(None, taxlevel))
	return ';'.join(taxlevel)	

def cleanse_sintax(taxon):
	if pd.isna(taxon): taxon = ''
	taxlevel = taxon.split(',')

	# remove no information levels
	for i in range(len(taxlevel)):
		if taxlevel[i].strip(' ') == '':
			taxlevel[i] = '?'
			#taxlevel = taxlevel[:i]
			#break

	# remove 'k:' in the front
	for l in range(len(taxlevel)):
		taxlevel[l] = taxlevel[l].split(':')[-1].lower().strip()
		taxlevel[l] = taxlevel[l].replace('_',' ').strip()
	
	# remove confidence in the back
	for l in range(len(taxlevel)):
		taxlevel[l] = taxlevel[l].split('(')[0].lower().strip()

	return removeNoInfo(taxlevel)

def cleanse_spingo(taxon, gecies):
	if taxon.lower() == 'ambiguous': taxon = 'ambiguous;ambiguous;ambiguous;ambiguous;ambiguous;ambiguous'
	taxlevel = taxon.split(';')

	# combine genus species with taxon
	if taxlevel[0].lower() == 'ambiguous' and gecies.lower() != 'ambiguous':
		taxlevel[-1] = gecies.split('_')[0]
		taxlevel.append(gecies.split('_')[1])
	else:
		taxlevel.append(gecies)
	
	if len(taxlevel) != 7: print(taxlevel)

	# remove no information levels
	for l in range(len(taxlevel)):
		taxlevel[l] = taxlevel[l].replace('_',' ').strip(' ')

	return removeNoInfo(taxlevel)


def cleanse_mothur(taxon):
	taxlevel = taxon.split(';')
	taxlevel[0] = taxlevel[0].split('::')[-1]

	# remove no information levels
	for i in range(len(taxlevel)):
		if taxlevel[i].endswith('__') or taxlevel[i].strip(' ') == '':
			taxlevel[i] = '?'
			#taxlevel = taxlevel[:i]
			#break
	
	# remove 'k__' in the front
	for l in range(len(taxlevel)):
		taxlevel[l] = taxlevel[l].split('__')[-1].lower()
		taxlevel[l] = taxlevel[l].replace('_',' ').strip()

	# remove confidence in the back
	for l in range(len(taxlevel)):
		taxlevel[l] = taxlevel[l].split('(')[0].lower().strip()


	return removeNoInfo(taxlevel)


def cleanse(taxon):
	taxlevel = taxon.split(';')
	taxlevel[0] = taxlevel[0].split('::')[-1]

	# remove no information levels
	for i in range(len(taxlevel)):
		if taxlevel[i].endswith('__') or taxlevel[i].strip(' ') == '':
			taxlevel[i] = '?'
	
	# remove 'k__' in the front
	for l in range(len(taxlevel)):
		taxlevel[l] = taxlevel[l].split('__')[-1].lower()
		taxlevel[l] = taxlevel[l].replace('_',' ').strip()

	return removeNoInfo(taxlevel)

def Preprocess(tsvFile='new_idb_taxonomy_12_ext.tsv'):
	taxa = pd.read_csv(tsvFile, sep='\t', header = 0)
	taxa = taxa.to_dict('records')
	taxa_dict = {}
	for i in taxa:
		#if i['Taxon'] in taxa_dict.keys():
		if i['Feature ID'] in taxa_dict.keys():
			print(f"{i['Feature ID']} exists !! ({len(taxa_dict[i['Feature ID']])})")
			taxa_dict[i['Feature ID']].append(cleanse(i['Taxon']))
		else:
			taxa_dict[i['Feature ID']] = cleanse(i['Taxon'])

	return taxa_dict

def taxaLevel(y):
	y = y.lower()
	level = len(y.split(';'))
	token = y.split(';')[:level]
	return token, level

def checkMap(taxon1, taxon2):
	if taxon1 == 'escherichia' or taxon1 == 'shigella':
		if taxon2.find(taxon1) != -1:
			return True
	if taxon1 in taxa_map.keys():
		if taxon2 in taxa_map[taxon1]:
			return True
	return False

def compare(taxon1, taxon2):
	taxon1, taxon2 = taxon1.lower().strip('[').strip(']').strip('\''), taxon2.lower().strip('[').strip(']').strip('\'')
	if taxon1 == taxon2:
		return True
	
	if checkMap(taxon1, taxon2) or checkMap(taxon2, taxon1):
		return True
	
	return False


def algos(name):
	f = open(name + '_algos.csv', 'w')
	print(f'Test case,Algorithm,Level,Unassigned,Accuracy,Precision,Recall,F1-Score,Weighted F1-Score', file = f)
	f.close()
	
	to_plot_am, to_plot_ac = [], []
	algos = ['sintax', 'spingo', 'mothur', 'qiime']
	classifier = 'itgdb_v1'

	
	for algo in range(len(algos)):
		print(f'--- Analyzing the result of {algos[algo]} on {name} ---')
		tsvFile = os.path.join(base_dir, classifier, classifier + '_taxonomy_' + name + '.' + algos[algo])	
		
		if algos[algo] == 'sintax':
			taxa = pd.read_csv(tsvFile, sep='\t', names = ['Feature ID', 'Conf', '+', 'Taxon'])
			taxa = taxa.to_dict('records')
			y = {}
			for i in taxa:
				if i['Feature ID'] in y.keys():
					print(f"{i['Feature ID']} exists !!")
					y[i['Feature ID']].append(cleanse_sintax(i['Conf']))
				else:
					y[i['Feature ID']] = cleanse_sintax(i['Conf'])

		elif algos[algo] == 'mothur':
			taxa = pd.read_csv(tsvFile, sep='\t', names = ['Feature ID', 'Conf'])
			taxa = taxa.to_dict('records')
			y = {}
			for i in taxa:
				if i['Feature ID'] in y.keys():
					print(f"{i['Feature ID']} exists !!")
					y[i['Feature ID']].append(cleanse_mothur(i['Conf']))
				else:
					y[i['Feature ID']] = cleanse_mothur(i['Conf'])

		elif algos[algo] == 'qiime':
			y = Preprocess(tsvFile)
		
		elif algos[algo] == 'spingo':
			taxa = pd.read_csv(tsvFile, sep='\t', names = ['Feature ID', '1', '2', '3', 'Taxon', '4', 'Gecies', '5'])
			taxa = taxa.to_dict('records')
			y = {}
			for i in taxa:
				if i['Feature ID'] in y.keys():
					print(f"{i['Feature ID']} exists !!")
					y[i['Feature ID']].append(cleanse_spingo(i['Taxon'], i['Gecies']))
				else:
					y[i['Feature ID']] = cleanse_spingo(i['Taxon'], i['Gecies'])


		if os.path.isfile(os.path.join(base_dir, name, name + '_taxonomy.txt')):
			print('with file')

			tgt_file = os.path.join(base_dir, name, name + '_taxonomy.txt')
			targets = Preprocess(tgt_file)

		else:
			print('with ids')

			tgt_tofile = open(os.path.join(base_dir, name + '_taxonomy.txt'), 'w')
			print(f'FeatureID\tTaxon', file=tgt_tofile)

			targets = {}
			for i in y.keys():
				print(f"{i.split('::')[0]}\t{i.split('::')[1]}", file=tgt_tofile)
				clean = cleanse(i)
				targets[i] = clean
			tgt_tofile.close()

		eachLevel, correct_level, unassigned, totalLevel, totalLevel, acc_frac_sum, over, correct_taxa, P, R, macroF1, wF1 = analysis_with_tgt(y, targets, name)


		f = open( name + '_algos.csv', 'a')
		labels = ['Family', 'Genus', 'Species']
		for lev in range(4,7):
			print(f'{output_name(name)},{output_name(algos[algo])},{labels[lev-4]},{unassigned/len(y.keys())*100:4.2f}%,{correct_level[lev]/len(y.keys())*100:4.2f}%,{P[lev]*100:4.2f},{R[lev]*100:4.2f},{macroF1[lev]*100:4.2f},{wF1[lev]*100:4.2f}', file = f)
		f.close()

		to_plot_ac.append(correct_level[-3:]/len(y.keys()))
		to_plot_am.append(eachLevel[-3:])
	
	plot(name, to_plot_ac, algos)
	plotPR(name, 'algos')


def sintax(name):
	f = open(name + '.csv', 'w')
	print(f'Test case,DB,Level,Unassigned,Accuracy,Precision,Recall,F1-Score,Weighted F1-Score', file = f)
	f.close()
	
	to_plot_am, to_plot_ac = [], []
	classifiers = ['rdp', 'slv', 'gg', 'grd', 'udb', 'gtdb', 'rsgidb', 'itgdb_v1']
	# classifiers = ['itgdb_v1']

	
	for classifier in range(len(classifiers)):
		print(f'--- Analyzing the result of {classifiers[classifier]} on {name} ---')
		tsvFile = os.path.join(base_dir, classifiers[classifier], classifiers[classifier] + '_taxonomy_' + name +'.sintax')	

		taxa = pd.read_csv(tsvFile, sep='\t', names = ['Feature ID', 'Conf', '+', 'Taxon'])
		taxa = taxa.to_dict('records')
		y = {}

		for i in taxa:
			if i['Feature ID'] in y.keys():
				print(f"{i['Feature ID']} exists !!")
				y[i['Feature ID']].append(cleanse_sintax(i['Conf'])) #
			else:
				y[i['Feature ID']] = cleanse_sintax(i['Conf']) #

		
		tgt_file = os.path.join(base_dir, name, name + '_taxonomy.txt')
		targets = Preprocess(tgt_file)

		eachLevel, correct_level, unassigned, totalLevel, totalLevel, acc_frac_sum, over, correct_taxa, P, R, macroF1, wF1 = analysis_with_tgt(y, targets, name)

		f = open( name + '.csv', 'a')
		labels = ['Family', 'Genus', 'Species']
		for lev in range(4,7):
			print(f'{output_name(name)},{output_name(classifiers[classifier])},{labels[lev-4]},{unassigned/len(y.keys())*100:4.2f}%,{correct_level[lev]/len(y.keys())*100:4.2f}%,{P[lev]*100:4.2f},{R[lev]*100:4.2f},{macroF1[lev]*100:4.2f},{wF1[lev]*100:4.2f}', file = f)
		f.close()

		to_plot_ac.append(correct_level[-3:]/len(y.keys()))
		to_plot_am.append(eachLevel[-3:])
	
	plot(name, to_plot_ac, classifiers)
	plotPR(name, 'db')

def analysis_with_tgt(y, targets, name):

	unassigned, totalLevel, acc_frac_sum = 0, 0, 0
	eachLevel = np.zeros(7)

	correct_level = np.zeros(7)
	over = 0
		
	correct_taxa=[]
	PR_dict_list = [{}, {}, {}, {}, {}, {}, {}]
	
	for i in y.keys():
				
		target = targets[i].split(';')
		token, level = taxaLevel(y[i])

		# count the number to each level
		for j in range(level): eachLevel[j] += 1
		# count the total level
		totalLevel += level
		# count the unassigned
		if level <= 1:
			unassigned += 1
			continue

		# overclassification
		if len(target) < level: over += 1
	
		totalAcc = False

		# totally accurate
		if level >= 6 and len(target) >= 6:		
			if compare(token[5], target[5]):
				cnt_idx = 6
				totalAcc = True

		if totalAcc and level >= 7 and len(target) >= 7:
			if compare(token[6], target[6]):
				cnt_idx = 7
				totalAcc = True

		if not totalAcc:
			# accuracy
			cnt_idx = 0
			for t in range(min(len(target), level)):
				if compare(token[t], target[t]):
					cnt_idx = t + 1
				else:
					break

		for t in range(4, max(len(token), len(target))):
			# answer no info but predict --> False positive
			if ((t < len(target) and target[t] == '?') or t >= len(target)) and t < level:
				# False positive
				if token[t] not in PR_dict_list[t].keys():
					in_dict, in_dict_taxa = False, ''
					for d in PR_dict_list[t].keys():
						if compare(d, token[t]):
							in_dict = True
							in_dict_taxa = d
							break
					if not in_dict:
						PR_dict_list[t][token[t]] = [0,1,0] # TP, FP, FN
					else:
						PR_dict_list[t][in_dict_taxa][1] += 1
						
				else:
					PR_dict_list[t][token[t]][1] += 1

				continue

			if t <= cnt_idx:
				# True positive
				if target[t] not in PR_dict_list[t].keys():
					in_dict, in_dict_taxa = False, ''
					for d in PR_dict_list[t].keys():
						if compare(d, target[t]):
							in_dict = True
							in_dict_taxa = d
							break
					if not in_dict:
						PR_dict_list[t][target[t]] = [1,0,0] # TP, FP, FN
					else:
						PR_dict_list[t][in_dict_taxa][0] += 1
						
				else:
					PR_dict_list[t][target[t]][0] += 1

			else:
				# False negative
				if target[t] not in PR_dict_list[t].keys():
					in_dict, in_dict_taxa = False, ''
					for d in PR_dict_list[t].keys():
						if compare(d, target[t]):
							in_dict = True
							in_dict_taxa = d
							break
					if not in_dict:
						PR_dict_list[t][target[t]] = [0,0,1] # TP, FP, FN
					else:
						PR_dict_list[t][in_dict_taxa][2] += 1
						
				else:
					PR_dict_list[t][target[t]][2] += 1

				if t < level: # no prediction --> False negative only
					# False positive
					if token[t] not in PR_dict_list[t].keys():
						in_dict, in_dict_taxa = False, ''
						for d in PR_dict_list[t].keys():
							if compare(d, token[t]):
								in_dict = True
								in_dict_taxa = d
								break
						if not in_dict:
							PR_dict_list[t][token[t]] = [0,1,0] # TP, FP, FN
						else:
							PR_dict_list[t][in_dict_taxa][1] += 1
							
					else:
						PR_dict_list[t][token[t]][1] += 1



		for t in range(cnt_idx):
			if t == 6:
				correct_taxa.append([t+1, i, target[t-1] + ' ' + target[t], ';'.join(target[:t+1])])
			else:
				correct_taxa.append([t+1, i, target[t], ';'.join(target[:t+1])])
		
		for t in range(cnt_idx): correct_level[t] += 1
		
		acc_frac = cnt_idx/len(target)
		acc_frac_sum += acc_frac


	macroF1 = [0, 0, 0, 0, 0, 0, 0]
	weightedF1 = [0, 0, 0, 0, 0, 0, 0]
	P = [0, 0, 0, 0, 0, 0, 0]
	R = [0, 0, 0, 0, 0, 0, 0]
	for t in range(4, 7):
		weight_sum, count = 0, 0
		for tax in PR_dict_list[t].keys():
			TP, FP, FN = PR_dict_list[t][tax][0], PR_dict_list[t][tax][1], PR_dict_list[t][tax][2]
			weight = TP+FN
			if weight != 0: # Answer has this taxa
				r = TP/(TP+FN)
				R[t] += r
				
				if TP+FP != 0:
					p = TP/(TP+FP)
					P[t] += p
					if p+r != 0:
						macroF1[t] += 2*p*r/(p+r)
						weightedF1[t] += 2*p*r/(p+r)*weight
				weight_sum += weight
			else:
				count += 1
	
		P[t] /= (len(PR_dict_list[t].keys()) - count)
		R[t] /= (len(PR_dict_list[t].keys()) - count)
		macroF1[t] /= (len(PR_dict_list[t].keys()) - count)
		weightedF1[t] /= weight_sum
		print(t , len(PR_dict_list[t].keys()) - count)

	print('Amount:', eachLevel, len(y.keys()))
	print('Amount:', eachLevel / len(y.keys()))
	print('Accuracy: ', correct_level / len(y.keys()))

	return eachLevel, correct_level, unassigned, totalLevel, totalLevel, acc_frac_sum, over, correct_taxa, P, R, macroF1, weightedF1

def plot(name, to_plot, classifiers):
	name = output_name(name)

	from random import randrange
	labels = ['family', 'genus', 'species']
	colors = ['darkgreen', 'red', 'blue', 'chocolate', 'gold',  'purple', 'blueviolet',  'lightseagreen', 'navy', 'yellowgreen', 'lightcoral', 'deepskyblue']
	x = np.arange(0,len(labels)*0.5, 0.5)
	width = 0.1

	fig, ax = plt.subplots()
	
	print('--- plot Accuracy ---')
	print(len(to_plot))

	for i in range(len(to_plot)):
		ax.bar(x + (i-len(to_plot)/2)*width/2, to_plot[i]*100, width*9/20, label=output_name(classifiers[i]), color = colors[i])
	
	ax.set_ylabel('Accuracy (%)')

	plt.ylim(0, 100)
	plt.yticks(np.arange(0, 110, 10))
	ax.set_title(name)
	ax.set_xticks(x)
	ax.set_xticklabels(labels)
	lg = ax.legend(bbox_to_anchor=(1.05, 1.), loc='upper left', fontsize='medium')
	plt.savefig(name + '_plotAccuracy' + ('_algos' if 'sintax' in classifiers else '') + '.jpeg', dpi=400, format='jpeg', bbox_extra_artists=(lg,), bbox_inches='tight')

def plotPR(name, task):
	print(f'--- plot PR of {name} comparing {task} ---')
	import plotly.express as px
	from PIL import Image
	
	pd.options.plotting.backend = "plotly"
	
	df = pd.read_csv(name + ('_algos' if task == 'algos' else '') + '.csv', sep=',')
	fig = px.scatter(df, x="Precision", y="Recall", color=("DB" if task == 'db' else "Algorithm"), symbol="Level")
	fig.update_yaxes(range = [0,100])
	fig.update_xaxes(range = [1,100])
	fig.update_layout(legend = dict(font = dict(size = 10)), \
	                  legend_title = dict(font = dict(size = 12)), \
					  height=615,\
					  yaxis = dict(tickmode = 'linear', tick0 = 0, dtick = 10), \
					  xaxis = dict(tickmode = 'linear', tick0 = 0, dtick = 10))
	
	fig.write_image(name + '_' + task + "_PRscatter.jpeg", scale=5)
	im = Image.open(name + '_' + task + "_PRscatter.jpeg")
	im.save(name + '_' + task + "_PRscatter.jpeg", dpi = (400, 400))


def output_name(name):
	if name == 'itgdb_v1': lbl = 'Taxa_ITGDB'
	elif name == 'rsgidb': lbl = 'Seq_ITGDB'
	elif name == 'gg': lbl = 'Greengenes'
	elif name == 'slv': lbl = 'SILVA'
	elif name == 'udb': lbl = '16S-UDb'
	elif name == 'gtdb': lbl = 'GTDB'
	elif name == 'qiime': lbl = 'QIIME2'
	elif name == 'mothur': lbl = 'Mothur'
	elif name.startswith('udb'): lbl = '16S-UDb'
	elif name.startswith('exclude'): lbl = 'Exclusion'
	elif name.startswith('final_mock'): lbl = 'Mock'
	elif name.startswith('le2overlap'):lbl = 'Less-3-overlap'
	elif name.startswith('union'): lbl = 'Union'
	elif name.startswith('intersect'): lbl = 'Intersection'
	elif name.startswith('mock'): lbl = 'Mock'
	else: lbl = name.upper()
	
	return lbl


if __name__ == "__main__":
	''' python taxaAnalysis.py [subject]'''
	subject = sys.argv[1]
	if subject == 'sintax': sintax(sys.argv[2])
	elif subject == 'algos': algos(sys.argv[2])
	elif subject == 'plotPR': plotPR(sys.argv[2], sys.argv[3])
	else: print(f'Error: Unrecognized subject {subject} !')

