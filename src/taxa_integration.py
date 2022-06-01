import copy
import multiprocessing as mp
import math
import os

def taxa_IntegrateTwo(dict_basis, dict_add, s_basis, s_add, n_threads):	#dict_add is the candidate database
	print("Integrating database: ", s_basis, "+", s_add)

	l_add = []	#l_add = [[s_ID, sd_microbime_obj],.....]
	for mkey,mval in dict_add.items():
		l_add.append([mkey,mval])

		
	n_group_mic_num = math.floor(len(l_add) / n_threads)	
	
	tasks = []
	
	for i in range(0, n_threads):	
		if(i < n_threads):
			tasks.append((l_add[i*(n_group_mic_num):((i+1)*n_group_mic_num)-1], dict_basis, s_basis, s_add))
		else:
			tasks.append((l_add[i*(n_group_mic_num):len(l_add)], dict_basis, s_basis, s_add))
			
	with mp.Pool(processes=n_threads) as pool:	
		
		try:
			l_results = pool.starmap(Integrate_content, tasks)	
			
			set_no_collect = set()	
			for rIter in l_results:
				for sIter in rIter:
					set_no_collect.add(sIter)	
			
			for mkey,mval in dict_add.items():
				if(mkey in set_no_collect):
					continue
				else:
					dict_basis[mkey] = mval
			
			return dict_basis
			
		except KeyboardInterrupt:
			print('got ctrl-C while pool mapping, terminating the pool')
			pool.terminate()
			print('pool is terminated')
	

def Integrate_content(l_add_partial, dict_basis, s_basis, s_add):	
	
	dict_cp_basis = copy.deepcopy(dict_basis)
	
	l_no_collect = []	
	
	for iter in l_add_partial:
		s_add_taxa = ""

		if(s_add == "slv"):	#The format of species name in SILVA is "genus name" + "species name".
			l_gs = iter[1].l_taxa[6].split('_')
			if(len(l_gs) < 2):
				s_add_taxa = ';'.join(iter[1].l_taxa)
			else:
				s_add_taxa = ';'.join(iter[1].l_taxa[0:4]) + ";" + l_gs[0] + ";" + l_gs[1]
		else:
			s_add_taxa = ';'.join(iter[1].l_taxa)
		
		b_is_collect = True
		
		for mkey,mval in dict_cp_basis.items():
			s_basis_taxa = ';'.join(mval.l_taxa)
			
			if(s_add_taxa.lower() == s_basis_taxa.lower()):
				if(iter[1].s_seq == mval.s_seq):
					b_is_collect = False
					break
					
		if(b_is_collect == True):
			dict_basis[iter[0]] = iter[1]
		else:
			l_no_collect.append(iter[0])
	
	return l_no_collect
	
def OutputITGDB(dict_itgdb, out):
	s_projectPath = os.environ["HOME"] + "/Project/itgdb/"

	wps = open(out + "_seq.fasta", 'w', encoding = "utf-8")
	wpt = open(out + "_tax.txt", 'w', encoding = "utf-8")

	for mkey,mval in dict_itgdb.items():
		wps.write(">" + mkey + "\n")
		wps.write(mval.s_seq + "\n")
		
		s_taxa = "k__" + mval.l_taxa[0] + "; p__" + mval.l_taxa[1] + "; c__" + mval.l_taxa[2] + "; o__" + mval.l_taxa[3] + "; f__" + mval.l_taxa[4] + "; g__" + mval.l_taxa[5] + "; s__" + mval.l_taxa[6]
		wpt.write(mkey + "\t" + s_taxa + "\n")
	

	wps.close()
	wpt.close()