import time
import copy


def seq_IntegrateTwo(candidate_dict, basis_dict):
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

def saveFile(Seq_ItgDB, out):
	Seq_ItgDBFasta = open( out + '_seq.fasta', 'w')
	Seq_ItgDBTxt = open( out + '_taxa.txt', 'w')

	for i in Seq_ItgDB.keys():
		for j in Seq_ItgDB[i]:
			ids = '>' + j[0]
			Seq_ItgDBFasta.write(ids+'\n')
			Seq_ItgDBFasta.write(j[1]+'\n')
			Seq_ItgDBTxt.write(ids[1:]+'\t'+i+'\n')

	Seq_ItgDBFasta.close()
	Seq_ItgDBTxt.close()
