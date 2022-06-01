import os
import argparse
from utils import *
from seq_Integration import *
from taxa_Integration import *


if __name__ == "__main__":
	task = sys.argv[1]
	if task != 'taxa' and task != 'seq': print('Invalid task. Please enter \"seq\" or \"taxa\" for sequence or taxonomy-based integration.')

	parser = argparse.ArgumentParser()
	parser.add_argument('--ggSeq'	,	type=str  , default='../gg/99_otus.fasta')
	parser.add_argument('--ggTaxa'	,	type=str  , default='../gg/99_otu_taxonomy.txt')
	parser.add_argument('--slvSeq'	,	type=str  , default='../slv/silva-138-99-seqs.fasta')
	parser.add_argument('--slvTaxa'	,	type=str  , default='../slv/silva-138-99-tax.txt')
	parser.add_argument('--rdpSeq'	,	type=str  , default='../rdp/RefOTUs.fasta')
	parser.add_argument('--rdpTaxa'	,	type=str  , default='../rdp/Ref_taxonomy.txt')
	parser.add_argument('--out'		, 	type=str  , default='ItgDB')
	args = parser.parse_args()

	Output_name = task + '_' + args.out

	if task == 'seq':
		dict_gg  = Preprocess(fastaFile=args.ggSeq, txtFile=args.ggTaxa)
		dict_slv = Preprocess(fastaFile=args.slvSeq, txtFile=args.slvTaxa)
		dict_rdp = Preprocess(fastaFile=args.rdpSeq, txtFile=args.rdpTaxa)

		# Change this order if needed
		addin_slv, replace_slv, dict_rsidb = seq_IntegrateTwo(dict_slv, dict_rdp)
		addin_gg, replace_gg, dict_seq_itgdb = seq_IntegrateTwo(dict_gg, dict_rsidb)
		
		print(f'add {addin_slv} from SILVA, add {addin_gg} from Greengenes')
		print(f'replace {replace_slv} with SILVA, replace {replace_gg} with Greengenes')
		
		saveFile(dict_seq_itgdb, out = Output_name)

	elif task == 'taxa':
		N_THREADS = 48

		dict_rdp = ParseFile(args.rdpSeq, args.rdpTaxa, "rdp")
		dict_slv = ParseFile(args.slvSeq, args.slvTaxa, "slv")
		dict_gg  = ParseFile(args.ggSeq, args.ggTaxa, "gg")

		#Integrate RDP, SILVA, and Greengenes databases.
		print(len(dict_rdp))
		dict_itgdb = taxa_IntegrateTwo(dict_rdp, dict_gg, "rdp", "gg", n_threads)
		print(len(dict_itgdb))
		dict_taxa_itgdb = taxa_IntegrateTwo(dict_itgdb, dict_slv, "itgdb", "slv", n_threads)
		print(len(dict_taxa_itgdb))

		OutputITGDB(dict_taxa_itgdb, out = Output_name)


