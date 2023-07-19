# Main program for running the whole script from commandline

import sys
import time
import argparse
from Bio import SeqIO
from Curation import Curation
from MolSeq import MolSeq


# Required argument: --query [QUERY FASTA]
# Optional argument: --flag [muts/ambig/ins/del/sub] (ie, the type of flags to return)
# Optional argument: --table6 [TABLE6 FILE PATH]
if __name__ == "__main__":

	parser = argparse.ArgumentParser()

	# Required argument
	parser.add_argument('--query', dest = 'query', type = str)
	
	# Optional arguments
	parser.add_argument('--flag', dest = 'flag', type = str)
	parser.add_argument('--table6', dest = 'table6', type = str)
	args = parser.parse_args()

	if (not args.query):
		sys.exit("\nERROR: No query sequence input\n")
	if (args.flag and (args.flag != 'mut' and args.flag != 'ambig' and args.flag != 'ins' and args.flag != 'del' and args.flag != 'sub')):
		sys.exit("\nERROR: Invalid flag argument\n --flag [all/ambig/ins/del/sub]\n")
	if (args.table6):
		table6 = args.table6
	else:
		table6 = 'outputs/Table6_Jan2019Release.txt'

	for seq_record in SeqIO.parse(args.query, 'fasta'):
		seq_id = str(seq_record.id)
		seq = str(seq_record.seq)
		seq_fasta = MolSeq(seq_id, seq).to_fasta()
		with open('query.fasta', 'w') as f:	f.write(seq_fasta)
		start = time.time()
		cur = Curation('query.fasta')
		if not args.flag:
			print("Accession:", cur.get_accession())
			print("Subtype:", cur.get_strain())
			print("Summary Flag:", cur.summary_flag())
			print("Ambiguity Flags:", cur.ambiguity_flags())
			print("Mutation Flags:\n", cur.mutation_flags())
			print('\n')
		elif args.flag == 'mut':
			print("Accession:", cur.get_accession())
			print("Subtype:", cur.get_strain())
			print("Mutation Flags:\n", cur.mutation_flags())
			print('\n')
		elif args.flag == 'ambig':
			print("Accession:", cur.get_accession())
			print("Subtype:", cur.get_strain())
			print("Summary Flag:", cur.summary_flag())
			print("Ambiguity Flags:", cur.ambiguity_flags())
			print('\n')
		elif args.flag == 'ins':
			print("Accession:", cur.get_accession())
			print("Subtype:", cur.get_strain())
			print("Summary Flag:", cur.summary_flag())
			print("Insertion Flags:\n", cur.insertion_flags())
			print('\n')
		elif args.flag == 'del':
			print("Accession:", cur.get_accession())
			print("Subtype:", cur.get_strain())
			print("Summary Flag:", cur.summary_flag())
			print("Deletion Flags:\n", cur.deletion_flags())
			print('\n')
		elif args.flag == 'sub':
			print("Accession:", cur.get_accession())
			print("Subtype:", cur.get_strain())
			print("Summary Flag:", cur.summary_flag())
			print("Substitution Flags:\n", cur.substitution_flags())
			print('\n')
		cur.update_table6(Table6 = table6)
		end = time.time()
		print("Compute time for sequence:", round(end - start, 3))
		print('\n')