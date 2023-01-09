# This script is used for the Blast database of influenza profile sequence.  We already have a collection of profile
# alignments in a folder called 'profiles', but we need to use profile alignments in this folder and the sequences that
# make up each alignment in order to build a comprehensive profile blast database for which we can run an in coming
# query sequence against

import os
import pandas as pd
import numpy as np
import subprocess
from Bio import SeqIO


sequences = set()
seq_dict = {}
for file in os.listdir('profiles'):
	if file.endswith('.afa') or file.endswith('.fasta'):
		profile_name = file
		for seq_record in SeqIO.parse('profiles/'+file, 'fasta'):
			seq = str(seq_record.seq)
			seq = seq.replace("-", "")
			seq = seq.replace("~", "")
			seq_id = str(seq_record.id)
			seq_name = seq_id+"|"+profile_name
			if not seq_dict.get(seq_name):
				if not seq in sequences:
					seq_dict[seq_name] = seq
					sequences.add(seq)


outfile = open("blast/flu_profiles_db.fasta", "w")

for key in seq_dict.keys():
	outfile.write(">" + key + "\n" + seq_dict[key] + "\n")

outfile.close()


			
