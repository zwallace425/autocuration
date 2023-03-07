# This program is meant for testing accuracy of the autocuration pipeline

# The program requires the input FASTA file to run the test on and the name
# of the output CSV file

import sys
import pandas as pd
import numpy as np
import Autocuration as ac
from Bio import SeqIO

test_file = sys.argv[1]
results = sys.argv[2]

count = 1
evaluation = []
for seq_record in SeqIO.parse(test_file, 'fasta'):
	seq_desc = str(seq_record.description)
	seq = str(seq_record.seq)

	actual = seq_desc.split("Alignment_flag:")[1]

	# No need to evaluate ambiguity sequences, very trivial, should always work
	if actual in ['HammDist', 'Excess-N', 'Excess-Ambig']:
		continue
	
	seq_fasta = ac.MolSeq(seq_desc, seq).to_fasta()
	with open('query.fasta', 'w') as f:	f.write(seq_fasta)
	curation = ac.Curation('query.fasta')
	
	accession = curation.get_accession()
	my_result = curation.mutation_flags()
	if not isinstance(my_result, pd.DataFrame):
		df = pd.DataFrame({"Accession": [accession], "Actual Flag": [actual], "My Flag": [my_result]})
	else:
		my_result = set(list(my_result['Flag']))
		if actual in my_result:
			df = pd.DataFrame({"Accession": [accession], "Actual Flag": [actual], "My Flag": [actual]})
		else:
			my_result = list(my_result)[0]
			df = pd.DataFrame({"Accession": [accession], "Actual Flag": [actual], "My Flag": [my_result]})

	evaluation.append(df)

	print("Done with seq:", count)
	count += 1

evaluation = pd.concat(evaluation, ignore_index = True)
evaluation.to_csv(results)


