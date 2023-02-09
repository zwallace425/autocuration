# -*- coding: utf-8 -*-

# This object should handle the curation of an input suquence, utilizing the profile alignments,
# lookup table, and boundary file

import os
import re
import pandas as pd
import numpy as np
import subprocess
from Bio import SeqIO
from Blast import Blast
from InDelSubs import InDelSubs
from MolSeq import MolSeq


class Curation(object):

	# This object expects a single query nucleotide sequence in fasta at minimum. Optionally, it can also take a 
	# strain name denoted as [Speicies]_[Segment #]_[Subtype].  The object can also take in a specified boundary 
	# file, lookup table file, and a directory name specifying the location of all the profile fasta files.  If
	# those required files are not inputted as arguments, the object navigates to the default location of these
	# files.
	def __init__(self, query, strainName = None, boundaryFile = "profiles/Flu_profile_boundaries_20181012.txt", 
		lookupTable = "profiles/Flu_profile_lookupTable_20181203.txt", profile_dir = "profiles", output_dir = "outputs"):

		# Grab accession number and nucleotide sequence string for query sequence in the fasta file
		accession, sequence = self.get_acc_and_seq(query)

		# Only BLAST if strain name not passed to init
		if strainName:
			# Get the appropriate profile for the profile_dir based on the strain name
			profile = [filename for filename in os.listdir(profile_dir) if filename.startswith(strainName)][0]
		else:
			# BLAST query to determine the appropriate profile and strain name
			b = Blast(query)
			profile = b.get_profile()
			identity = b.get_identity()
			strainName = b.get_strain()

		# Compute the alignment of the query to the profile using MAFFT
		alignment = profile_dir+"/precomputed_alignment.fasta"
		cmd = "mafft --add "+query+" --globalpair --maxiterate 1000 "+profile_dir+"/"+profile+" > "+alignment
		subprocess.call(cmd, shell = True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

		# Parse boundary and lookup table text files
		boundary_df = self.parse_boundaryFile(strainName, boundaryFile)
		lookup_df = self.parse_lookupTable(strainName, lookupTable)

		muts = InDelSubs(alignment)
		
		# Get the Flags
		del_flags =	muts.deletion_flags(boundary_df, lookup_df)
		ins_flags = muts.insertion_flags(boundary_df)
		sub_flags = muts.substitution_flags(boundary_df)

		flags = del_flags + ins_flags + sub_flags

		# Determine the ambiguity flags, if any
		ambig_flags = []
		molseq = MolSeq(accession, sequence)
		if molseq.get_n_content() > 0.005:
			ambig_flags.append("Excess-N")
		if molseq.get_ambig_content() > 0.005:
			ambig_flags.append("Excess-Ambig")
		if identity < 0.5:
			ambig_flags.append("Excess-Dist")


		# Only save alignment if no insertion flags
		if len(ins_flags) == 0:
			self.save_alignment(accession, alignment, output_dir)

		self.flags = flags
		self.del_flags = del_flags
		self.ins_flags = ins_flags
		self.sub_flags = sub_flags
		self.ambig_flags = ambig_flags
		self.profile = profile
		self.strainName = strainName
		self.accession = accession

	# Return a table with all the NCR/CDS curation information about the sequence, otherwise return 'NO FLAGS'
	def curation_table(self):

		if self.flags == []:
			return("PASS")
		else:
			df = pd.DataFrame(self.flags, columns = ['Flag', 'Profile Position', 'Query Position', 'Variant', 'Length'])
			return(df)

	# Return a ambiguity flag(s) for input query sequence
	def ambiguity_flags(self):

		if self.ambig_flags == []:
			return("PASS")
		else:
			return(self.ambig_flags)

	# Update Dr. Macken's Table 6 for curation bookeeping
	def update_table6(self, Table6 = 'outputs/Table6_Jan2019Release.txt', ):

		table6 = pd.read_csv(Table6, sep = '\t')

		# If no flags, abort function and return table 6 back
		if self.flags == []:
			return

		# Loop through all detected flags
		for flag in self.flags:
			if (flag[0] == "5'NCR-ext") or (flag[0] == "3'NCR-ext"):
				table6_profile = table6[(table6['PROFILE_NAME'] == self.profile) & (table6['FLU_SUBTYPE'] == self.strainName) & 
				(table6['AUTO_ALIGNMENT_ISSUE'] == flag[0])]
				# Check if flag returns something already in the table
				if not table6_profile.empty:
					index = table6_profile.index[0]
					acc_list = set(str(table6['ACCESSION_LIST'][index]).split(","))
					if self.accession not in acc_list:
						acc_list.add(self.accession)
						acc_list = list(acc_list)
						acc_list.sort()
						acc_list = ",".join(acc_list)

						table6.at[index, 'STATUS'] = str("Updated")
						table6.at[index, 'ACCESSION_COUNT'] = int(table6['ACCESSION_COUNT'][index])+1
						table6.at[index, 'ACCESSION_LIST'] = str(acc_list)

				else:
					table6_profile = pd.DataFrame({'PROFILE_NAME': [self.profile], 'STATUS': ["New"], 'FLU_SUBTYPE': [self.strainName], 
						'AUTO_ALIGNMENT_ISSUE': [flag[0]], 'POS_PROFILE': [""], 'MUTATION_SUM'	: [""], 'ACCESSION_COUNT': [1],
						'ACCESSION_LIST': [self.accession]})
					table6 = pd.concat([table6, table6_profile], axis = 0)

			else:
				table6_profile = table6[(table6['PROFILE_NAME'] == self.profile) & (table6['FLU_SUBTYPE'] == self.strainName) & 
				(table6['AUTO_ALIGNMENT_ISSUE'] == flag[0]) & (table6['POS_PROFILE'] == flag[1])]
				if not table6_profile.empty:
					index = table6_profile.index[0]
					acc_list = set(str(table6['ACCESSION_LIST'][index]).split(","))
					if self.accession not in acc_list:
						acc_list.add(self.accession)
						acc_list = list(acc_list)
						acc_list.sort()
						acc_list = ",".join(acc_list)

						table6.at[index, 'STATUS'] = str("Updated")
						table6.at[index, 'ACCESSION_COUNT'] = int(table6['ACCESSION_COUNT'][index])+1
						table6.at[index, 'ACCESSION_LIST'] = str(acc_list)

						if flag[0] == "5'CTS-mut" or flag[0] == "3'CTS-mut":
							mut_sum = str(table6['MUTATION_SUM'][index])
							mut_sum = dict(item.split(":") for item in mut_sum.split(","))
							if mut_sum.get(flag[3]):
								mut_sum[flag[3]] = str(int(mut_sum[flag[3]])+1)
							else:
								mut_sum[flag[3]] = str(1)
							mut_sum = ",".join([key+':'+val for key, val in mut_sum.items()])
							table6.at[index, 'MUTATION_SUM'] = mut_sum
				
				else:
					mut_sum = ""
					if flag[0] == "5'CTS-mut" or flag[0] == "3'CTS-mut":
						mut_sum = flag[3]+':'+str(1)
					table6_profile = pd.DataFrame({'PROFILE_NAME': [self.profile], 'STATUS': ["New"], 'FLU_SUBTYPE': [self.strainName], 
						'AUTO_ALIGNMENT_ISSUE': [flag[0]], 'POS_PROFILE': [flag[1]], 'MUTATION_SUM'	: [mut_sum], 'ACCESSION_COUNT': [1],
						'ACCESSION_LIST': [self.accession]})
					table6 = pd.concat([table6, table6_profile], axis = 0)

		table6 = table6.sort_values(by = ['PROFILE_NAME', 'ACCESSION_COUNT'], ascending = [True, False]).reset_index(drop=True)
		table6.to_csv(Table6, sep = '\t', index = False)


	# Return just a table of the deletion flags, if any
	def deletion_flags(self):

		if self.del_flags == []:
			return("PASS")
		else:
			df = pd.DataFrame(self.del_flags, columns = ['Flag', 'Profile Position', 'Query Position', 'Variant', 'Length'])
			return(df)

	# Return just a table of the insertion flags, if any
	def insertion_flags(self):

		if self.ins_flags == []:
			return("PASS")
		else:
			df = pd.DataFrame(self.ins_flags, columns = ['Flag', 'Profile Position', 'Query Position', 'Variant', 'Length'])
			return(df)

	# Return just a table of the substitution flags, if any
	def substitution_flags(self):

		if self.sub_flags == []:
			return("PASS")
		else:
			df = pd.DataFrame(self.sub_flags, columns = ['Flag', 'Profile Position', 'Query Position', 'Variant', 'Length'])
			return(df)

	# Return the file name of the profile
	def get_profile(self):

		return(self.profile)
	
	# Return the strain name of the virus denoted as [Speicies]_[Segment #]_[Subtype]
	def get_strain(self):

		return(self.strainName)

	# Return the accession number of the query sequence
	def get_accession(self):

		return(self.accession)

	# Meant for grabbing accession number and nucleotide sequence string from the query sequence fasta file
	@staticmethod
	def get_acc_and_seq(query):
		
		for seq_record in SeqIO.parse(query, 'fasta'):
			metadata = str(seq_record.id)
			sequence = str(seq_record.seq).strip()
		
		acc = re.split(r'[^a-zA-Z0-9\s\w]+', metadata)
		if re.search(r'[a-zA-Z]+', acc[0]) and re.search(r'[0-9]+', acc[0]):
			accession = acc[0]
		else:
			accession = acc[1]

		return(accession, sequence)

	# Meant for parsing the boundary file to know the CTS, NCR, and CDS start and end regions of the profile
	@staticmethod
	def parse_boundaryFile(strainName, boundaryFile):

		boundary_types = open(boundaryFile, 'r').readlines()
		boundary = ""
		for line in boundary_types:	
			if (line.startswith(strainName)):
				boundary = line
				break

		boundary = boundary.split("|")
		boundary.pop(0)
		boundaries = {}
		for coord in boundary:
			feat = coord.split("=")[0]
			pos = int(coord.split("=")[1])
			boundaries[feat] = pos

		boundary_df = pd.DataFrame({"Region": ["CTS5", "NCR5", "CDS", "NCR3", "CTS3"],
			"Start": [boundaries['START'], boundaries['CTS5']+1, boundaries['ATG'], boundaries['STOP']+1, boundaries['CTS3']],
			"End": [boundaries['CTS5'], boundaries['ATG']-1, boundaries['STOP'], boundaries['CTS3']-1, boundaries['END']]})		

		return(boundary_df)

	# Meant to parse the lookup table so that we can use information for a specific strain.
	@staticmethod
	def parse_lookupTable(strainName, lookupTable):

		lookup_open = open(lookupTable, 'r', encoding = "ISO-8859-1")
		lookup_allowed = lookup_open.readlines()
		strain_lookup = []
		for line in lookup_allowed:	
			if (line.startswith(strainName)):
				line_attribs = line.split('\t')
				flag = line_attribs[1]
				start_end = line_attribs[2].split("..")
				if (len(start_end) == 2):
					start = int(start_end[0])
					end = int(start_end[1])
				else:
					start = int(start_end[0])
					end = int(start_end[0])
				lookup_df = pd.DataFrame({'Flag': [flag], 'Start': [start], 'End': [end]})
				strain_lookup.append(lookup_df)

		lookup_df = pd.concat(strain_lookup, ignore_index = True)

		return(lookup_df)

	# Save the mafft alignment if there were no insertions
	@staticmethod
	def save_alignment(query_acc, alignment, output_dir):

		with open(alignment, 'r') as file1:
			with open(output_dir+"/"+query_acc+"_aligned.fasta", 'w') as file2:
				for line in file1:
					file2.write(line)
