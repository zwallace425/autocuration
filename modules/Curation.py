# -*- coding: utf-8 -*-

# This object should handle the curation of an input suquence, utilizing the profile alignments,
# lookup table, and boundary file

import os
import re
import pandas as pd
import numpy as np
import subprocess
from Bio import SeqIO
from datetime import datetime
from Blast import Blast
from InDelSubs import InDelSubs
from MolSeq import MolSeq


class Curation(object):

	# This object expects a SINGLE QUERY nucleotide sequence in fasta at minimum. Optionally, it can also take a 
	# strain name denoted as [Speicies]_[Segment #]_[Subtype].  The object can also take in a specified boundary 
	# file, lookup table file, and a directory name specifying the location of all the profile fasta files.  If
	# those required files are not inputted as arguments, the object navigates to the default location of these
	# files.
	def __init__(self, query, strain_name = None, boundary_file = "profiles/Flu_profile_boundaries_20181012.txt", 
		lookup_table = "profiles/Flu_profile_lookupTable_20181203.txt", profile_dir = "profiles", output_dir = "outputs"):
			
		# Grab accession number and nucleotide sequence string for query sequence in the fasta file
		accession, sequence = self.get_acc_and_seq(query)

		# Get the profile. Only BLAST if strain name not passed to init
		if strain_name:
			# Get the appropriate profile for the profile_dir based on the strain name
			try:
				profile = [filename for filename in os.listdir(profile_dir) if filename.startswith(strain_name)][0]
			except:
				raise Exception("\nERROR: Invalid profiles directory or strain name\n")
		else:
			# BLAST query to determine the appropriate profile and strain name
			b = Blast(query)
			profile = b.get_profile()
			identity = b.get_identity()
			strain_name = b.get_strain()

		# Check if BLAST result found profile. If not, abort the rest of the pipeline
		if profile == "Unknown":
			self.profile = "Unknown"
			self.strain_name = "Unknown"
			self.accession = accession
			return

		# Parse boundary and lookup table text files
		boundary_df = self.parse_boundary_file(strain_name, boundary_file)
		lookup_df = self.parse_lookup_table(profile, lookup_table)

		# Determine the ambiguity flags, if any
		ambig_flags = []
		molseq = MolSeq(accession, sequence)
		if molseq.get_n_content() > 0.005:
			ambig_flags.append("Excess-N")
		if molseq.get_ambig_content() > 0.005:
			ambig_flags.append("Excess-Ambig")
		if identity < 0.95:
			ambig_flags.append("Excess-Dist")

		# Compute the alignment of the query to the profile using MUSCLE
		alignment = profile_dir+"/precomputed_alignment.fasta"
		cmd = './muscle -maxiters 1000 -profile -in1 '+profile_dir+'/'+profile+' -in2 '+query+' -out '+alignment
		subprocess.call(cmd, shell = True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

		# Initialize InDelSubs object to identify mutations from the alignment
		muts = InDelSubs(alignment)
		
		# Get the Flags
		del_flags =	muts.deletion_flags(boundary_df, lookup_df)
		ins_flags = muts.insertion_flags(boundary_df)
		sub_flags = muts.substitution_flags(boundary_df)

		mut_flags = del_flags + ins_flags + sub_flags

		# Only save alignment if no insertion flags
		if len(ins_flags) == 0:
			self.save_alignment(accession, alignment, output_dir)

		self.mut_flags = mut_flags
		self.del_flags = del_flags
		self.ins_flags = ins_flags
		self.sub_flags = sub_flags
		self.ambig_flags = ambig_flags
		self.profile = profile
		self.strain_name = strain_name
		self.accession = accession

	# Return a table with all mutation flags occuring in the sequence, otherwise return 'Pass' if no flags
	# However, if no profile alignment was found in BLAST step, return 'Unknown'
	def mutation_flags(self):

		if self.profile == "Unknown":
			return("Unknown")

		if self.mut_flags == []:
			return("Pass")
		else:
			df = pd.DataFrame(self.mut_flags, columns = ['Flag', 'Profile Position', 'Query Position', 'Variant', 'Length'])
			return(df)

	# Return a ambiguity flag(s) for input query sequence
	# If no profile alignment was found in BLAST step, return 'Excess-Dist'
	def ambiguity_flags(self):

		if self.profile == "Unknown":
			return("Excess-Dist")

		if self.ambig_flags == []:
			return("Pass")
		else:
			return(self.ambig_flags)

	# Update Dr. Macken's Table 6 for curation bookeeping
	def update_table6(self, Table6 = 'outputs/Table6_Jan2019Release.txt'):

		# If no profile was found from BLAST (no strong hits), abort function
		# There is nothing to update for table 6
		if self.profile == "Unknown":
			return

		# If no flags, abort function
		# There is nothing to update for table 6
		if self.mut_flags == []:
			return

		try:
			table6 = pd.read_csv(Table6, sep = '\t')
		except:
			raise Exception("\nERROR: Invalid table6 directory and/or file\n")

		# Revert STATUS_THIS_MONTH to Unchanged for those that were updated last month
		# Revert CURRENT_MONTH_INCREASE count to zero if new month approached
		# Update PAST_MONTH_INCREASE count to the value of CURRENT_MONTH_INCREASE
		today = datetime.today()
		current_month = datetime.today().replace(day=1)
		table6['LAST_UPDATED'] = pd.to_datetime(table6['LAST_UPDATED'], errors='coerce').dt.floor('d')
		past_months = table6[(table6['LAST_UPDATED'] < current_month) & (table6['STATUS_THIS_MONTH'] != "Unchanged")]
		for ind in past_months.index:
			table6.at[ind, 'STATUS_THIS_MONTH'] = "Unchanged"
			table6.at[ind, 'PAST_MONTH_INCREASE'] = int(table6['CURRENT_MONTH_INCREASE'][ind])
			table6.at[ind, 'CURRENT_MONTH_INCREASE'] = 0

		# Loop through all detected flags
		for flag in self.mut_flags:
			if (flag[0] == "5'NCR-ext") or (flag[0] == "3'NCR-ext"):
				table6_profile = table6[(table6['PROFILE_NAME'] == self.profile) & (table6['FLU_SUBTYPE'] == self.strain_name) & 
				(table6['AUTO_ALIGNMENT_ISSUE'] == flag[0])]
				# Check if flag returns something already in the table
				if not table6_profile.empty:
					index = table6_profile.index[0]
					status = str(table6['STATUS_THIS_MONTH'][index])
					acc_list = set(str(table6['ACCESSION_LIST'][index]).split(","))
					if self.accession not in acc_list:
						acc_list.add(self.accession)
						acc_list = list(acc_list)
						acc_list.sort()
						acc_list = ",".join(acc_list)

						if status == "Unchanged":
							table6.at[index, 'STATUS_THIS_MONTH'] = str("Updated")
							table6.at[index, 'LAST_UPDATED'] = today
						table6.at[index, 'ACCESSION_TOTAL'] = int(table6['ACCESSION_TOTAL'][index])+1
						table6.at[index, 'CURRENT_MONTH_INCREASE'] = int(table6['CURRENT_MONTH_INCREASE'][index])+1
						table6.at[index, 'ACCESSION_LIST'] = str(acc_list)

				else:
					table6_profile = pd.DataFrame({'PROFILE_NAME': [self.profile], 'STATUS_THIS_MONTH': ["New"], 'LAST_UPDATED': [today], 
						'FLU_SUBTYPE': [self.strain_name], 'AUTO_ALIGNMENT_ISSUE': [flag[0]], 'POS_PROFILE': [""], 'MUTATION_SUM': [""], 
						'ACCESSION_TOTAL': [1], 'CURRENT_MONTH_INCREASE': [1], 'ACCESSION_LIST': [self.accession]})
					table6 = pd.concat([table6, table6_profile], axis = 0)

			else:
				table6_profile = table6[(table6['PROFILE_NAME'] == self.profile) & (table6['FLU_SUBTYPE'] == self.strain_name) & 
				(table6['AUTO_ALIGNMENT_ISSUE'] == flag[0]) & (table6['POS_PROFILE'] == flag[1])]
				if not table6_profile.empty:
					index = table6_profile.index[0]
					status = str(table6['STATUS_THIS_MONTH'][index])
					acc_list = set(str(table6['ACCESSION_LIST'][index]).split(","))
					if self.accession not in acc_list:
						acc_list.add(self.accession)
						acc_list = list(acc_list)
						acc_list.sort()
						acc_list = ",".join(acc_list)

						if status == "Unchanged":
							table6.at[index, 'STATUS_THIS_MONTH'] = str("Updated")
							table6.at[index, 'LAST_UPDATED'] = today
						table6.at[index, 'ACCESSION_TOTAL'] = int(table6['ACCESSION_TOTAL'][index])+1
						table6.at[index, 'CURRENT_MONTH_INCREASE'] = int(table6['CURRENT_MONTH_INCREASE'][index])+1
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
					table6_profile = pd.DataFrame({'PROFILE_NAME': [self.profile], 'STATUS_THIS_MONTH': ["New"], 'LAST_UPDATED': [today], 
						'FLU_SUBTYPE': [self.strain_name], 'AUTO_ALIGNMENT_ISSUE': [flag[0]], 'POS_PROFILE': [flag[1]], 'MUTATION_SUM': [mut_sum], 
						'ACCESSION_TOTAL': [1], 'CURRENT_MONTH_INCREASE': [1], 'ACCESSION_LIST': [self.accession]})
					table6 = pd.concat([table6, table6_profile], axis = 0)

		table6 = table6.sort_values(by = ['PROFILE_NAME', 'ACCESSION_TOTAL'], ascending = [True, False]).reset_index(drop=True)
		table6.to_csv(Table6, sep = '\t', index = False)


	# Return just a table of the deletion flags, if any
	def deletion_flags(self):

		if self.profile == "Unknown":
			return("Unknown")

		if self.del_flags == []:
			return("Pass")
		else:
			df = pd.DataFrame(self.del_flags, columns = ['Flag', 'Profile Position', 'Query Position', 'Variant', 'Length'])
			return(df)

	# Return just a table of the insertion flags, if any
	def insertion_flags(self):

		if self.profile == "Unknown":
			return("Unknown")

		if self.ins_flags == []:
			return("Pass")
		else:
			df = pd.DataFrame(self.ins_flags, columns = ['Flag', 'Profile Position', 'Query Position', 'Variant', 'Length'])
			return(df)

	# Return just a table of the substitution flags, if any
	def substitution_flags(self):

		if self.profile == "Unknown":
			return("Unknown")

		if self.sub_flags == []:
			return("Pass")
		else:
			df = pd.DataFrame(self.sub_flags, columns = ['Flag', 'Profile Position', 'Query Position', 'Variant', 'Length'])
			return(df)

	# Return the file name of the profile
	def get_profile(self):

		return(self.profile)
	
	# Return the strain name of the virus denoted as [Speicies]_[Segment #]_[Subtype]
	def get_strain(self):

		return(self.strain_name)

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
	def parse_boundary_file(strain_name, boundary_file):

		try:
			boundary_types = open(boundary_file, 'r').readlines()
		except:
			raise Exception("\nERROR: Invalid boundary_file directory and/or file\n")
		boundary = ""
		for line in boundary_types:	
			if (line.startswith(strain_name)):
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
	def parse_lookup_table(profile, lookup_table):

		try:
			lookup_open = open(lookup_table, 'r', encoding = "ISO-8859-1")
		except:
			raise Exception("\nERROR: Invalid lookup_table directory and/or file\n")
		lookup_allowed = lookup_open.readlines()
		profile_lookup = []
		for line in lookup_allowed:	
			if (line.startswith(profile)):
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
				profile_lookup.append(lookup_df)

		if profile_lookup != []:
			lookup_df = pd.concat(profile_lookup, ignore_index = True)
		else:
			lookup_df = pd.DataFrame({'Flag': ['XXX'], 'Start': [-1], 'End': [-1]})

		return(lookup_df)

	# Save the MUSCLE alignment if there were no insertions
	@staticmethod
	def save_alignment(query_acc, alignment, output_dir):

		with open(alignment, 'r') as file1:
			with open(output_dir+"/"+query_acc+"_aligned.fasta", 'w') as file2:
				for line in file1:
					file2.write(line)