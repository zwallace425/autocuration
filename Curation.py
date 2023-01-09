# -*- coding: utf-8 -*-

# This object should handle the curation of an input suquence, utilizing the profile alignments,
# lookup table, and boundary file

import os
import pandas as pd
import numpy as np
import subprocess
from Blast import Blast
from InDelSubs import InDelSubs


class Curation(object):

	# This object expecta a query nucleotide sequence in fasta and a strain name denoted as [Speicies]_[Segment #]_[Subtype]
	# at minimum.  The object can also take in a specified boundary file, lookup table file, and a directory name specifying
	# the location of all the profile fasta files.
	def __init__(self, query, strainName = None, boundaryFile = "profiles/Flu_profile_boundaries_20181012.txt", 
		lookupTable = "profiles/Flu_profile_lookupTable_20181203.txt", profile_dir = "profiles", output_dir = "outputs"):

		# Only BLAST if strain name not passed to init
		if strainName:
			# Get the appropriate profile for the profile_dir based on the strain name
			profile = [filename for filename in os.listdir(profile_dir) if filename.startswith(strainName)][0]
		else:
			# BLAST query to determine the appropriate profile and strain name
			b = Blast(query)
			profile = b.get_profile()
			strainName = b.get_strain()

		# Compute the alignment of the query to the profile using MAFFT
		alignment = profile_dir+"/precomputed_alignment.fasta"
		cmd = "mafft --add "+query+" --globalpair --maxiterate 1000 "+profile_dir+"/"+profile+" > "+alignment
		subprocess.call(cmd, shell = True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

		# Parse boundary and lookup table text files
		boundary_df = Curation.parse_boundaryFile(strainName, boundaryFile)
		lookup_df = Curation.parse_lookupTable(strainName, lookupTable)

		muts = InDelSubs(alignment)
		
		# Get the Flags
		del_flags =	muts.deletion_flags(boundary_df, lookup_df)
		ins_flags = muts.insertion_flags(boundary_df)
		sub_flags = muts.substitution_flags(boundary_df)

		flags = del_flags + ins_flags + sub_flags

		#if len(ins_flags) > 0:
			# Save the pre-computed alignment to a output_dir folder only if there were no insertions

		self.flags = flags
		self.del_flags = del_flags
		self.ins_flags = ins_flags
		self.sub_flags = sub_flags
		self.profile = profile
		self.strainName = strainName


	# Return a table with all the curation information about the sequence, otherwise return 'NO FLAGS'
	def curation_table(self):

		if self.flags == []:
			return("NO FLAGS")
		else:
			df = pd.DataFrame(self.flags, columns = ['Flag', 'Profile Position', 'Query Position', 'Variant', 'Length'])
			return(df)

	# Update Dr. Macken's Table 6 for curation bookeeping
	def update_table6(self, table6 = 'outputs/Table6_Jan2019Release.txt'):

		print("Code to update table6")

	# Return just a table of the deletion flags, if any
	def deletion_flags(self):

		if self.del_flags == []:
			return("NO DELETION FLAGS")
		else:
			df = pd.DataFrame(self.del_flags, columns = ['Flag', 'Profile Position', 'Query Position', 'Variant', 'Length'])
			return(df)

	# Return just a table of the insertion flags, if any
	def insertion_flags(self):

		if self.ins_flags == []:
			return("NO INSERTION FLAGS")
		else:
			df = pd.DataFrame(self.ins_flags, columns = ['Flag', 'Profile Position', 'Query Position', 'Variant', 'Length'])
			return(df)

	# Return just a table of the substitution flags, if any
	def substitution_flags(self):

		if self.sub_flags == []:
			return("NO SUBSTITUTION FLAGS")
		else:
			df = pd.DataFrame(self.sub_flags, columns = ['Flag', 'Profile Position', 'Query Position', 'Variant', 'Length'])
			return(df)
		

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
