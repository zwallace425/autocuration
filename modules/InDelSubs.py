# This object is meant for detecting the positions of any insertions, deletions, deletions, or substitutions,
# in the query sequence upon alignment to the profile and then determines the artifact flag  due to
# those mutations

import pandas as pd
import numpy as np
from Bio import SeqIO
from itertools import groupby
from operator import itemgetter
from collections import Counter


class InDelSubs(object):

	# This object just expects the profile alignment in fasta, computed in a separate module.
	# This alignment is not assumed to be have preserved length of the original profile because
	# it needs to identify potential insertions.  Hence, the primary purpose of this init object
	# is to grab ahold of "non-keep-length" insertions and deletion positions for later computing
	# all mutations.
	def __init__(self, alignment):

		# Get all the sequences from the alignment
		sequences = []
		for seq_record in SeqIO.parse(alignment, 'fasta'):
			seq = str(seq_record.seq)
			sequences.append(seq)

		query_seq = sequences[len(sequences)-1]
		profile_seqs = sequences[:len(sequences)-1]

		# Deletion positions, denoted as gaps, within the query sequence.
		# These positions are with respect to the Non-Keeplength alginment
		# nkdp = Non-Keeplength Deletion Positions
		nkdp = set([i for i, ltr in enumerate(query_seq) if ltr == "-"])

		first_seq_dashes = set([i for i, ltr in enumerate(profile_seqs[0]) if ltr == "-"])
		profile_seqs_minus = profile_seqs[1:len(profile_seqs)]

		profile_dash_union = first_seq_dashes
		for seq in profile_seqs_minus:
			seq_dashes = set([i for i, ltr in enumerate(seq) if ltr == "-"])
			profile_dash_union = profile_dash_union.union(seq_dashes)


		profile_dash_intsct = first_seq_dashes
		for seq in profile_seqs_minus:
			seq_dashes = set([i for i, ltr in enumerate(seq) if ltr == "-"])
			profile_dash_intsct = profile_dash_intsct.intersection(seq_dashes)

		# Insertion positions, denoted as gaps conserved accross entire profile but
		# not in the query.  These positions are with respect to the Non-Keeplength alignment
		# nkip = Non-Keeplength Insertion Positions
		nkip = set(profile_dash_intsct - nkdp)

		# Accepted deletion positions are gaps/dashes that exist in the profile sequences.
		# These positions are with respect to the Non-Keeplength alignment.
		# nkadp = Non-Keeplenght Accepted Deletion Positions
		nkadp = set(profile_dash_union - nkip)

		# Group together non-keep-length sequential deletions as a list of lists
		del_groups = []
		flag_dels = list(set(set(nkdp) - set(nkadp)))
		flag_dels = np.array(flag_dels)
		flag_dels = np.sort(flag_dels)
		for k, g in groupby(enumerate(flag_dels), lambda ix: ix[0] - ix[1]):
			del_groups.append(list(map(itemgetter(1), g)))

		# Group together non-keep-length sequential insertions as a list of lists
		ins_groups = []
		nkip = np.array(list(nkip))
		nkip = np.sort(nkip)
		for k, g in groupby(enumerate(nkip), lambda ix: ix[0] - ix[1]):
			ins_groups.append(list(map(itemgetter(1), g)))

		self.nkdp = nkdp
		self.nkip = nkip
		self.nkadp = nkadp
		self.query_seq = query_seq
		self.profile_seqs = profile_seqs
		self.del_groups = del_groups
		self.ins_groups = ins_groups

	# Reports the flags triggered by deletions. Requires the processed boundary file and
	# lookup table file as pandas dataframes.
	def deletion_flags(self, boundary_df, lookup_df):

		# Needed for adjusting profile positions to query positions
		all_profile_ins = [(j-len([i for i in self.nkip if i <= j]))+1 for j in self.nkip]
		all_profile_del = [(j-len([i for i in self.nkip if i < j]))+1 for j in self.nkdp]

		flags = []
		for pos in self.del_groups:
			
			# Profile and query deletions in the group. profile_del has exact profile positions
			# but query_del has the positions in the query proceeding deletions event
			profile_del = [(j-len([i for i in self.nkip if i < j]))+1 for j in pos]
			query_del = [(j-len([i for i in self.nkdp if i <= j]))+1 for j in pos]

			# This will get ALL regions for which the deletion gaps may be overlapping
			regions = boundary_df[(boundary_df['Start'] <= profile_del[len(profile_del)-1]) & (boundary_df['End'] >= profile_del[0])]
			regions = regions.to_dict('records')
			
			# Loop through all possible regions for which gaps may be overlapping
			for reg in regions:
				region = reg['Region']
				
				# Annotated region start end with respect to profile
				start_p = reg['Start']
				end_p = reg['End']
				
				# Annotated region start end with respect to query
				# Take into account insertions and deletions occuring before the profile position
				start_q = start_p + len([i for i in all_profile_ins if i < start_p]) - len([j for j in all_profile_del if j <= start_p])
				end_q = end_p + len([i for i in all_profile_ins if i < end_p]) - len([j for j in all_profile_del if j <= end_p])

				# Positions of deletions in the annotated region with respect to profile and query
				region_dels_p = [i for i in profile_del if start_p <= i <= end_p]
				region_dels_q = [i for i in query_del if start_q <= i <= end_q]

				# This will not allow flagging of sequences that are simply shorter;
				# ie, gaps at the beginning or end of sequence relative to profile
				if (region_dels_q[0] == 0) or (region_dels_q[0] == (len(self.query_seq)-len(self.nkdp))):
					continue

				if len(pos) > 1:
					profile_pos = str(region_dels_p[0])+".."+str(region_dels_p[len(region_dels_p)-1])
				else:
					profile_pos = str(region_dels_p[0])
				query_pos = str(region_dels_q[0])+".."+str(region_dels_q[0]+1)

				# Go through each of the possible flags
				if (region == 'CTS5'):
					flags.append(tuple(("5'CTS-del", profile_pos, query_pos, "del", len(region_dels_p))))
				elif (region == 'CTS3'):
					flags.append(tuple(("3'CTS-del", profile_pos, query_pos, "del", len(region_dels_p))))
				elif (region == 'NCR5'):
					if self.check_lookup(lookup_df, "5'NCR-del", region_dels_p[0], region_dels_p[len(region_dels_p)-1]) == "FAIL":
						flags.append(tuple(("5'NCR-del", profile_pos, query_pos, "del", len(region_dels_p))))
				elif (region == 'NCR3'):
					if self.check_lookup(lookup_df, "3'NCR-del", region_dels_p[0], region_dels_p[len(region_dels_p)-1]) == "FAIL":
						flags.append(tuple(("3'NCR-del", profile_pos, query_pos, "del", len(region_dels_p))))
				elif (region == 'CDS'):
					if (len(pos) % 3) == 0:
						if self.check_lookup(lookup_df, "CDS-3Xdel", region_dels_p[0], region_dels_p[len(region_dels_p)-1]) == "FAIL":
							flags.append(tuple(("CDS-3Xdel", profile_pos, query_pos, "del", len(region_dels_p))))
					else:
						if self.check_lookup(lookup_df, "CDS-del", region_dels_p[0], region_dels_p[len(region_dels_p)-1]) == "FAIL":
							flags.append(tuple(("CDS-del", profile_pos, query_pos, "del", len(region_dels_p))))

		return flags

	# Reports the flags triggered by insertions. Requires the boundary file as pandas dataframe.
	def insertion_flags(self, boundary_df):

		profile_length = boundary_df['End'][4]

		flags = []
		for pos in self.ins_groups:

			ins_muts = self.query_seq[pos[0]:pos[len(pos)-1]+1].upper()

			profile_ins = [(j-len([i for i in self.nkip if i <= j]))+1 for j in pos]
			query_ins = [(j-len([i for i in self.nkdp if i < j]))+1 for j in pos]

			if len(pos) > 1:
				query_pos = str(query_ins[0])+".."+str(query_ins[len(query_ins)-1])
			else:
				query_pos = str(query_ins[0])
			profile_pos = str(profile_ins[0])+".."+str(profile_ins[0]+1)

			if (profile_ins[0] == 0):
				flags.append(tuple(("5'NCR-ext", profile_pos, query_pos, ins_muts, len(pos))))
				continue
			if (profile_ins[0] == profile_length):
				profile_pos = str(profile_ins[0])+".."
				flags.append(tuple(("3'NCR-ext", profile_pos, query_pos, ins_muts, len(pos))))
				continue

			region = boundary_df[(boundary_df['Start'] <= profile_ins[len(profile_ins)-1]) & (boundary_df['End'] >= profile_ins[0])]
			region = list(region['Region'])[0]

			# Go through each of the possible flags
			if (region == 'CTS5'):
				flags.append(tuple(("5'CTS-ins", profile_pos, query_pos, ins_muts, len(pos))))
			elif (region == 'CTS3'):
				flags.append(tuple(("3'CTS-ins", profile_pos, query_pos, ins_muts, len(pos))))
			elif (region == 'NCR5'):
				flags.append(tuple(("5'NCR-ins", profile_pos, query_pos, ins_muts, len(pos))))
			elif (region == 'NCR3'):
				flags.append(tuple(("3'NCR-ins", profile_pos, query_pos, ins_muts, len(pos))))
			elif (region == 'CDS'):
				if (len(pos) % 3) == 0:
					flags.append(tuple(("CDS-3Xins", profile_pos, query_pos, ins_muts, len(pos))))
				else:
					flags.append(tuple(("CDS-ins", profile_pos, query_pos, ins_muts, len(pos))))

		return flags

	# Reports the flags triggered by CTS substitutions. Requires the boundary file as pandas dataframe.
	def substitution_flags(self, boundary_df):

		CTS5 = boundary_df[(boundary_df['Region'] == 'CTS5')]
		CTS3 = boundary_df[(boundary_df['Region'] == 'CTS3')]

		# Get 5'/3' CTS start/end and adjust zero basing
		CTS5_start = list(CTS5['Start'])[0] - 1
		CTS3_start = list(CTS3['Start'])[0] - 1
		CTS5_end = list(CTS5['End'])[0] - 1
		CTS3_end = list(CTS3['End'])[0] - 1

		# Adjust the 5'/3' CTS start/end to accommidate for the non-keep length alignment
		# by factoring in the positions of sequential deletions from the ins_groups
		CTS5_start_adj = CTS5_start
		CTS3_start_adj = CTS3_start + sum([len(self.ins_groups[i]) for i in range(len(self.ins_groups)) if self.ins_groups[i][0] < CTS3_start])
		CTS5_end_adj = CTS5_end + sum([len(self.ins_groups[i]) for i in range(len(self.ins_groups)) if self.ins_groups[i][0] < CTS5_end])
		CTS3_end_adj = CTS3_end + sum([len(self.ins_groups[i]) for i in range(len(self.ins_groups)) if self.ins_groups[i][0] < CTS3_end])

		# Get the CTS regions of the query sequence
		CTS5_query = self.query_seq[CTS5_start_adj:CTS5_end_adj+1]
		CTS3_query = self.query_seq[CTS3_start_adj:CTS3_end_adj+1]

		cts5_muts = [i for i in range(len(CTS5_query)) if CTS5_query[i] not in [seq[CTS5_start_adj:CTS5_end_adj+1][i]
			for seq in self.profile_seqs] and not ((i in self.nkip) or (i in self.nkdp))]

		cts3_muts = [i+CTS3_start_adj for i in range(len(CTS3_query)) if CTS3_query[i] not in [seq[CTS3_start_adj:CTS3_end_adj+1][i] 
			for seq in self.profile_seqs] and not ((i+CTS3_start_adj in self.nkip) or (i+CTS3_start_adj in self.nkdp))]

		cts5_positions = []
		for k, g in groupby(enumerate(cts5_muts), lambda ix: ix[0] - ix[1]):
			cts5_positions.append(list(map(itemgetter(1), g)))

		cts3_positions = []
		for k, g in groupby(enumerate(cts3_muts), lambda ix: ix[0] - ix[1]):
			cts3_positions.append(list(map(itemgetter(1), g)))

		flags = []

		#Check for 5'CTS mutations:
		for pos in cts5_positions:
			subs = self.query_seq[pos[0]:pos[len(pos)-1]+1].upper()
			# Check if the current substituions are just N characters, continue loop if so and do not flag
			if len(set(list(subs))) == 1 and list(set(list(subs)))[0] == 'N':
				continue
			profile_sub = [(j-len([i for i in self.nkip if i < j]))+1 for j in pos]
			query_sub = [(j-len([i for i in self.nkdp if i < j]))+1 for j in pos]

			if len(pos) > 1:
				profile_pos = str(profile_sub[0])+".."+str(profile_sub[len(profile_sub)-1])
				query_pos = str(query_sub[0])+".."+str(query_sub[len(query_sub)-1])
			else:
				profile_pos = str(profile_sub[0])
				query_pos = str(query_sub[0])

			flags.append(tuple(("5'CTS-mut", profile_pos, query_pos, subs, len(pos))))

		#Check for 3'CTS mutations:
		for pos in cts3_positions:
			subs = self.query_seq[pos[0]:pos[len(pos)-1]+1].upper()
			if len(set(list(subs))) == 1 and list(set(list(subs)))[0] == 'N':
				continue
			profile_sub = [(j-len([i for i in self.nkip if i < j]))+1 for j in pos]
			query_sub = [(j-len([i for i in self.nkdp if i < j]))+1 for j in pos]

			if len(pos) > 1:
				profile_pos = str(profile_sub[0])+".."+str(profile_sub[len(profile_sub)-1])
				query_pos = str(query_sub[0])+".."+str(query_sub[len(query_sub)-1])
			else:
				profile_pos = str(profile_sub[0])
				query_pos = str(query_sub[0])

			flags.append(tuple(("3'CTS-mut", profile_pos, query_pos, subs, len(pos))))

		return flags

	# Check the lookup table to determine if the flag gets accepted.
	@staticmethod
	def check_lookup(lookup_df, flag, start, end):

		table = lookup_df[(lookup_df['Flag'] == flag) & (lookup_df['Start'] <= start) & (lookup_df['End'] >= end)]

		if len(table) > 0:
			return("PASS")
		else:
			return("FAIL")