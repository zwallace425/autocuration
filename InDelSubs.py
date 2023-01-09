# This object is meant for detecting the positions of any deletions in the query sequence
# upon alignment to the profile

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

		self.nkdp = nkdp
		self.nkip = nkip
		self.nkadp = nkadp
		self.query_seq = query_seq
		self.profile_seqs = profile_seqs


	# Reports the flags triggered by deletions. Requires the process boundary file and
	# lookup table file.
	def deletion_flags(self, boundary_df, lookup_df):

		# Extract deletions not accepted in the profile alignment
		flag_dels = list(set(set(self.nkdp) - set(self.nkadp)))

		# Group together non-keep-length sequential deletions
		positions = []
		flag_dels = np.array(flag_dels)
		flag_dels = np.sort(flag_dels)
		for k, g in groupby(enumerate(flag_dels), lambda ix: ix[0] - ix[1]):
			positions.append(list(map(itemgetter(1), g)))

		flags = []
		for pos in positions:
			
			profile_del = [(j-len([i for i in self.nkip if i < j]))+1 for j in pos]
			query_del = [(j-len([i for i in self.nkdp if i <= j]))+1 for j in pos]

			region = boundary_df[(boundary_df['Start'] <= profile_del[0]) & (boundary_df['End'] >= profile_del[len(profile_del)-1])]
			region = list(region['Region'])[0]

			if len(pos) > 1:
				profile_pos = str(profile_del[0])+".."+str(profile_del[len(profile_del)-1])
			else:
				profile_pos = str(profile_del[0])
			query_pos = str(query_del[0])+".."+str(query_del[0]+1)

			# Go through each of the possible flags
			if (region == 'CTS5'):
				flags.append(tuple(("5'CTS-del", profile_pos, query_pos, "del", len(pos))))
			elif (region == 'CTS3'):
				flags.append(tuple(("3'CTS-del", profile_pos, query_pos, "del", len(pos))))
			elif (region == 'NCR5'):
				if InDelSubs.check_lookup(lookup_df, "5'NCR-del", profile_del[0], profile_del[len(profile_del)-1]) == "FAIL":
					flags.append(tuple(("5'NCR-del", profile_pos, query_pos, "del", len(pos))))
			elif (region == 'NCR3'):
				if InDelSubs.check_lookup(lookup_df, "3'NCR-del", profile_del[0], profile_del[len(profile_del)-1]) == "FAIL":
					flags.append(tuple(("3'NCR-del", profile_pos, query_pos, "del", len(pos))))
			elif (region == 'CDS'):
				if (len(pos) % 3) == 0:
					if InDelSubs.check_lookup(lookup_df, "CDS-3Xdel", profile_del[0], profile_del[len(profile_del)-1]) == "FAIL":
						flags.append(tuple(("CDS-3Xdel", profile_pos, query_pos, "del", len(pos))))
				else:
					if InDelSubs.check_lookup(lookup_df, "CDS-del", profile_del[0], profile_del[len(profile_del)-1]) == "FAIL":
						flags.append(tuple(("CDS-del", profile_pos, query_pos, "del", len(pos))))

		return flags

	# Reports the flags triggered by insertions. Requires the boundary file.
	def insertion_flags(self, boundary_df):

		profile_length = boundary_df['End'][4]

		# Group together non-keep-length sequential insertions
		positions = []
		nkip = np.array(list(self.nkip))
		nkip = np.sort(nkip)
		for k, g in groupby(enumerate(nkip), lambda ix: ix[0] - ix[1]):
			positions.append(list(map(itemgetter(1), g)))

		flags = []
		for pos in positions:

			ins_muts = self.query_seq[pos[0]:pos[len(pos)-1]+1]

			profile_ins = [(j-len([i for i in self.nkip if i <= j]))+1 for j in pos]
			query_ins = [(j-len([i for i in self.nkdp if i < j]))+1 for j in pos]

			if len(pos) > 1:
				query_pos = str(query_ins[0])+".."+str(query_ins[len(query_ins)-1])
			else:
				query_pos = str(query_ins[0])
			profile_pos = str(profile_ins[0])+".."+str(profile_ins[0]+1)

			if (profile_ins[0] == 0) or (profile_ins[0] == 1) :
				flags.append(tuple(("5'NCR-ext", profile_pos, query_pos, ins_muts, len(pos))))
				continue
			if (profile_ins[0] == profile_length) or (profile_ins[0] == profile_length-1):
				flags.append(tuple(("3'NCR-ext", profile_pos, query_pos, ins_muts, len(pos))))
				continue

			region = boundary_df[(boundary_df['Start'] <= profile_ins[0]) & (boundary_df['End'] >= profile_ins[0])]
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

	# Reports the flags triggered by CTS substitutions. Requires the boundary file
	def substitution_flags(self, boundary_df):

		CTS5 = boundary_df[(boundary_df['Region'] == 'CTS5')]
		CTS3 = boundary_df[(boundary_df['Region'] == 'CTS3')]

		# Get 5'/3' CTS start/end and adjust zero basing
		CTS5_start = list(CTS5['Start'])[0] - 1
		CTS3_start = list(CTS3['Start'])[0] - 1
		CTS5_end = list(CTS5['End'])[0] - 1
		CTS3_end = list(CTS3['End'])[0] - 1

		# Adjust the 5'/3' CTS start/end to accommidate for the non-keep length alignment
		# by factoring in the nkip values
		CTS5_start_adj = CTS5_start + len([i for i in self.nkip if i < CTS5_start])
		CTS3_start_adj = CTS3_start + len([i for i in self.nkip if i < CTS3_start])
		CTS5_end_adj = CTS5_end + len([i for i in self.nkip if i < CTS5_end])
		CTS3_end_adj = CTS3_end + len([i for i in self.nkip if i < CTS3_end])

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
			subs = self.query_seq[pos[0]:pos[len(pos)-1]+1]
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
			subs = self.query_seq[pos[0]:pos[len(pos)-1]+1]
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







