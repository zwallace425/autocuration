# This object is meant for BLASTing the query sequence against the profiles so that we can 
# ulimately deterimine the strain name ([Species]_[Segment]_[Subtype]) of the Influenza query
# sequence and us that to choose the profile and guide the rest of auto-curation.

from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

class Blast(object):

	# Intialize the blasting object by passing in the query sequence and optional blast database.
	# The blast database should be pre-computed and stored somewhere, otherwise, a blast database
	# can be passed into the init function.  The default parameters are the default locations of the
	# blast database and blast command line program.
	def __init__(self, query, blast_db = 'blast/flu_profiles_db.fasta', blastn_cmd = 'blast/blastn', 
		blast_result = 'test_files/blast_result.txt'):

		# Blast query against database of profile sequences
		cmdline = NcbiblastnCommandline(cmd=blastn_cmd, query=query, db=blast_db, outfmt=5, out=blast_result)
		stdout, stderr = cmdline()

		result_handle = open(blast_result)
		for result in NCBIXML.parse(result_handle):
			if len(result.alignments) > 0:
				profile = result.alignments[0].title.split("|")[3]
				identity = float(result.alignments[0].hsps[0].identities)/float(result.alignments[0].hsps[0].align_length)
			else:
				profile = "Unknown"
				identity = 0
			break

		self.profile = profile
		self.identity = identity
	
	# Return the profile name mapping to the query sequence
	def get_profile(self):
		
		return(self.profile)

	# Return the percent identity of the query blast result
	def get_identity(self):

		return(self.identity)

	# Return the strain name of the query
	def get_strain(self):

		if self.profile == "Unknown":
			strain = "Unknown"
		else:
			strain = self.profile.split("_")
			strain = strain[:len(strain)-1]
			strain = "_".join(strain)

		return(strain)




       
                    
