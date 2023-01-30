# This object is meant for BLASTing the query sequence against the profiles so that we can 
# ulimately deterimine the strain name ([Species]_[Segment]_[Subtype]) of the Influenza query
# sequence and us that to choose the profile and guide the rest of auto-curation.

from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

class Blast(object):

	# Intialize the blasting object by passing in the query sequence and optional blast database.
	# The blast database should be pre-computed and stored somewhere, otherwise, a blast database
	# can be passed into the init function.
	def __init__(self, query, blast_db = 'blast/flu_profiles_db.fasta', blastn_cmd = 'blast/blastn', 
		blast_result = 'blast/blast_result.xml'):

		# Blast query against database of profile sequences
		cmdline = NcbiblastnCommandline(cmd=blastn_cmd, query=query, db=blast_db, outfmt=5, out=blast_result)
		stdout, stderr = cmdline()

		result_handle = open(blast_result)
		profile = [blast_result.alignments[0].title for blast_result in NCBIXML.parse(result_handle)][0]

		profile = profile.split("|")[3]

		self.profile = profile
	
	# Return the profile name mapping to the query sequence
	def get_profile(self):
		
		return(self.profile)

	# Return the strain name of the query
	def get_strain(self):

		strain = self.profile.split("_")
		strain = strain[:len(strain)-1]
		strain = "_".join(strain)

		return(strain)




       
                    
