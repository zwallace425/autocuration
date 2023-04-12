# Autocuration of Influenza Sequences and Saving Pre-Computed Alignments

This is a pipeline for curating Influenza genome sequences and designating artifact flags. The artifact
flags are meant for alerting descrepancies that have been discovered in the genome sequence and triggering
the need for correcting the sequencing.  The methods for discovering artifact flags in the seqeunce were
developed by JCVI/BV-BRC in collaboration with Dr. Catherine Macken at the Unversity of Aukland.

## Pipeline Dependencies

(1) python >=3.x is required

(2) MAFFT is a required.  Please download MAFFT version 7.490 from https://mafft.cbrc.jp/alignment/software/
and place in your absolute PATH or use conda install by running
	
	conda install -c bioconda mafft

(3) MUSCLE can be used in lieu of the MAFFT alignment.  MUSCLE exits in this repo as an executable file
and will get called if the pipeline runs MUSCLE instead of MAFFT.  However, this excecutable has been
downloaded for Mac OSX 64 bit.  To download MUSCLE v3.8.31 for the appropriate OS,see that downloads 
at https://drive5.com/muscle/downloads_v3.htm and follow instructions at
https://drive5.com/muscle/manual/install.html

(3) To install the rest of the necessary dependencies, run the following in command line
	
	source env_setup.sh

## Running the Pipeline and Generating Outputs

After navigating the cloned directory, the prefered way to run the pipeline via the command line is with

	python Autocuration.py --query [Influenza FASTA sequence(s)]

The above command will run the pipeline using MAFFT incorporating iterative refinement with both WSP
and consistency scores.  To run the pipeline using MUSCLE, execute the following command

	python Autocuration.py --query [Influenza FASTA sequence(s)] --align_alg muscle  

In running those commands, the pipeline will output a curation report with designated artifact flags,
if any exist.  Additionally, running the pipeline will automatically update Dr. Macken's 'Table 6',
which is a dynamically growing table book keeping all past and present autocuration flag results.
That original table is found in the 'outputs' folder as Table6_Jan2019Release.txt but gets updated
with runs of the pipeline and left in the same folder.  

Along with autocuration, a key component of this pipeline is saving a pre-computed alignment of the 
inputted query sequence.  This alignment from MAFFT or MUSCLE will get saved ONLY IF the sequence had 
no insertions. This pre-computed alignment is saved in 'outputs' as 'ACCESSION_aligned.fasta'.

## Profile Alignments, Lookup Table, and Boundary File

The pipeline depends on the Influenza subtype-specific profile alignments for ultimately curating
the sequences, as query sequences get aligned to the profile to determine what the invalid mutations
are.  These profiles are curated alignments specific to the influenza subtype and are all found in the
'profiles' folder. Additionally, a 'lookup table' and 'boundary file' are key files used in the curation
process, the former for serving as a lookup dictionary for valid deletions and the latter marking the 
start and end loctions of CTS, NCR, and CDS regions of the profile.  The profile alignments, lookup table,
and boundary file are all in the 'profiles' folder.  If these need to get updated, which will very likely
be the case, they need to kept into the same folder with the older versions removed.

## BLAST Database

A key step in this pipeline is determining the subtype of the query sequence (Species/Type/Segment), and
this is done by BLASTing all incoming query sequences against a database of all sequences making up the
profiles.  This BLAST database of profile sequences in stored in the 'blast' folder, along with the BLAST
command line executables for making a BLAST database and running the BLAST job.  This said, whenever the
profile alignments are updated, the BLAST database MUST be updated.

To recreate a BLAST database due to a new collection of profile alignments in the profiles folder, run

	python build_blast_db.py

## Pipeline Performance

Disjoint sequence batches for Influenza A-D tagged with artifact flags in the FASTA metadata from the previous 
version of the Autocuration pipeline were downloaded from the legacy Influenza Research Database.  The
sequences from these four Influenza FASTA files were ran against our pipeline in 'performance/test_accuracy.py', 
and the curation results were used to analyze precision, recall, and accuracy of our version of the Autocuration 
pipeline relative to the results from the previous pipeline.  There are several precision, recall, and accuracy 
metrics that need to used to evaluate this pipeline, and a thorough evaluation of this pipeline performance with 
each of these metrics can be viewed in the notebook 'AnalyzeAccuracy.ipynb'.

## Pipeline Broken Down by Modules

The pipeline can also be broken into indvidual scriptable modules useful as pluggins. These are all
located in the 'modules' folder.  Note, all of these modules are just python classes and they are all
aggregated into one script in 'Autocuration.py'.

main.py - used for running all the modules as pipeline

Curation.py - primary module for curating sequence and updating Table 6

InDelSubs.py - disovers all the mutation artifact flags induced by insertions, deletions, or substitutions

Blast.py - used for BLASTing query sequence against profile sequence to determine subtype of the sequence

MolSeq.py - used for counting indeterminant nucleotides or irregular characters, computing ambiguity



