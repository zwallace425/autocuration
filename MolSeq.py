# This is code is modified from the originally developed version by Christian Zmasek at JCVI.  It's used to handle a 
# query sequence, identify ambiguity in a sequence (indeterminate nucleotides or irregular characters), and write sequences
# back to a fasta

class MolSeq(object):

    # This object expects some the accession number of a sequence and the nucleotide sequence as a string of NAs.
    # Sequence in fasta in not supported for this object.
    def __init__(self, accession, sequence):
        
        self.seq_id = str(accession).strip()
        self.seq = str(sequence).strip()

    def get_seq_id(self):
        
        return self.seq_id

    def set_seq_id(self, id):
        
        self.seq_id = id

    def get_seq(self):
        
        return self.seq

    def get_length(self):
        
        return len(self.seq)

    def to_fasta(self):
        
        return ">{}\n{}".format(self.get_seq_id(), self.get_seq())

    # Counting the regular nucleotides of the nucelic acid sequence
    def count_regular_nucs(self):
        
        a = self.get_seq().count('a') + self.get_seq().count('A')
        c = self.get_seq().count('c') + self.get_seq().count('C')
        g = self.get_seq().count('g') + self.get_seq().count('G')
        t = self.get_seq().count('t') + self.get_seq().count('T')
        return a + c + g + t

    # Count the indeterminate, N's, of the NA sequence
    def count_indeterminate_nucs(self):
        
        return self.get_seq().count('n') + self.get_seq().count('N')

    # Return the percentage of Ns in the sequence
    def get_n_content(self):
        
        return (self.count_indeterminate_nucs() / self.get_length())

    # Return the percentage of ambiguity in the sequence
    def get_ambig_content(self):
        
        return ((self.get_length() - (self.count_regular_nucs() + self.count_indeterminate_nucs())) /
                self.get_length())
