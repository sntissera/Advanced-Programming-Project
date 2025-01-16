from Bio.Seq import Seq
from Bio.SeqUtils import nt_search

class MitochondrialDna: 
    def __init__ (self, sequence:str): 
        self._sequence = Seq(sequence.upper())
    
    def extract_seq(self, start, end):
        '''Extracts subsequences from genomic sequences'''
        if start < 1 or end > len(self._sequence) or start > end:
            return 'Invalid start or end position.' 
         
        return self._sequence[start-1:end]
        

    def gc_content(self, sequence = None): 
       '''Calculates the GC content in a given sequence'''
        seq = sequence if sequence != None else self._sequence
        valid_bases = ['G', 'C', 'A', 'T']
        filtered_seq = [base for base in seq if base in valid_bases]
        filtered_length = len(filtered_seq)
    
        if filtered_length == 0:
            return "Invalid sequence or no valid bases found."

        g = seq.count('G')
        c = seq.count('C')
        gc_percentage = ((g+c)/filtered_length)*100
        return f"{round(gc_percentage,2)} %"
    
    def seq_len(self, sequence = None):
        '''Calculates the length of a given sequence'''
        seq = sequence if sequence != None else self._sequence
        valid_bases = ['G', 'C', 'A', 'T']
        filtered_seq = [base for base in seq if base in valid_bases]
        filtered_length = len(filtered_seq)
    
        if filtered_length == 0:
            return "Invalid sequence or no valid bases found."

        return f"{filtered_length} bp"


class GenomicMotif (MitochondrialDna): 
    '''Represents sequence motifs'''

    def __init__ (self, motif:str, sequence:str):
        super().__init__(sequence)
        self.motif = motif.upper()
        
    def search_motif(self): 
        '''Searches for motifs within mitochondrial DNA'''
        sequence_str = str(self._sequence)
        positions = nt_search(sequence_str, self.motif)
        
        if len(positions) <= 1:
            return 'No motifs found in the selected sequence'
        else:
            # motif_count = len(positions[1:]) #Counts motif occurrences
            seq_len = len(sequence_str)
            # distribution = (motif_count/seq_len)*100 if seq_len > 0 else 0 #analyse motif distribution across sequences
            location = positions[1:]
            return location
