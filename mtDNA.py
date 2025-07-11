from Bio.Seq import Seq
from Bio.SeqUtils import nt_search
 
class MitochondrialDna: 
    def __init__ (self, sequence:str): 
        self._sequence = Seq(sequence)
    
    def extract_seq(self, start, end):
        '''Extracts subsequences from genomic sequences'''
        if start < 1 or end > len(self._sequence) or start > end:
            return ("Invalid start or end position.")
        sub_sequence =  str(self._sequence[start-1:end])
        modified_seq = '\n'.join([sub_sequence[i:i+10] for i in range (0, len(sub_sequence), 10)])
        return modified_seq
        

    def gc_content(self, sequence = None): 
        '''Calculates the GC content in a given sequence'''

        seq = sequence if sequence != None else self._sequence
        g = seq.count('G')
        c = seq.count('C')
        gc_percentage = (g+c)/len(seq)*100
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

    def __init__ (self, sequence:str, motif:str):
        super().__init__(sequence)
        self.set_motif(motif)

    def set_motif(self, motif):
        for i in motif: 
            if i not in ['G', 'C', 'A', 'T']:
                raise ValueError("Motif not valid")
            else: 
                self.motif = motif

    def search_motif(self): 
        '''Searches for motifs within mitochondrial DNA'''
        sequence_str = str(self._sequence)
        positions = nt_search(sequence_str, self.motif)
        
        if len(positions) <= 1:
            return 'No motifs found in the selected sequence'
        else:
            return positions[1:]

    def distribution(self):
        '''analyse motif distribution across sequences'''
        
        sequence_str = str(self._sequence)
        positions = nt_search(sequence_str, self.motif)
        motif_count = len(positions[1:]) 
        seq_len = len(sequence_str)
        distribution = (motif_count/seq_len)*100 if seq_len > 0 else 0 
        return f'{round(distribution,2)} %'

 
