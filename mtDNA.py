from Bio.Seq import Seq
from Bio.SeqUtils import nt_search

class MitochondrialDna: 
    def __init__ (self, sequence:str): 
        self._sequence = Seq(sequence.upper())
    
    def extract_seq(self, start, end):
        '''Extracts subsequences from genomic sequences'''
        if start < 1 or end > len(self._sequence) or start > end:
            print ("Invalid start or end position.")
        return self._sequence[start-1:end]
        

    def gc_content(self, sequence = None): 
        '''Calculates the GC content in a given sequence'''

        seq = sequence if sequence else self._sequence
        g = seq.count('G')
        c = seq.count('C')
        gc_percentage = (g+c)/len(seq)*100
        return {'G count: ', g, 
                'C count: ', c,
                'GC Percentage: ',gc_percentage}
    
    def seq_len(self, sequence = None): 
        '''Calculates the length of a given sequence'''

        seq = sequence if sequence else self._sequence
        return len(seq)


class GenomicMotif (MitochondrialDna): 
    '''Represents sequence motifs'''

    def __init__ (self, motif:str, sequence:str):
        super().__init__(sequence)
        self.motif = motif.upper()
        
    def search_motif(self): 
        '''Searches for motifs within mitochondrial DNA'''
        positions = nt_search(str(self._sequence),self.motif)
        if len(positions) <= 1:
            return 'No motifs found in the selected sequence'
        else:
            motif_count = len(positions[1:]) #Counts motif occurrences
            seq_len = len(self._sequence)
            distribution = (motif_count/seq_len)*100 if seq_len > 0 else 0 #analyse motif distribution across sequences
            return positions[1:]
            return f"Motif count: {motif_count}\n Percentage distribution: {distribution}"
