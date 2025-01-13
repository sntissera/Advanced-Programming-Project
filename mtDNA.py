from Bio.Seq import Seq
from Bio.SeqUtils import nt_search

class MitochondrialDna: 
    def __init__ (self, sequence:str): 
        self._sequence = Seq(sequence)
    
    def extract_seq(self, start, end):
        '''Extracts subsequences from genomic sequences'''
        subseq = self._sequence[start-1:end]
        return subseq

    def gc_content(self, sequence = None): 
        '''Calculates the GC content in a given sequence'''

        seq = sequence if sequence else self._sequence
        g = seq.count('G')
        c = seq.count('C')
        print('G Count: ',g)
        print('C count: ',c)
        print('GC Percentage: ', ((g+c)/len(seq))*100)
        return f"{g}\n{c}\n{((g+c)/len(seq))*100}"
    
    def seq_len(self, sequence = None): 
        '''Calculates the length of a given sequence'''

        seq = sequence if sequence else self._sequence
        print('Sequence length: ', len(seq))
        return len(seq)


class GenomicMotif (MitochondrialDna): 
    '''Represents sequence motifs'''

    def __init__ (self, motif:str, sequence:str):
        super().__init__(sequence)
        self.motif = motif
        
    def search_motif(self): 
        '''Searches for motifs within mitochondrial DNA'''
        position = nt_search(str(self._sequence),self.motif)
        if len(position) <= 1:
            raise ValueError ('No motifs found in the selected sequence')
        return position

    def count_motif(self): 
        '''Counts motif occurrences and analyse their distribution across sequences'''
        position = nt_search(str(self._sequence),self.motif)
        motif_count = len(position[1:])
        seqlen = len(self._sequence)
        
        print ('Motif count: ', motif_count)
        print('Distribution through out the sequence(%): ', motif_count/seqlen*100)
        return f"{motif_count}\n{motif_count/seqlen*100}"




