from Bio.Seq import Seq
from Bio.SeqUtils import nt_search

class mitochondrial_dna: 
    def __init__ (self, sequence:str): 
        self._sequence = Seq(sequence)
    
    def extract_seq(self):
        '''Extracts subsequences from genomic sequences'''
        start = int(input('Position of the first nucleotide of the sequence: '))
        end = int(input('Position of the last nucleotide of the sequence: '))
        subseq = self._sequence[start-1:end]
        return (subseq)

    def gc_content(self, sequence = None): 
        '''Calculates the GC content in a given sequence'''

        seq = sequence if sequence else self._sequence
        g = seq.count('G')
        c = seq.count('C')
        print('G Count: ',g)
        print('C count: ',c)
        print('GC Percentage: ', ((g+c)/len(seq))*100)
    
    def seq_len(self, sequence = None): 
        '''Calculates the length of a given sequence'''

        seq = sequence if sequence else self._sequence
        print('Sequence length: ', len(seq))


class genomic_motif (mitochondrial_dna): 
    '''Represents sequence motifs'''

    def __init__ (self, motif:str, sequence:str):
        super().__init__(sequence)
        self.motif = motif
        
    def search_motif(self): 
        '''Searches for motifs within mitochondrial DNA'''
        position = nt_search(str(self._sequence),self.motif)
        if len(position) <= 1:
            raise ValueError ('No motifs found in the selected sequence')
        return(position)

    def count_motif(self): 
        '''Counts motif occurrences and analyse their distribution across sequences'''
        position = nt_search(str(self._sequence),self.motif)
        motif_count = len(position[1:])
        seqlen = len(self._sequence)
        
        print ('Motif count: ', motif_count)
        print('Distribution through out the sequence(%): ', motif_count/seqlen*100)


#For testing the code    
seq1 = mitochondrial_dna('''CCGGCAGGAACTGCTAACACATAATCGGCGTTTTGAAGTTCGCAAGGAGAGTGCCTTCCGCGTTCCCGGCTTCGTATCAA
GATATGGGACAAGGAGGGACTGTGACATTTACATGTTACGACCCCCAGTAGCTTAAATACAACGGTTTAACTACAATAGTCGGCCGGGCAGGGTTCGGATAGGTTTAC
AACTACTTTTCAACATTCTCGTGATTACTAAAGCAGTCAGGCAGAAGTGATCGAAGCGCTCTTAGATAGTGCGCCAGACCCGCTGAGCCCGCGGCATACTAGCAGTGA
ACCACGTTAGTACTGCTTGTATCGATTCGATACCCCTGAGGCCCGGACAATCTTCAGCCTTTCATAGAGAGTAAGTCTCATTTGAATTATAAACCTCGTTTATCCGTA
TCCCCTGGTCCTTTAGACCCCTGTTAAGTCTCCGGCTCGTTAGCTAGCTCGTAAGTTACTGTATTACTGGGTAGCGTTGTATAGATTTTTCGTGAGCCGTTTGGCTGC
CTAGCCACAGGAACGTCAAGGCGACGGTCCCTGATGATATGTAGAGTTGTGCTTCAGTGGGGGTCGTGGTTGACCAACAGGACCTCCTCATGCATTACACGGTGTAGTA
GAAACAGTTAAAAGGTTTATTCAGAAACCCCTTATGGTGGCTGCTTCGCAAGCCAGTGCCCATAGTTCACGTGAGCCTAATAATGGAATATAGGCCTGTGTGACTATGG
CCGGTTCCCTTCTCAAAGGAAGGTATCTACCATCGCCGTTACGTTCACCATAAGTTAACGTAAGGCATGCGTTCTACTCATTTTCCGCGAGGTAGCATAAGTTTTCTTA
ATCCACACATTGGACGCACTGGTGGGGTGAGCCAGAACCACGCATAGTTATAGGAGCGTGCGATTAAGCCTGCGACGACTCTTCGACCTATAACTCAAATTCATGCCAC
TGCTATCCAATCATTCTCCAATTGCGCAACATGGGATTTTACGCCGATTCTC''')

seq1.gc_content()
seq1.seq_len()
subseq = seq1.extract_seq()
seq1.gc_content(subseq)
seq1.seq_len(subseq)
motif1 = genomic_motif('ATG',subseq)
print(motif1.search_motif())
motif1.count_motif()

