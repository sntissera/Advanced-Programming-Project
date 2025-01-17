import io
from Bio.Seq import Seq
from Bio.SeqUtils import nt_search
import matplotlib.pyplot as plt 

class MitochondrialDna: 
    def __init__ (self, sequence): 
        self._sequence = Seq(sequence)
    
    def extract_seq(self, start, end):
        '''Extracts subsequences from genomic sequences'''
        if start < 1 or end > len(self._sequence) or start > end:
            return ("Invalid start or end position.")
        return str(self._sequence[start-1:end])
        

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
        self.motif = motif.upper()
       
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

    def visualize_motif(self):
        '''Visualize the location of a motif in the genome'''
        
        sequence_str = str(self._sequence)
        positions = self.search_motif()
        positions = [pos-1 for pos in positions]

        #empty plot
        plt.figure(figsize=(10,2))
        plt.plot([0,len(sequence_str)],[1,1], color = 'gray', lw=2, label = 'Genome' )

        #vertical lines
        for pos in positions:
            plt.axvline(pos,color='red', linestyle = '--',lw=1,label='Motif Location' if pos == positions[0] else "" )
            plt.text(pos,1.05,str(pos+1),color = 'blue', fontsize = 8, ha='center',rotation=90)
            
        #plot
        plt.title('Motif Location in Genome')
        plt.xlabel('Position in Genome')
        plt.ylabel('Presence of Motif')
        plt.yticks([])
        plt.legend(loc='upper right')
        plt.tight_layout()

        img = io.BytesIO()
        plt.savefig(img, format='png')
        img.seek(0)
        plt.close()
        return img

