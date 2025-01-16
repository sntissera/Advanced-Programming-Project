from Bio.Seq import Seq
from Bio.SeqUtils import nt_search
from Bio.Align import PairwiseAligner

from mtDNA_parser import MitochondrialDNAParser
from mtDNA import MitochondrialDna, GenomicMotif

class ComparativeAnalysis:
    def __init__(self, sequences):
        """
        Initializes with a dictionary of all the sequences to be analyzed;
        keys are the sequences' IDs.
        The sequences are passed as instances of MitochondrialDna.
        """
        self.sequences = sequences

    def summary(self):
        """
         Summarizes length and GC content of all sequences.
         """
        summaries ={}
        for seq_id, seq in self.sequences.items():
            sequence = MitochondrialDna(seq)
            seq_gc = sequence.gc_content()
            seq_length = sequence.seq_len()
            summaries[seq_id] = {
                "Length": seq_length,
                "GC Content": seq_gc,
            }
        return summaries
        
    def summarize_findings(self):
        """
        Prints summary findings such as longest/shortest sequence and highest GC content.
        """
        summaries = self.summary()
        longest = max(summaries.items(), key=lambda x: x[1]['Length'])
        shortest = min(summaries.items(), key=lambda x: x[1]['Length'])
        highest_gc = max(summaries.items(), key=lambda x: x[1]['GC Content'])
        lowest_gc = min(summaries.items(), key=lambda x: x[1]['GC Content'])

        print(f"The longest sequence is {longest[0]} with {longest[1]['Length']} bases.")
        print(f"The shortest sequence is {shortest[0]} with {shortest[1]['Length']} bases.")
        print(f"The sequence with the highest GC content is {highest_gc[0]} with {highest_gc[1]['GC Content']:.2f}% GC content.")
        print(f"The sequence with the lowest GC content is {lowest_gc[0]} with {highest_gc[1]['GC Content']:.2f}% GC content.")

class ConservedMotifs:
    def __init__(self, sequences):
        """
        Initializes with a dictionary of sequences where keys are sequence IDs.
        The sequences are passed as instances of MitochondrialDna.
        """
        self.sequences = sequences

    def conserved_motifs(self, motif):
        """
        Identifies conserved motifs across multiple sequences.
        Takes as input a motif and gives as output a dictionary where:
        - keys are IDs of sequences in which the motif is found;
        - values are lists of positions in which the motif is found in that sequence
        """
        conserved_positions = {}
        for seq_id, seq in self.sequences.items():
            seq_analyzed = GenomicMotif(motif, seq)
            positions = seq_analyzed.search_motif()
            conserved_positions[seq_id] = {
                "Positions": positions,
                "Motif count": len(positions)
            }     
        return conserved_positions


class AlignmentAnalysis:

    def __init__(self, sequences):
        self.sequences = sequences

    def pairwise_alignment(self, seq_id_1, seq_id_2):
        """Perform pairwise alignment using Biopython and visually represent the alignment."""
        if seq_id_1 not in self.sequences or seq_id_2 not in self.sequences:
            raise ValueError("One or both sequence IDs not found.")

        seq1 = str(self.sequences[seq_id_1]._sequence)
        seq2 = str(self.sequences[seq_id_2]._sequence)

        # create the aligner
        aligner = PairwiseAligner()
        aligner.mode = 'global'
        alignments = aligner.align(seq1, seq2)

        # get the first alignment (best alignment)
        alignment = alignments[0]
        aligned_seq1 = alignment.aligned[0]
        aligned_seq2 = alignment.aligned[1]

        # Generate the visualization
        visual_seq1, visual_seq2, matches = [], [], []
        for (start1, end1), (start2, end2) in zip(aligned_seq1, aligned_seq2):
            for i in range(start1, end1):
                visual_seq1.append(seq1[i])
            for i in range(start2, end2):
                visual_seq2.append(seq2[i])

        for base1, base2 in zip(visual_seq1, visual_seq2):
            if base1 == base2:
                matches.append("|")
            else:
                matches.append(" ")

        # convert lists to strings for printing
        aligned_seq1_str = "".join(visual_seq1)
        matches_str = "".join(matches)
        aligned_seq2_str = "".join(visual_seq2)

        # print the alignment
        print(aligned_seq1_str)
        print(matches_str)
        print(aligned_seq2_str)

        return aligned_seq1_str, matches_str, aligned_seq2_str

