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
    
        summary_comparison = {}
        summary_comparison["longest"] = max(summaries.items(), key=lambda x: x[1]['Length'])
        summary_comparison["shortest"] = min(summaries.items(), key=lambda x: x[1]['Length'])
        summary_comparison["highest_gc"] = max(summaries.items(), key=lambda x: x[1]['GC Content'])
        summary_comparison["lowest_gc"] = min(summaries.items(), key=lambda x: x[1]['GC Content'])

        return summary_comparison

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
            seq_analyzed = GenomicMotif(seq, motif)
            positions = seq_analyzed.search_motif()
            conserved_positions[seq_id] = {
                "Positions": ', '.join(map(str, positions)),
                "Motif count": len(positions)
            }     
        return conserved_positions


class AlignmentAnalysis:
    def __init__(self, sequences):
        self.sequences = sequences

    def pairwise_alignment(self, seq_id_1, seq_id_2):
        """Perform pairwise alignment using Biopython and visually represent the alignment."""
        if seq_id_1 not in self.sequences or seq_id_2 not in self.sequences:
            return {
                "Error": "One or both sequence IDs not found. Please provide valid sequence IDs."
                }
            
        seq1 = str(self.sequences[seq_id_1])
        seq2 = str(self.sequences[seq_id_2])

        # create the aligner
        aligner = PairwiseAligner()
        aligner.mode = "global"
        
        # perform alignment
        alignments = aligner.align(seq1, seq2)
        alignment = alignments[0]

        # extract aligned sequences
        aligned_seq1, aligned_seq2 = alignment.aligned

        # generate the visualization
        visual_seq1, visual_seq2, matches = [], [], []
        for (start1, end1), (start2, end2) in zip(aligned_seq1, aligned_seq2):
            for i in range(start1, end1):
                visual_seq1.append(seq1[i])
            for i in range(start2, end2):
                visual_seq2.append(seq2[i])

        for base1, base2 in zip(visual_seq1, visual_seq2):
            matches.append("|" if base1 == base2 else " ")

        # convert lists to strings
        aligned_seq1_str = "".join(visual_seq1)
        matches_str = "".join(matches)
        aligned_seq2_str = "".join(visual_seq2)

        # split into blocks of defined width
        def split_str(sequence, width=80):
            return [sequence[i : i + width] for i in range(0, len(sequence), width)]
        
        seq1_lines = split_str(aligned_seq1_str, line_width)
        matches_lines = split_str(matches_str, line_width)
        seq2_lines = split_str(aligned_seq2_str, line_width)
        
        formatted = []
        for line_1, line_match, line_2 in zip(seq1_lines, matches_lines, seq2_lines):
            formatted.append(line_1)
            formatted.append(line_match)
            formatted.append(line_2)
            formatted.append("")
            
        # return the alignment
        return "\n".join(formatted)
