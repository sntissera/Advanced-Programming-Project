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
            motif_results = seq_analyzed.search_motif()
            if isinstance(motif_results, list):
                conserved_positions[seq_id] = {
                    "Positions": ', '.join(map(str, motif_results)),
                    "Motif count": len(motif_results)
                }
            else:
                conserved_positions[seq_id] = {
                    "Positions": motif_results,
                    "Motif count": 0
                }
        return conserved_positions


class AlignmentAnalysis:
    def __init__(self, sequences: dict[str, str]) -> None:
        """Initialize with a dictionary of sequences."""
        self.sequences = sequences

    def pairwise_alignment(self, seq_id_1, seq_id_2, width: int = 80):
        """Perform pairwise alignment using Biopython and visually represent the alignment."""
        
        seq1 = str(self.sequences[seq_id_1])
        seq2 = str(self.sequences[seq_id_2])
        aligner = PairwiseAligner()

        score = aligner.score(seq1, seq2)

        
        # return {
        #     "score": score,
        #     "alignment vizualization": ,
        # }

        # aligner = PairwiseAligner()
        # aligner.mode = "global"

        # # Perform alignment
        # alignments = aligner.align(seq1, seq2)
        # alignment = alignments[0]

        # # Extract aligned sequences with gaps
        # aligned_seq1 = []
        # aligned_seq2 = []
        # for (start1, end1), (start2, end2) in zip(alignment.aligned[0], alignment.aligned[1]):
        #     aligned_seq1.append(seq1[start1:end1])
        #     aligned_seq2.append(seq2[start2:end2])

        # # Join fragments to create full aligned sequences
        # aligned_seq1 = ''.join(aligned_seq1)
        # aligned_seq2 = ''.join(aligned_seq2)

        # # Generate matches visualization including gaps
        # matches = []
        # for base1, base2 in zip(aligned_seq1, aligned_seq2):
        #     if base1 == "-" or base2 == "-":
        #         matches.append("-")  # Represent gaps with space
        #     else:
        #         matches.append("|" if base1 == base2 else "")
        # print(matches)
        # # Create the visualized matches string
        # matches_str = "".join(matches)

        # # Split into blocks of defined width
        # def split_str(sequence, width=80):
        #     return [sequence[i : i + width] for i in range(0, len(sequence), width)]
        
        # seq1_lines = split_str(aligned_seq1, width)
        # matches_lines = split_str(matches_str, width)
        # seq2_lines = split_str(aligned_seq2, width)
        
        # # Format the output
        # formatted = []
        # for line_1, line_match, line_2 in zip(seq1_lines, matches_lines, seq2_lines):
        #     formatted.append(line_1)
        #     formatted.append(line_match)
        #     formatted.append(line_2)
        #     formatted.append("")  # Add a blank line between blocks

        # # Return the alignment and score
        # return {
        #     "score": alignment.score,
        #     "alignment vizualization": "\n".join(formatted),
        # }

# Test data
sequences = {
    "seq1": "ACTGACTGACTG",
    "seq2": "ACTGACCGACTG",
}

seq_analysis = AlignmentAnalysis(sequences)
print(seq_analysis.pairwise_alignment("seq1","seq2"))