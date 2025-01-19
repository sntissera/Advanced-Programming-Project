import pandas as pd

class MitochondrialDNAParser:

    def __init__(self, file_path):
        self.file_path = file_path
        self.data = None
        self._parse_fasta()

    def _parse_fasta(self): 
        sequence_ids = []
        descriptions = []
        sequences = []

        with open(self.file_path, 'r') as file:
            current_sequence = []
            sequence_id = ""
            description = ""

            for line in file:
                line = line.strip()
                if line.startswith(">"):
                    if current_sequence:  
                        sequences.append("".join(current_sequence))
                        current_sequence = []

                    header_line = line[1:].strip().split(maxsplit=1)
                    sequence_id = header_line[0]
                    description = header_line[1] if len(header_line) > 1 else ""
                    sequence_ids.append(sequence_id)
                    descriptions.append(description)
                else:
                    current_sequence.append(line)

            if current_sequence:
                sequences.append("".join(current_sequence))

        self.data = pd.DataFrame({
            "Sequence ID": sequence_ids,
            "Description": descriptions,
            "Sequence": sequences
        })

    def get_data(self):
        return self.data

    def get_sequence_by_id(self, sequence_id):
        result = self.data[self.data["Sequence ID"] == sequence_id]
        if not result.empty:
            return result["Sequence"].values[0]
        return None

    def get_all_sequences(self):
        return self.data["Sequence ID"].tolist()

