#!/usr/bin/env python

class ProcessFasta(object):
    """Perform sliding window of size x and """
    def __init__(self):
        self.fasta_sequences = {}

    def write_training_data(self, window_size, keep_N=False, sliding_window=True, output_file="sliding_window_training.tsv"):
        """Writes fasta sequences to tsv using sliding window and does not keep masked sequences by default."""
        with open(output_file, 'w') as output_handler:
            if keep_N == True:
                #for cur_sequence_identifier, cur_fasta_sequence in self.fasta_sequences.items():
                pass 
            # Don't write reads containing "N" to tsv file
            else:
                for cur_sequence_identifier, cur_fasta_sequence in self.fasta_sequences.items():
                    if sliding_window == True:
                        # Sliding window, write to tsv if no N
                        for index in range(0, len(cur_fasta_sequence) - window_size + 1):
                            cur_spliced_sequence = cur_fasta_sequence[index:index + window_size]
                            if "N" not in cur_spliced_sequence and len(cur_spliced_sequence) == window_size:
                                output_handler.write("\t".join([cur_sequence_identifier, cur_spliced_sequence]) + "\n")
                    else:
                        # Cut fasta sequences into fragments with step size of window_size if multifasta is too large
                        for index in range(0, len(cur_fasta_sequence) - window_size + 1, window_size):
                            cur_spliced_sequence = cur_fasta_sequence[index:index + window_size]
                            if "N" not in cur_spliced_sequence and len(cur_spliced_sequence) == window_size:
                                output_handler.write("\t".join([cur_sequence_identifier, cur_spliced_sequence]) + "\n")

    def load_multifasta(self, fasta_file):
        self.fasta_sequences = {}
        with open(fasta_file, 'r') as input_file:
            for index, line in enumerate(input_file):
                if index == 0:
                    cur_key = line.strip()
                    cur_sequence = ""
                elif ">" in line:
                    self.fasta_sequences[cur_key] = cur_sequence
                    cur_key = line.strip()
                    cur_sequence = ""
                else:
                    cur_sequence = cur_sequence + line.upper().strip()
            # Add final sequence to dictionary
            self.fasta_sequences[cur_key] = cur_sequence

def main():
    fasta_processor = ProcessFasta()

    # Load first set of Alu sequences from BLAST Alu database
    fasta_processor.load_multifasta("C:\\Users\\Kellie\\Desktop\\Machine_Learning\\Alus\\alu.n")
    fasta_processor.write_training_data(window_size=100, output_file="alu_sliding_1.tsv")

    # Create second set of alu sequences tsv from NIH Alu database
    fasta_processor.load_multifasta("C:\\Users\\Kellie\\Desktop\\Machine_Learning\\Alus\\Alu_seqs.fasta")
    fasta_processor.write_training_data(window_size=100, sliding_window=False, output_file="alu_cut_2.tsv")

if __name__ == '__main__':
    main()
