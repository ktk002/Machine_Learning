#!/usr/bin/env python
import os

class ProcessFasta(object):
    """Perform sliding window of size x and generate tsv with sequence id and sequence of length x."""
    def __init__(self):
        self.fasta_sequences = {}

    def write_training_data(self, window_size=100, keep_N=False, sliding_window=True, output_file="sliding_window_training.tsv"):
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
                            # Only keep fragments which contain characters "ATCG"
                            if not any(i not in 'ACTG' for i in cur_spliced_sequence) and len(cur_spliced_sequence) == window_size:
                                output_handler.write("\t".join([cur_sequence_identifier, cur_spliced_sequence]) + "\n")
                    else:
                        # Cut fasta sequences into fragments with step size of window_size if multifasta is too large
                        for index in range(0, len(cur_fasta_sequence) - window_size + 1, window_size):
                            cur_spliced_sequence = cur_fasta_sequence[index:index + window_size]
                            if not any(i not in 'ACTG' for i in cur_spliced_sequence) and len(cur_spliced_sequence) == window_size:
                                output_handler.write("\t".join([cur_sequence_identifier, cur_spliced_sequence]) + "\n")

    def load_multifasta(self, fasta_file, clear_fasta_dict=True):
        if clear_fasta_dict == True:
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

    def combine_fastas(self, in_dir):
        """Pre-processing: Merge all fasta files into fasta sequences dict"""
        for filename in os.listdir(in_dir):
            if filename.endswith(".fasta") or filename.endswith(".fa") or filename.endswith(".fna"):
                print("name of fasta is: ", filename)
                abs_path = os.path.join(in_dir, filename)
                # Load fasta sequence into self.fasta_sequences attribute
                self.load_multifasta(abs_path, clear_fasta_dict=False)

    def write_fasta_dict(self, output_file="merged_1.fasta"):
        file_count = 1

        # Check for existing file, increment counter until file doesn't exist
        while os.path.exists(output_file):
            # Grab current file number and increment
            tokens = output_file.split('.')
            base_name = tokens[0]
            extension = tokens[1]
            output_file = base_name + '_' + str(file_count) + '.' + extension

        print("Writing to new fasta file: ", output_file)
        with open(output_file, "w") as output_fasta:
            # Process all fasta sequences in specified directory
            for sequence_identifier, sequence in self.fasta_sequences.items():
                output_fasta.write(sequence_identifier + "\n")
                output_fasta.write(sequence + "\n")


def main():
    fasta_processor = ProcessFasta()

    # Create negative training data
    # Merge all negative genomes into 1 fasta file
#    fasta_processor.combine_fastas("C:\\Users\\Kellie\\Desktop\\EC2\kellie\\Machine_Learning\\negative_genomes",output_file="merged_1.fasta")
 #   fasta_processor.write_fasta_dict(output_fasta="all_negative.fasta")
 #   fasta_processor.write_training_data(window_size=100, sliding_window=False, output_file="all_negative.tsv")

    # Load first set of Alu sequences from BLAST Alu database
 #   fasta_processor.load_multifasta("C:\\Users\\Kellie\\Desktop\\Machine_Learning\\Alus\\alu.n")
 #   fasta_processor.write_training_data(window_size=100, output_file="alu_sliding_1.tsv")

    # Create second set of alu sequences tsv from NIH Alu database
 #   fasta_processor.load_multifasta("C:\\Users\\Kellie\\Desktop\\Machine_Learning\\Alus\\Alu_seqs.fasta")
 #   fasta_processor.write_training_data(window_size=100, sliding_window=False, output_file="alu_cut_2.tsv")

    # Create positive training set for intron sequences

    # Create negative testing dataset
    # fasta_processor.combine_fastas("C:\\Users\\Kellie\\Desktop\\Negative_genomes\\Level_5_genomes")
    # fasta_processor.write_fasta_dict(output_file="testing_genomes_set.fasta")
    # fasta_processor.write_training_data(window_size=100, output_file="5_negative_test_genomes.tsv")


if __name__ == '__main__':
    main()
