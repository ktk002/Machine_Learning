#!/usr/bin/env python
import os
import random

# Set seed to random arbitrary value for reproducibility
random.seed(100)

MITOCHONDRIA_DATA_PATH = "C:\\Users\\Kellie\\Desktop\\Mitochondrial_DNA"
MITOCHONDRIA_FILE_NAME = "NC_012920.1.fasta"
INTRON_DATA_PATH = "C:\\Users\\Kellie\\Desktop\\Machine_Learning\\Intron_DNA"
SINE_DATA_PATH = "C:\\Users\\Kellie\\Desktop\\Machine_Learning\\SINE_DNA"
ALU_DATA_PATH = "C:\\Users\\Kellie\\Desktop\\Machine_Learning\\Alus"
EXON_DATA_PATH = "C:\\Users\\Kellie\\Desktop\\Machine_Learning\\Exon_DNA"
MITOCHONDRIAL_DATA_PATH = "C:\\Users\\Kellie\\Desktop\\Machine_Learning\\Mitochondrial_DNA"

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

    def get_random_indices(self, total_num_sequences, num_to_sample):
        """Returns a list of num_to_sample randomly selected indices from total_num_sequences"""
        random_sampled_list = list(random.sample(range(0, total_num_sequences), num_to_sample))
        return random_sampled_list

    def get_sampled_indices(self, list_to_sample, random_indices):
        """Returns list of sampled items from input list."""
        sampled_elements = [list_to_sample[x] for x in random_indices]
        return sampled_elements


def main():
    fasta_processor = ProcessFasta()
    # 1) # Create negative training data
    #  Merge all negative genomes into 1 fasta file
    # fasta_processor.combine_fastas("C:\\Users\\Kellie\\Desktop\\EC2\kellie\\Machine_Learning\\negative_genomes",output_file="merged_1.fasta")
    # fasta_processor.write_fasta_dict(output_fasta="all_negative.fasta")
    # fasta_processor.write_training_data(window_size=100, sliding_window=False, output_file="all_negative.tsv")
    #
    # 2) Load first set of Alu sequences from BLAST Alu database
    # fasta_processor.load_multifasta("C:\\Users\\Kellie\\Desktop\\Machine_Learning\\Alus\\alu.n")
    # fasta_processor.write_training_data(window_size=100, output_file="alu_sliding_1.tsv")
    #
    # 3) Create second set of alu sequences tsv from NIH Alu database
    # fasta_processor.load_multifasta("C:\\Users\\Kellie\\Desktop\\Machine_Learning\\Alus\\Alu_seqs.fasta")
    # fasta_processor.write_training_data(window_size=100, sliding_window=False, output_file="alu_cut_2.tsv")

    #4) Create positive training set for intron sequences by randomly picking 500 sequences and writing to: "500_positive_test_introns.tsv"
    # Write remaining sequences to: "training_intron_sequences.tsv"
    # ALL_INTRON_PATH = os.path.join(INTRON_DATA_PATH, "all_intron_sequences.tsv")
    # TEST_INTRONS = os.path.join(INTRON_DATA_PATH, "500_positive_test_introns.tsv")
    # TRAINING_INTRONS = os.path.join(INTRON_DATA_PATH, "training_intron_sequences.tsv")
    # with open(ALL_INTRON_PATH, "r") as file_reader:
    #     all_intron_list = file_reader.readlines()
    # total_num_seq = len(all_intron_list)
    # num_sequences_to_sample = 500
    # all_indices = [index for index in range(total_num_seq)]
    # random_indices = fasta_processor.get_random_indices(total_num_seq, num_sequences_to_sample)
    # random_sequences = fasta_processor.get_sampled_indices(all_intron_list, random_indices)
    # # Write new results to tsv: "500_positive_test_introns.tsv"
    # with open(TEST_INTRONS, "w") as file_writer:
    #     for cur_line in random_sequences:
    #         file_writer.write(cur_line)
    #
    # # Write remaining results not found in 500_positive_test_introns.tsv to "training_intron_sequences"
    # # Take difference in remaining indices
    # training_indices = list(set(all_indices) - set(random_indices))
    # training_sequences = fasta_processor.get_sampled_indices(all_intron_list, training_indices)
    # with open(TRAINING_INTRONS, "w") as file_writer:
    #     for cur_line in training_sequences:
    #         file_writer.write(cur_line)

    #5) Create negative testing dataset
    # fasta_processor.combine_fastas("C:\\Users\\Kellie\\Desktop\\Negative_genomes\\Level_5_genomes")
    # fasta_processor.write_fasta_dict(output_file="testing_genomes_set.fasta")
    # fasta_processor.write_training_data(window_size=100, output_file="5_negative_test_genomes.tsv")

    #6) Select n number of random samples from negative training set (all_negative.tsv) and write to new tsv
    # Open tsv to sample from
    # with open("all_negative.tsv", "r") as file_reader:
    #     all_negative_list = file_reader.readlines()
    # total_num_seq = len(all_negative_list)
    # num_sequences_to_sample = 100
    # random_indices = fasta_processor.get_random_indices(total_num_seq, num_sequences_to_sample) # Returns list of strings containing randomly chosen sequences
    # random_sequences = fasta_processor.get_sampled_indices(all_negative_list, random_indices)
    # print("Length of random_sequences is: ", str(len(random_sequences)))
    # print("First line is: ", str(random_sequences[0]))
    # # Write sampled results to a new tsv
    # with open("subsampled_negative_sequences.tsv", "w") as file_writer:
    #     for cur_line in random_sequences:
    #         file_writer.write(cur_line)
    #
    # Convert exon fasta to tsv with all exon sequences
    # fasta_processor.load_multifasta("C:\\Users\\Kellie\\Desktop\\Machine_Learning\\exons.fasta")
    # fasta_processor.write_training_data(window_size=100, output_file="all_exon_sequences.tsv")

    # 7) Create positive (human DNA) exon training dataset of 100 samples
    # Select n number of random samples from exon positive training data and write to new tsv
    # with open("all_exon_sequences.tsv", "r") as file_reader:
    #     all_exon_list = file_reader.readlines()
    # total_num_seq = len(all_exon_list)
    # num_sequences_to_sample = 100
    # random_indices = fasta_processor.get_random_indices(total_num_seq, num_sequences_to_sample) # Returns list of strings containing randomly chosen sequences
    # random_sequences = fasta_processor.get_sampled_indices(all_exon_list, random_indices)
    # # Write new results to tsv
    # with open("subsampled_exons_100.tsv", "w") as file_writer:
    #     for cur_line in random_sequences:
    #         file_writer.write(cur_line)

    # 8) Create mitochondria positive training dataset
    # MITOCHONDRIA_FASTA_PATH = os.path.join(MITOCHONDRIA_DATA_PATH, MITOCHONDRIA_FILE_NAME)
    # # Convert fasta to tsv
    # fasta_processor.load_multifasta(MITOCHONDRIA_FASTA_PATH)
    # OUTPUT_MITOCHONDRIA_PATH = os.path.join(MITOCHONDRIA_DATA_PATH, "all_mitochondrial_sequences.tsv")
    # fasta_processor.write_training_data(window_size=100, output_file=OUTPUT_MITOCHONDRIA_PATH)

    # 9) Randomly select 500 mitochondrial sequences and write them to a file called "500_positive_test_mitochondria.tsv"
    # Write the remaining sequences to a file called "training_mitochondrial_sequences.tsv"
    # ALL_MITOCHONDRIA_PATH = os.path.join(MITOCHONDRIA_DATA_PATH, "all_mitochondrial_sequences.tsv")
    # TEST_MITOCHONDRIA = os.path.join(MITOCHONDRIA_DATA_PATH, "500_positive_test_mitochondria.tsv")
    # TRAINING_MITOCHONDRIA = os.path.join(MITOCHONDRIA_DATA_PATH, "training_mitochondrial_sequences.tsv")
    # with open(ALL_MITOCHONDRIA_PATH, "r") as file_reader:
    #     all_mitochondria_list = file_reader.readlines()
    # total_num_seq = len(all_mitochondria_list)
    # num_sequences_to_sample = 500
    # all_indices = [index for index in range(total_num_seq)]
    # random_indices = fasta_processor.get_random_indices(total_num_seq, num_sequences_to_sample)
    # random_sequences = fasta_processor.get_sampled_indices(all_mitochondria_list, random_indices)
    # # Write new results to tsv: "500_positive_test_mitochondria.tsv"
    # with open(TEST_MITOCHONDRIA, "w") as file_writer:
    #     for cur_line in random_sequences:
    #         file_writer.write(cur_line)
    #
    # # Write remaining results not found in 500_positive_test_mitochondria.tsv to "training_mitochondrial_sequences"
    # # Take difference in remaining indices
    # training_indices = list(set(all_indices) - set(random_indices))
    # training_sequences = fasta_processor.get_sampled_indices(all_mitochondria_list, training_indices)
    # with open(TRAINING_MITOCHONDRIA, "w") as file_writer:
    #     for cur_line in training_sequences:
    #         file_writer.write(cur_line)

    # 10) Convert SINE fasta file into 100 bp sequences in tsv format
    # Convert fasta to tsv
    # fasta_processor.load_multifasta("C:\\Users\\Kellie\\Desktop\\Machine_Learning\\SINE_DNA\\SINEs.bnk")
    # fasta_processor.write_training_data(window_size=100, sliding_window=True, output_file="C:\\Users\\Kellie\\Desktop\\Machine_Learning\\SINE_DNA\\all_sines_sliding_window.tsv")

    # 11) Randomly select 500 SINE sequences to be used as testing data
    # Write the remaining sequences to a file called "training_sine_sequences.tsv"
    # ALL_SINE_PATH = os.path.join(SINE_DATA_PATH, "all_sines_sliding_window.tsv")
    # TEST_SINE = os.path.join(SINE_DATA_PATH, "500_positive_test_sines.tsv")
    # TRAINING_SINE = os.path.join(SINE_DATA_PATH, "training_sine_sequences.tsv")
    # with open(ALL_SINE_PATH, "r") as file_reader:
    #     all_sine_list = file_reader.readlines()
    # total_num_seq = len(all_sine_list)
    # num_sequences_to_sample = 500
    # all_indices = [index for index in range(total_num_seq)]
    # random_indices = fasta_processor.get_random_indices(total_num_seq, num_sequences_to_sample)
    # random_sequences = fasta_processor.get_sampled_indices(all_sine_list, random_indices)
    # # Write new results to tsv: "500_positive_test_sine.tsv"
    # with open(TEST_SINE, "w") as file_writer:
    #     for cur_line in random_sequences:
    #         file_writer.write(cur_line)
    #
    # # Write remaining results not found in 500_positive_test_sine.tsv to "training_sine_sequences"
    # # Take difference in remaining indices
    # training_indices = list(set(all_indices) - set(random_indices))
    # training_sequences = fasta_processor.get_sampled_indices(all_sine_list, training_indices)
    # with open(TRAINING_SINE, "w") as file_writer:
    #     for cur_line in training_sequences:
    #         file_writer.write(cur_line)
    #
    # 12) Concatenate all positive training data into a single file called "multiclass_training_data.tsv" with added column
    # containing single word to represent data type (ex row: ">alu_sequence1 ACT alu")
    # intron_path = os.path.join(INTRON_DATA_PATH, "all_intron_sequences.tsv")
    # exon_path = os.path.join(EXON_DATA_PATH, "all_exon_sequences.tsv")
    # alu_path = os.path.join(ALU_DATA_PATH, "all_alu.tsv")
    # sine_path = os.path.join(SINE_DATA_PATH, "all_sines.tsv")
    # mitochondria_path = os.path.join(MITOCHONDRIAL_DATA_PATH, "all_mitochondrial_sequences.tsv")
    # multiclass_training_path = "multiclass_merged_training_data.tsv"
    # data_type_list = ['intron', 'exon', 'alu', 'sine', 'mitochondria']
    # with open(multiclass_training_path, 'w') as outfile:
    #     for training_type in data_type_list:
    #         with open(intron_path, 'r') as intron_infile, \
    #                 open(exon_path, 'r') as exon_infile, \
    #                 open(alu_path, 'r') as alu_infile, \
    #                 open(sine_path, 'r') as sine_infile, \
    #                 open(mitochondria_path, 'r') as mitochondria_infile:
    #             intron_lines = intron_infile.readlines()
    #             exon_lines = exon_infile.readlines()
    #             alu_lines = alu_infile.readlines()
    #             sine_lines = sine_infile.readlines()
    #             mitochondrial_lines = mitochondria_infile.readlines()
    #
    #             for line in intron_lines:
    #                 tokens = line.strip().split()
    #                 tokens.append('intron')
    #                 new_line = "\t".join(tokens) + "\n"
    #                 outfile.write(new_line)
    #
    #             for line in exon_lines:
    #                 tokens = line.strip().split()
    #                 tokens.append('exon')
    #                 new_line = "\t".join(tokens) + "\n"
    #                 outfile.write(new_line)
    #
    #             for line in alu_lines:
    #                 tokens = line.strip().split()
    #                 tokens.append('alu')
    #                 new_line = "\t".join(tokens) + "\n"
    #                 outfile.write(new_line)
    #
    #             for line in sine_lines:
    #                 tokens = line.strip().split()
    #                 tokens.append('sine')
    #                 new_line = "\t".join(tokens) + "\n"
    #                 outfile.write(new_line)
    #
    #             for line in mitochondrial_lines:
    #                 tokens = line.strip().split()
    #                 tokens.append('mitochondria')
    #                 new_line = "\t".join(tokens) + "\n"
    #                 outfile.write(new_line)


if __name__ == '__main__':
    main()
