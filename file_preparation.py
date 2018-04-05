#!/usr/bin/env python
import os
import random
import logging
from argparse import ArgumentParser
import numpy as np

# Initialize seed for data reproducibility
random.seed(1) 

class FilePrep():
    """Constructor"""
    def __init__(self):
        self.reads = [] # [(header,sequence)]
        self.non_masked_reads = [] # [(header,sequence)] after removing hard and(or) soft masked reads
        self.seqs_to_remove = [] # headers of the sequences which will be removed at the end from the metagenome
        self.filtered_reads = [] # filtered lists containing (header,sequence), of remaining reads
        self.positive = [] # data for positive training classifications --> [(header,sequence)]
        self.negative = [] # data for negative training classifications --> [(header,sequence)]

    """Load positive and negative datasets into numpy arrays attributes"""
    def load_training_data(self,positive_dataset,negative_dataset):
        with open(positive_dataset,"U") as input_positive, open(negative_dataset,"U") as input_negative:
            for line in input_positive:
                tokens = line.split()
                header = tokens[0]
                sequence = tokens[1]
                self.positive.append((header,sequence))
            for line in input_negative:
                tokens = line.split()
                header = tokens[0]
                sequence = tokens[1]
                self.negative.append((header,sequence))

    """Convert ATCG from positive/negatives datasets to 0123 to one hot encoding"""
    def one_hot(self,tuple_to_convert):
        sequence_to_convert = tuple_to_convert[1] # index 0 = header, index 1 = sequence
        # Make sure that all nucleotides are upper case (may not be if user wished to keep soft masked reads)
        upper_sequence_to_convert = sequence_to_convert.upper()

        # Access Python 3's static str maketrans function to create translation from "ACGT" to "0123" numeric sequence
        numeric_string = upper_sequence_to_convert.translate(str.maketrans("ACGT","0123"))
        numeric_list = list(map(int,[digit for digit in numeric_string])) # List of 0,1,2,3 

        # Returns a one hot encoded the numpy numeric list with 4 class labels for A,C,G,T
        return keras.utils.np_utils.to_categorical(numeric_list,num_classes=4)

    def new_one_hot(self,tuple_to_convert):
        CHARS = 'ACGT'
        CHARS_COUNT = len(CHARS)
        sequence_to_convert = tuple_to_convert[1]
        upper_sequence_to_convert = sequence_to_convert.upper()
        max_len = len(sequence_to_convert)
        res = np.zeros((1, CHARS_COUNT * max_len),dtype=np.uint8) 
        seqlen = len(sequence_to_convert) 
        arr = np.chararray((seqlen,), buffer=seq)
        for index,char in enumerate(CHARS):
            res[index*seqlen:(ii+1)*seqlen][arr == char] = 1

    """Pre-processing: Load fasta file into reads property"""
    def load_fasta(self,fasta_file): 
        logging.debug("Loading fasta file: %s" % (fasta_file))
        # Reset reads being held in case using function to process many fasta files
        self.reads = []
        with open(fasta_file,'U') as input_file:
            tokens = input_file.readlines()
        tokens = map(str.strip,tokens) 
        cur_header, cur_seq = "",""
        index = 0
        for line in tokens:
            # First line
            if ">" in line and (index == 0):
                cur_header = line
                index += 1
            elif ">" in line:
                self.reads.append((cur_header,cur_seq))
                cur_header = line
                cur_seq = ""
            else:
                cur_seq += line
        # Add last sequence of file
        self.reads.append((cur_header,cur_seq))
        logging.debug("Finished loading fasta file: %s" % (fasta_file))

    """Pre-processing: Cut all sequences in each fasta file in in_dir into kmers of length 100 (or specified length) and write to a new file"""
    def write_kmers(self,negative_genomes_in_dir,output_file="merged_1.fasta",kmer_size=100):
        # Check for existing file
        if os.path.exists(output_file):
            # Grab current file number and increment
            new_file_num = int(output_file.split(".")[0].split("_")[1]) + 1
            output_file = "".join(["merged_",str(new_file_num),".fasta"]) 

        logging.debug("Writing kmers to output file: %s" % (output_file))
        with open(output_file,"w") as output_file:
            # Process all fasta sequences in specified directory
            for filename in os.listdir(negative_genomes_in_dir):
                if filename.endswith(".fasta") or filename.endswith(".fa") or filename.endswith(".fna"):
                    abs_path = os.path.join(in_dir,filename)
                    self.load_fasta(abs_path)
                    for cur_tuple in self.reads:
                        header = cur_tuple[0]
                        seq = cur_tuple[1]
                        # Loop over sequences to extract kmers with step size of kmer length
                        for index in range(0,len(seq)-kmer_size,kmer_size):
                            kmer = seq[index:index+kmer_size]
                            output_file.write("\t".join([header,kmer])+"\n") 
        logging.debug("Finished merging fasta files to: %s" % (output_file))

    """Randomly sample x number of DNA kmers of length k from a fasta file and write to a new file"""
    def sample_kmers(self,input_file,output_file="",kmer_size=100,num_seqs=200000,x=100000):
        counter = 0 # Will stop when 
        # Select random number of indices to include
        random.randint(0,len() - kmer_size)
         
    """Pre-processing: Remove any reads containing hard or soft masked regions from reads attribute by default"""
    def remove_masked_reads(self,soft=True):
        logging.debug("Removing masked reads...")
        reads_to_remove = set() # indices of reads to remove from reads attribute

        # Check if removing both hard and soft or only hard
        if soft == True:
            for index,cur_tuple in enumerate(self.reads):
                header = cur_tuple[0]
                seq = cur_tuple[1]

                # Remove read if it contains any hard masked (N's), soft masking (lowercase), or length < 100
                if("N" in seq or any(nucleotide for nucleotide in seq if nucleotide.islower()) or len(seq) < 100):
                    reads_to_remove.add(index)

        # Only remove hard masked reads
        else:
            for index,cur_tuple in enumerate(self.reads):
                header = cur_tuple[0]
                sequence = cur_tuple[1]
                if("N" in seq or len(seq) < 100):
                    reads_to_remove.add(index)

        # Store only non-masked reads in non_masked_reads attribute
        self.non_masked_reads = [cur_tuple for index,cur_tuple in enumerate(self.reads) if index not in reads_to_remove]
        logging.debug("Finished removing masked reads.")

    """Randomly samples x number of sequences from the training data file and outputs those sequences to a new file"""
    def random_sample_reads(self,infile,x=100000):
        # Write to a new file with added prefix "random"
        outfile = "random_" + infile
        # Make sure x is smaller than the number of sequences in the file
        with open(infile,"U") as input_file:
            tokens = input_file.readlines()
        if x > len(tokens):
            raise ValueError("Number of specified sequences is larger than number of lines in file.")
        random_indices = random.sample(xrange(len(tokens)),x)
        # Selects randomly selected sequences
        with open(outfile,"w") as output_file:
            for index in random_indices:
                output_file.write(tokens[index])

    """Post-processing: Use classification labels to remove sequences classified as human by deep learning."""
    def filter_reads(self,classification_labels):    
        logging.debug("Removing reads classified as human DNA...")
        for cur_tuple,classification in zip(self.read_list,classification_labels):
            # Sequence labeled as non-human, keep in filtered list
            if classification == 0:
                self.filtered_reads.append(cur_tuple)
        logging.debug("Finished removing reads classified as human DNA.")

    """Post-processing: Writes new fasta file based on filtered_list"""
    def write_fasta(self,output_file):
        filtered_list = self.filtered_reads
        logging.debug("Writing final filtered fasta file to: %s..." % (output_file))
        with open(output_file,'w') as output_file:
            for cur_tuple in filtered_list:
                header = cur_tuple[0]
                sequence = cur_tuple[1]
                output_file.write("\n".join([header,sequence]))
        logging.debug("Finished writing final filtered fasta file to: %s" % (output_file))

def main():
    parser = ArgumentParser(description="Meta-Sweeper is a program to filter human DNA from metagenomic reads.")
    parser.add_argument('-i','--input_fasta',help="Specify input fasta file name to be filtered",
                        required=True,type=str)
    parser.add_argument('-o','--output_fasta',help="Specify final output fasta file name with filtered reads",
                        required=True,type=str)
    parser.add_argument('-d','--debug',help="Enable debug logging if flag is specified.",
                        required=False,action="store_true",dest="bool_switch_option",default=False)
    args = parser.parse_args()   

    # Configure log file destination if running in debug mode
    if args.bool_switch_option == True:
        logging.basicConfig(filename="Meta_Sweeper.log",
                            filemode='w',
                            format='%(asctime)s %(levelname)s %(message)s',
                            level=logging.DEBUG)

    # Build neural network models
    NN = NeuralNetwork()
#    NN.load_fasta(args.input_fasta)
  #  NN.load_training_data("random_positive_training.tsv","random_negative_training.tsv")  
    # Randomly sample x number of reads for positive and negative training
#    NN.random_sample_reads("positive_training.tsv")
 #   NN.random_sample_reads("negative_training.tsv")

    # Write source file containing 100bp long reads in merged_1.fasta file 
#    NN.write_kmers("C:\\Users\\Kellie\\Desktop\\EC2\kellie\\Machine_Learning\\negative_genomes",output_file="merged_1.fasta",kmer_size=100)

if __name__ == "__main__":
    main()