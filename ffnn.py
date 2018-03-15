#!/usr/bin/env python
import logging
import keras.utils
import numpy as np
from argparse import ArgumentParser
from keras.models import Sequential
from keras.layers import Dense, Activation
#from string import maketrans

class NeuralNetwork():
    """Constructor"""
    def __init__(self):
        self.reads = [] # [(header,sequence)]
	self.non_masked_reads = [] # [(header,sequence)] after removing hard and(or) soft masked reads
        self.seqs_to_remove = [] # headers of the sequences which will be removed at the end from the metagenome
        self.filtered_reads = [] # filtered lists containing (header,sequence), of remaining reads
	self.one_hot_matrix = None # one hot encoded data for training

    """Create feed forward neural network models""" 
    def build_FFNN(self,model_name):
        log.debug("Building feedforward neural network: %s" % (model_name))
        # 1. Rectified linear activation function model
        model1 = Sequential()
        # Configure first layer to have output arrays of shape (*,32) 
        # Configure first layer to take input arrays of shape (*,16),
        # Note: no need to specify dimension of input in additional layers  
        model1.add(Dense(32,input_dim=784))
        # Rectified linear activation function
        model1.add(Activation('relu'))
        # Compile model to configure learning process
        model1.compile(optimizer='rmsprop',
                        loss='binary_crossentropy',
                        metrics=['accuracy'])

        # 2. Sigmoid (logistic) activation function model
        model2 = Sequential()
        model2.add(Dense(32,input_dim=784))
        model2.add(Activation('sigmoid'))
        model2.compile()

        # 3. Tanh activation function model
        model3 = Sequential()
        model3.add(Dense(32,input_dim=784,activation='tanh'))
        model3.compile()

        # 4. Linear activation function model
        model4 = Sequential()
        model4.add(Dense(,activation='linear'))
        model4.compile()

        # 5. Scaled Exponential Linear Unit model(Klambauer et al., 2017)
        model5 = Sequential()
        model5.add(Dense(,activation='selu'))
        model5.compile()

        # Creates 10 row vectors, each with 5 features  
        #data = np.random.random((10,1))
        my_list = ["ACT","CTA"]
        data = np.array(my_list)
        print(data)
        # Makes 10 row vectors, each a dimension of 1, with a value in between 0 to 2, not including 2
        labels = np.random.randint(4, size=(10, 1)) 
        print(labels)
        one_hot_labels = keras.utils.to_categorical(labels,num_classes=2)
        #print(one_hot_labels)
        logging.debug("Finished building feedforward neural network: %s" % (model_name))

    """Convert ATCG to 1234 to one hot encoding"""
    def one_hot(self,sequence_to_convert):
	# Make sure that all nucleotides are upper case (may not be if user wished to keep soft masked reads)
	upper_sequence_to_convert = sequence_to_convert.upper()

        # Access Python 3's static str maketrans function to create translation from "ACTG" to "1234" numeric sequence
        numeric_string = upper_sequence_to_convert.translate(str.maketrans("ACTG","1234"))
        numeric_list = list(map(int,[digit for digit in numeric_string]))
	
	# One hot encode the numeric list
	keras.utils.to_categorical(numeric_list,num_classes=2)

    """Pre-processing: Load fasta file into reads property"""
    def load_fasta(self,fasta_file): 
        logging.debug("Loading fasta file: %s" % (fasta_file))
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
                self.read_dict.append((cur_header,cur_seq))
                cur_header = line
                cur_seq = ""
            else:
                cur_seq += line
        # Add last sequence of file
        self.read_dict.append((cur_header,cur_seq))
        logging.debug("Finished loading fasta file: %s" % (fasta_file))

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
				reads_to_remove.append(index)

	# Only remove hard masked reads
	else:
		for index,cur_tuple in enumerate(self.reads):
			header = cur_tuple[0]
			sequence = cur_tuple[1]
			if("N" in seq or len(seq) < 100):
				reads_to_remove.append(index)

	# Store only non-masked reads in non_masked_reads attribute			
	self.non_masked_reads = [cur_tuple for index,cur_tuple in enumerate(self.reads) if index not in reads_to_remove]
        logging.debug("Finished removing masked reads.")

    """Post-processing: Use classification labels to remove sequences classified as human by deep learning."""
    def filter_reads(self,classification_labels):    
        logging.debug("Removing reads classified as human DNA...")
        for cur_tuple,classification in zip(self.read_list,classification_labels):
            # Sequence labeled as non-human, keep in filtered list
            if classification == 0:
                self.filtered_reads.append(cur_tuple)
        logging.debug("Finished removing reads classified as human DNA.")

    """Post-processing: Writes new fasta file based on filtered_list"""
    def write_fasta(self,output_file,filtered_list=self.filtered_reads):
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
    NN.load_fasta(args.input_fasta)
    NN.build_FFNN()

if __name__ == "__main__":
    main()
