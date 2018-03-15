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
        self.seqs_to_remove = [] # headers of the sequences which will be removed
        self.filtered_reads = [] # filtered lists containing (header,sequence), of remaining reads

    """Create feed forward neural network models""" 
    def build_FFNN(self):
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

    """Convert ATCG to 1234 to one hot encoding"""
    def one_hot(self,string_to_convert):
        # Access Python 3's static str maketrans function
        numeric_string = string_to_convert.translate(str.maketrans("ACTG","1234"))
        new_list = list(map(int,[digit for digit in numeric_string]))

    """Load fasta file into read_dict property"""
    def load_fasta(self,fasta_file): 
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

    """Use classification labels to remove sequences classified as human by deep learning."""
    def filter_reads(self,classification_labels):    
        for cur_tuple,classification in zip(self.read_list,classification_labels):
            # Sequence labeled as non-human, keep in filtered list
            if classification == 0:
                self.filtered_reads.append(cur_tuple)
   
    """Writes new fasta file based on filtered_list"""
    def write_fasta(self,output_file,filtered_list=self.filtered_reads):
        with open(output_file,'w') as output_file:
            for cur_tuple in filtered_list:
                header = cur_tuple[0]
                sequence = cur_tuple[1]
                output_file.write("\n".join([header,sequence]))

def main():
    parser = ArgumentParser(description="Meta-Sweeper is a program to filter human DNA from metagenomic reads.")
    parser.add_argument('-i','--input_fasta',help="Specify input fasta file name to be filtered",
                        required=True,type=str)
    parser.add_argument('-o','--output_fasta',help="Specify final output fasta file name with filtered reads",
                        required=True,type=str)
    parser.add_argument('-d','--debug',help="Specify boolean to enable debug logging",
                        required=False,default=False,type=bool)
    args = parser.parse_args()   
    instance = NeuralNetwork()
    instance.load_fasta(args.input_fasta)
    #instance.build_FFNN()

if __name__ == "__main__":
    main()
