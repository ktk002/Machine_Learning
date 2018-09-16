#!/usr/bin/env/ python

import tensorflow as tf
import numpy as np


class NeuralNetwork(object):
    training_data_types = ['intron', 'alu']
    NEGATIVE_CLASSIFICATION_LABEL = 0
    POSITIVE_CLASSIFICATION_LABEL_COUNTER = 0

    def __init__(self):
        self.x_train_list = [] # List of lists of ATCG values converted to 0-3 (ex: [[0, 1, 2, 3], [1, 1, 1, 1]])
        self.x_train_normalized = [] # Numpy list of lists containing normalized ATCG values (ex: [[0, 0.3333, 1]])
        self.y_train = [] # List of labels for every sample in x_train_list
        self.x_test = []
        self.y_test = []
        self.classification_labels = [] # Labels ranging from 0 to number of classifications - 1

    @staticmethod
    def dna_to_number(dna):
        """Convert DNA string to list of numbers ranging from 0-3."""
        return list(dna.translate(str.maketrans("ATCG", "0123")))

    def load_training_data(self, input_file, type="positive", classification="binary"):
        """Loads data from input_file and adds as positive or negative data assuming binary classification."""
        # Increment to next largest classification label if classification is multi
        if classification == "multi":
            NeuralNetwork.POSITIVE_CLASSIFICATION_LABEL_COUNTER = NeuralNetwork.POSITIVE_CLASSIFICATION_LABEL_COUNTER + 1

        # Read in training data file
        with open(input_file, 'r') as input_file_reader:
            for line in input_file_reader:
                tokens = line.strip().split("\t")  # Split on tab, sequence is in second column
                sequence = tokens[1]  # Training sequence to be added to list
                # Add sequence to training data, regardless of label
                self.x_train_list.append(NeuralNetwork.dna_to_number(sequence))  # Add list containing numbers from 0-3

                # Negative data will be loaded the same for binary and multiclass classification
                if type == "negative":
                    self.y_train.append(NeuralNetwork.NEGATIVE_CLASSIFICATION_LABEL)  # Add label negative training data
                # Classification is multiclass, label of 0 is non-human dna, non zero is human dna
                elif type == "positive":
                    if classification == "binary":
                        self.y_train.append(1) # Add 1 to y_train
                    elif classification == "multi":
                        # Add training data to negative list
                        self.y_train.append(NeuralNetwork.POSITIVE_CLASSIFICATION_LABEL_COUNTER)
                    else:
                        print("Classification type invalid.")
                else:
                    print("Training data type invalid.")

    def build_nn(self, num_hidden_layers=1, num_hidden_neurons=267, activation="relu", classification="binary"):
        pass


def main():
    ##################################
    ###### Binary classification #####
    ##################################
    # Loop over each type of positive training data set

    # Construct neural networks for each type of training data
    for training_index, training_data_type in enumerate(NeuralNetwork.training_data_types):
        # Initialize new neural network
        model1 = NeuralNetwork()
        # First load positive and negative training data into neural network
        model1.load_training_data(input_file="all_alu.tsv", type="positive")
#        model1.load_training_data_binary(input_file="", type="negative")

        # Vary number of hidden layers
        # 1) 1 hidden layer, 267 hidden neurons
        model1.build_nn(num_hidden_layers=1, num_hidden_neurons=267)
        # 2) 2 hidden layers, 267 hidden neurons
        model1.build_nn(num_hidden_layers=2, num_hidden_neurons=267)
        # 3) 3 hidden layers, 267 hidden neurons
        model1.build_nn(num_hidden_layers=3, num_hidden_neurons=267)
        # 4) 4 hidden layers, 267 hidden neurons
        model1.build_nn(num_hidden_layers=4, num_hidden_neurons=267)
        # 5) 5 hidden layers, 267 hidden neurons
        model1.build_nn(num_hidden_layers=5, num_hidden_neurons=267)
        # 6) 10 hidden layers, 267 hidden neurons
        model1.build_nn(num_hidden_layers=10, num_hidden_neurons=267)
        # 7) 25 hidden layers, 267 hidden neurons
        model1.build_nn(num_hidden_layers=25, num_hidden_neurons=267)
        # 8) 50 hidden layers, 267 hidden neurons
        model1.build_nn(num_hidden_layers=50, num_hidden_neurons=267)
        # 9) 75 hidden layers, 267 hidden neurons
        model1.build_nn(num_hidden_layers=75, num_hidden_neurons=267)
        # 10) 100 hidden layers, 267 hidden neurons
        model1.build_nn(num_hidden_layers=100, num_hidden_neurons=267)

        # Vary number of hidden neurons
        # 1) 1 hidden layer, 200 hidden neurons
        model1.build_nn(num_hidden_layers=1, num_hidden_neurons=200)
        # 2) 1 hidden layer, 267 hidden neurons
        model1.build_nn(num_hidden_layers=1, num_hidden_neurons=267)
        # 3) 1 hidden layer, 300 hidden neurons
        model1.build_nn(num_hidden_layers=1, num_hidden_neurons=300)
        # 4) 1 hidden layer, 400 hidden neurons
        model1.build_nn(num_hidden_layers=1, num_hidden_neurons=400)
        # 5) 1 hidden layer, 500 hidden neurons
        model1.build_nn(num_hidden_layers=1, num_hidden_neurons=500)

        # Vary activation function
        # 1) 1 hidden layer, 267 hidden neurons, rectified linear unit activation function (relu)
        model1.build_nn(num_hidden_layers=1, num_hidden_neurons=267, activation="relu")
        # 2) 1 hidden layer, 267 hidden neurons, scaled exponential linear activation function (selu)-self normalizing
        model1.build_nn(num_hidden_layers=1, num_hidden_neurons=267, activation="selu")
        # 3) 1 hidden layer, 267 hidden neurons, hyperbolic tangent activation function (tanh)
        model1.build_nn(num_hidden_layers=1, num_hidden_neurons=267, activation="tanh")
        # 4) 1 hidden layer, 267 hidden neurons, sigmoid activation function (sigmoid)
        model1.build_nn(num_hidden_layers=1, num_hidden_neurons=267, activation="sigmoid")
        # 5) 1 hidden layer, 267 hidden neurons, linear activation function (linear)
        model1.build_nn(num_hidden_layers=1, num_hidden_neurons=267, activation="linear")

    ##################################
    ### Multiclass classification ####
    ##################################


if __name__ == '__main__':
    main()


