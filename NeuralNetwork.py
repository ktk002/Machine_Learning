#!/usr/bin/env/ python

import tensorflow as tf
import numpy as np

training_data_types = ['intron', 'alu']

class NeuralNetwork(object):
    def __init__(self):
        self.x_train_list = [] # List of lists of ATCG values converted to 0-3 (ex: [[0, 1, 2, 3], [1, 1, 1, 1]])
        self.x_train_normalized = DNA_to_ # List of lists containing normalized ATCG values (ex: [[0, 0.3333, 1],[1, 1 ,1]])
        self.y_train = [] # List of labels for every sample in x_train_list
        self.x_test = []
        self.y_test = []
        self.build_nn()

    @staticmethod
    def DNA_to_number(DNA):
        """Convert DNA string to list of numbers ranging from 0-3."""
        return list(DNA.translate(str.maketrans("ATCG", "0123")))

    def load_training_data(self, input_file, type="positive"):
        """Loads data from input_file and adds as positive or negative data."""
        with open(input_file, 'r') as input_file_reader:
            for line in input_file_reader:
                tokens = line.strip().split()
                identifier = tokens[0]
                sequence = tokens[1]

                # Add sequence depending on label type and add label
                if type == "positive":
                    self.x_train_list.append(DNA_to_number(sequence))
                    self.y

    def build_nn(self, num_hidden_layers=1, num_hidden_neurons=267, activation="relu", classification="binary"):
        pass

def main():
    ##################################
    ###### Binary classification #####
    ##################################
    # Loop over each type of positive training data set

    # Construct neural networks for each type of training data
    for training_index, training_data_type in enumerate(training_data_types):
        # Initialize new neural network
        model1 = NeuralNetwork()
        # First load positive and negative training data into neural network
        model1.load_training_data(input_file=, type="positive")
        model1.load_training_data(input_file=e, type="negative")

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
        # 2) 1 hidden layer, 267 hidden neurons, scaled exponential linear activation function (selu) --> self normalizing
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


