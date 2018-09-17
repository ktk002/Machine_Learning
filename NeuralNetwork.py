#!/usr/bin/env/ python

__author__ = 'Kellie Kim'

import tensorflow as tf
import numpy as np


class NeuralNetwork(object):
    training_data_types = ['intron', 'alu']
    NEGATIVE_CLASSIFICATION_LABEL = 0
    POSITIVE_CLASSIFICATION_LABEL_COUNTER = 0

    def __init__(self):
        self.x_train_list = []  # List of lists of ATCG values converted to 0-3 (ex: [[0, 1, 2, 3], [1, 1, 1, 1]])
 #       self.x_train_normalized = []  # Numpy list of lists containing normalized ATCG values (ex: [[0, 0.3333, 1]])
        self.y_train = []  # List of labels for every sample in x_train_list
        self.x_test = []
        self.y_test = []
        self.model = tf.keras.Sequential()
        self.predictions = []  # Classification predictions for unknown test data

    @staticmethod
    def dna_to_number(dna):
        """Convert DNA string to list of numbers ranging from 0-3."""
        return list(dna.translate(str.maketrans("ATCG", "0123")))

    def normalize_dna(self):
        """Takes x_train_list attribute and normalizes values in between 0 and 1."""
        self.x_train_list = tf.keras.utils.normalize(np.array(self.x_train_list), axis=1)

    def convert_to_numpy(self):
        """Converts training data and labels to numpy arrays."""
        self.x_train_list = np.array(self.x_train_list)
        self.y_train = np.array(self.y_train)

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

    def add_input_layer(self, num_features=100, num_neurons=267):
        """Add first layer of neural network and explicitly specify input shape of data."""
        self.model.add(tf.keras.layers.Dense(num_hidden_neurons=num_neurons, input_shape=(num_features,)))  # Add first input layer to Sequential model object

    def add_layers(self, num_layers=1, num_neurons=267, activation_function=tf.nn.relu):
        """Adds {num_layers} to model with {num_neurons} in each hidden layer and {activation_function}."""
        for index in range(num_layers):
            self.model.add(tf.keras.layers.Dense(num_neurons, activation=activation_function))

    def add_output_layer(self, activation_function=tf.nn.softmax, classification="binary"):
        """Adds final layer of neural network using number of classification labels."""
        if classification == "binary":
            num_nodes = 1  # Only 1 output required to make binary decision
            self.model.add(tf.keras.layers.Dense(num_nodes, activation=activation_function))
        else:
            #  Multiclassification, add nodes equal to number of classifications
            num_nodes = NeuralNetwork.POSITIVE_CLASSIFICATION_LABEL_COUNTER
            self.model.add(tf.keras.layers.Dense(num_nodes, activation=activation_function))

    def compile_model(self, optimizer='adam', loss='binary_crossentropy'):
        """Adds optimizer and loss function of neural network."""
        self.model.compile(optimizer=optimizer, loss=loss, metrics=['accuracy'])

    def fit_model(self, epochs=3):
        """Trains neural network model using {epochs} number of epochs."""
        self.model.fit(self.x_train_list, self.y_train, epochs=epochs)

    def print_validation_loss_accuracy(self):
        """Prints validation loss and validation accuracy from the neural network model with test data sets as input."""
        val_loss, val_acc = self.model.evaluate(self.x_test, self.y_test)
        print("validation loss is: ", val_loss)
        print("validation accuracy is: ", val_acc)

    def predict_unknown_datasets(self):
        """Use trained model to apply classification labels to unknown test data."""
        self.predictions = self.model.predict([self.x_test])

    def convert_classification_labels(self):
        """Takes argmax of each prediction to obtain final list of classification labels."""
        new_predictions = []
        for sample_number in range(self.predictions):
            new_predictions.append(np.argmax(self.predictions[sample_number]))
        self.predictions = new_predictions

    def build_nn(self, num_hidden_layers=1, num_hidden_neurons=267, activation="relu", classification="binary"):
        # Convert numerical values between 0 and 1 instead of 0-3
        self.normalize_dna()
        # After loading training data and normalizing, convert training data and labels to numpy arrays
        self.convert_to_numpy()
        self.add_input_layer()
        self.add_layers()
        self.add_output_layer()
        self.compile_model()
        self.fit_model()


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
        model1.load_training_data(input_file="", type="negative")

        # Vary number of hidden layers cases:
        # 1) 1 hidden layer, 267 hidden neurons
        # 2) 2 hidden layers, 267 hidden neurons
        # 3) 3 hidden layers, 267 hidden neurons
        # 4) 4 hidden layers, 267 hidden neurons
        # 5) 5 hidden layers, 267 hidden neurons
        # 6) 10 hidden layers, 267 hidden neurons
        # 7) 25 hidden layers, 267 hidden neurons
        # 8) 50 hidden layers, 267 hidden neurons
        # 9) 75 hidden layers, 267 hidden neurons
        # 10) 100 hidden layers, 267 hidden neurons
        num_hidden_layers_list = [1, 2, 3, 4, 5, 10, 25, 50, 75, 100]
        for hidden_layer_num in num_hidden_layers_list:
            model1.build_nn(num_hidden_layers=hidden_layer_num, num_hidden_neurons=267)

        # Vary number of hidden neurons cases:
        # 1) 1 hidden layer, 200 hidden neurons
        # 2) 1 hidden layer, 267 hidden neurons
        # 3) 1 hidden layer, 300 hidden neurons
        # 4) 1 hidden layer, 400 hidden neurons
        # 5) 1 hidden layer, 500 hidden neurons
        num_hidden_neurons_list = [200, 267, 300, 400, 500]
        for hidden_neuron_num in num_hidden_neurons_list:
            model1.build_nn(num_hidden_layers=1, num_hidden_neurons=hidden_neuron_num)

        # Vary activation function
        # 1) 1 hidden layer, 267 hidden neurons, rectified linear unit activation function (relu)
        # 2) 1 hidden layer, 267 hidden neurons, scaled exponential linear activation function (selu)-self normalizing
        # 3) 1 hidden layer, 267 hidden neurons, hyperbolic tangent activation function (tanh)
        # 4) 1 hidden layer, 267 hidden neurons, sigmoid activation function (sigmoid)
        # 5) 1 hidden layer, 267 hidden neurons, linear activation function (linear)
        activation_type_list = ["relu", "selu", "tanh", "sigmoid", "linear"]
        for activation_type in activation_type_list:
            model1.build_nn(num_hidden_layers=1, num_hidden_neurons=267, activation=activation_type)

    ##################################
    ### Multiclass classification ####
    ##################################


if __name__ == '__main__':
    main()


