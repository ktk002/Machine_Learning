#!/usr/bin/env/ python

__author__ = 'Kellie Kim'

import tensorflow as tf
import numpy as np
import random

from sklearn.metrics import confusion_matrix
from tensorflow.python.keras.optimizers import SGD

class NeuralNetwork(object):
  #  training_data_types = ['intron', 'alu']
    training_data_types = ['intron']
    NEGATIVE_CLASSIFICATION_LABEL = 0

    def __init__(self):
        self.x_train_list = []  # List containing list of lists of ATCG values converted to 0-3 (ex: [[0, 1, 2, 3], [1, 1, 1, 1]])
 #       self.x_train_normalized = []  # Numpy list of lists containing normalized ATCG values (ex: [[0, 0.3333, 1]])
        self.y_train = []  # Numpy array of labels for every sample in x_train_list; length = # of total samples
        self.x_test = []
        self.y_test = []
        self.model = tf.keras.Sequential()
        self.predictions = []  # Classification predictions for unknown test data

    @staticmethod
    def dna_to_number(dna):
        """Convert DNA string to list of numbers ranging from 0-3. Assumes only 4 letters in sequence"""
        return list(map(int, dna.translate(str.maketrans("ATCG", "0123"))))

    def get_shortest_list_length(self, list_of_lists):
        """Returns the list of the shortest length from a set of lists."""
        return min(map(len, list_of_lists))

    def normalize_dna(self, data_type="training"):
        """Takes x_train_list attribute and normalizes values in between 0 and 1."""
        if data_type == "training":
            #print("normalizing: ", self.x_train_list)
            self.x_train_list = tf.keras.utils.normalize(np.array(self.x_train_list), axis=1)
        else:
            print("normalizing: ", self.x_test)
            self.x_test = tf.keras.utils.normalize(np.array(self.x_test), axis=1)

    def convert_to_numpy(self):
        """Converts training data and labels to numpy arrays."""
        print(self.x_train_list)
        self.x_train_list = np.array(self.x_train_list)
        self.y_train = np.array(self.y_train)

    def load_data(self, input_file, training_type="negative", data_type="training", classification="binary"):
        """
        Loads data from input_file and adds as positive or negative data assuming binary classification.
        Assumes that negative data will always be loaded first and associated with label 0.
        Positive human data is labeled with values other than 0 such as 1.
        Assumes data_type will always only be either training or testing.
        """
        cur_list = []
        print("loading training type: ", training_type)
        print("loading data type: ", data_type)

        # Read in training data file
        with open(input_file, 'r') as input_file_reader:
            for line in input_file_reader:
                tokens = line.strip().split("\t")  # Split on tab, sequence identifer in first column, sequence is in second column
                sequence = tokens[1]  # Training sequence to be added to list
                # Add sequence to training data, regardless of label
                cur_list.append(NeuralNetwork.dna_to_number(sequence))  # Add list containing numbers from 0-3

        if data_type == "training":
            self.x_train_list.append(cur_list)
        else:
            self.x_test.append(cur_list)

    def equalize_training_data(self):
        """Sets negative and positive training data lengths to be equal or every classification to be of equal length to minimize training bias."""
        # Set length of all classifications to be equal to shortest list
        # Get length of smallest list
        shortest_length = self.get_shortest_list_length(self.x_train_list)
        x_train_copy = []

        # Make all classification lists equal to the shortest length by randomly choosing n elements from list.
        for training_list in self.x_train_list:
            # Cut list to shortest list length
            random_sample = random.sample(training_list, shortest_length)
            x_train_copy.append(random_sample)
        self.x_train_list = x_train_copy

    def set_labels(self, data_type="training"):
        """Assigns labels to y_train list attribute based on lengths of lists in x_train or x_test.
            Default data_type is training.
            Assumes data_type can only be training or testing.
        """
        y_copy = []
        x_list = None

        # Smallest number of samples in x_train_list attribute
#        min_num_samples = self.get_shortest_list_length(self.x_train_list)
        if data_type == "training":
            x_list = self.x_train_list
        else:
            x_list = self.x_test

        # Iterate over either training or testing list and add either 0 or 1 for the number of elements in the list
        for index, training_list in enumerate(x_list):
            print("adding number of training samples: ", len(training_list))
            y_copy = y_copy + [index] * len(training_list) # Use full number of negative or positive training samples, don't cut to smaller amount
  #          y_train_copy = y_train_copy + [index] * min_num_samples # equal number of negative and positive training samples

        if data_type == "training":
            self.y_train = np.array(y_copy)
        else:
            self.y_test = np.array(y_copy)

    def merge_lists(self, data_type="training"):
        """Merge list of lists in x training data into a single list after assigning all labels."""
        new_list = []

        if data_type == "training":
            x_list = self.x_train_list
        else:
            x_list = self.x_test
        for cur_list in x_list:
            new_list = new_list + [inner_list for inner_list in cur_list]

        if data_type == "training":
            self.x_train_list = new_list
        else:
            self.x_test = new_list

    def add_input_layer(self, num_features=100, num_neurons=267):
        """Add first layer of neural network and explicitly specify input shape of data."""
        self.model.add(tf.keras.layers.Dense(num_neurons, input_shape=(num_features,)))  # Add first input layer to Sequential model object

    def add_layers(self, num_layers=1, num_neurons=267, activation_function=tf.nn.relu):
        """Adds {num_layers} to model with {num_neurons} in each hidden layer and {activation_function}."""
        for index in range(num_layers):
            self.model.add(tf.keras.layers.Dense(num_neurons, activation=activation_function))

    def add_output_layer(self, activation_function=tf.nn.softmax, classification="binary"):
        """Adds final layer of neural network using number of classification labels."""
        if classification == "binary":
            num_nodes = 1  # Only 1 output required to make binary decision
            print("using sigmoid")
            self.model.add(tf.keras.layers.Dense(num_nodes, activation='sigmoid'))
        else:
            #  Multiclassification, add nodes equal to number of classifications
            num_nodes = NeuralNetwork.POSITIVE_CLASSIFICATION_LABEL_COUNTER
            self.model.add(tf.keras.layers.Dense(num_nodes, activation=activation_function))

    def compile_model(self, optimizer='adam', loss='binary_crossentropy', classification="binary"):
        """Adds optimizer and loss function of neural network."""
        if classification == "binary":
            self.model.compile(optimizer=optimizer, loss=loss, metrics=['accuracy'])
        elif classification == "multi":
            self.model.compile(optimizer=optimizer, loss='categorical_crossentropy', metrics=['accuracy'])
        else:  # custom compilation
            print("Invalid classification type during model compilation.")

    def fit_model(self, epochs=3):
        """Trains neural network model using {epochs} number of epochs."""
        self.model.fit(self.x_train_list, self.y_train, epochs=epochs)


    # Functions to calculate statistics
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
        new_predictions = [round(x[0]) for x in self.predictions]
#        for sample_number in range(len(self.predictions)):
#            new_predictions.append(np.argmax(self.predictions[sample_number]))
        print("predictions list is: ")
        print(new_predictions)
        print("old predictions list: ")
        print(self.predictions)
        self.predictions = np.array(new_predictions) # Convert list of python labels to a numpy array

    def get_confusion_matrix(self, expected_labels, predicted_labels):
        """
        Takes array of expected (self.y_test) and array of actual classification labels (self.predictions) of the same length and returns confusion matrix.
        Returns values: true positive, true negative, false positive, and false negative
        """
        print("type and length of expected labels is: ", type(expected_labels), len(expected_labels))
        print("type of predicted labels is: ", type(predicted_labels), len(predicted_labels))
        calculated_matrix = confusion_matrix(expected_labels, predicted_labels)
        TP = calculated_matrix[1][1] # Values of 1 in the matrix
        TN = calculated_matrix[0][0] # Values of 0 in the matrix
        FP = calculated_matrix[0][1]
        FN = calculated_matrix[1][0]

        print("TP is: ", TP)
        print("TN is: ", TN)
        print("FP is: ", FP)
        print("FN is: ", FN)
        return TP, TN, FP, FN

    def calculate_metrics(self, TP, TN, FP, FN):
        """
        Returns 4 values: sensitivity and precision, given values from confusion matrix.
        Equations used from: https://www.wikihow.com/Calculate-Sensitivity%2C-Specificity%2C-Positive-Predictive-Value%2C-and-Negative-Predictive-Value
        """
        sensitivity = (TP/(TP + FN)) # How likely test will detect positive human DNA (1)
        precision = (TN/(FP + TN)) # How likely test will detect negative human DNA (0)
        positive_predictive_value = (TP/(TP + FP)) # How likely DNA will be human if label is positive (1)
        negative_predictive_value = (TN/ (TN + FN)) # How likely DNA will not be human if label is negative (0)

        return sensitivity, precision, positive_predictive_value, negative_predictive_value

    def build_nn(self, num_hidden_layers=1, num_hidden_neurons=267, activation=tf.nn.relu, classification="binary"):
        # Convert numerical values between 0 and 1 instead of 0-3
        self.normalize_dna()
        # After loading training data and normalizing, convert training data and labels to numpy arrays
        self.convert_to_numpy()
        print("adding input layer")
        self.add_input_layer()
        self.add_layers(num_layers=num_hidden_layers, num_neurons=num_hidden_neurons, activation_function=activation)
        self.add_output_layer(classification=classification)

        # Set optimizer and learning rate
        self.compile_model(optimizer=SGD(lr=0.001))
        self.fit_model()


def main():
    ##################################
    ###### Binary classification #####
    ##################################
    # Loop over each type of positive training data set
    # Construct neural networks for each type of training data (options are introns and alus)
    for training_index, training_data_type in enumerate(NeuralNetwork.training_data_types):
        # Initialize new neural network
        model1 = NeuralNetwork()
        # First load positive and negative training data into neural network
  #      model1.load_data(input_file="alu_sliding_1.tsv", training_type="positive")
  #      model1.load_data(input_file="merged_alus.tsv", training_type="positive")
        model1.load_data(input_file="introns_hg38.tsv", training_type="positive", data_type="training")
   #     model1.load_data(input_file="fake_negative.tsv", training_type="negative")
        model1.load_data(input_file="all_negative.tsv", training_type="negative", data_type="training")

        # Load testing data sets
        model1.load_data(input_file="500_positive_training.tsv", training_type="positive", data_type="testing")
        model1.load_data(input_file="500_negative_training.tsv", training_type="negative", data_type="testing")

        # After loading all training data, cut all lists in x_train_list to all be the length of shortest list
#        model1.equalize_training_data()

        # Add classification labels to y_train and y_test attributes
        model1.set_labels(data_type="training")
        model1.set_labels(data_type="testing")

        # Merge x_train_list and x_test into single list of lists to be normalized and processed
        model1.merge_lists(data_type="training")
        model1.merge_lists(data_type="testing")

        # Normalize DNA in x_train_list to be values between 0 and 1
        model1.normalize_dna(data_type="training")
        model1.normalize_dna(data_type="testing")

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
 #       num_hidden_layers_list = [1, 2, 3, 4, 5, 10, 25, 50, 75, 100]
        num_hidden_layers_list = [1, 10, 20] # 1, 10, 100
        for hidden_layer_num in num_hidden_layers_list:
            print("number of hidden layers: ", hidden_layer_num)
            model1.build_nn(num_hidden_layers=1, num_hidden_neurons=267, activation="relu")
            model1.print_validation_loss_accuracy()
            model1.predict_unknown_datasets()
            model1.convert_classification_labels()
            TP, TN, FP, FN = model1.get_confusion_matrix(model1.y_test, model1.predictions)
            sensitivity, precision, positive_predictive_value, negative_predictive_value = model1.calculate_metrics(TP, TN, FP, FN)
            print("sensitivity is: ", sensitivity)
            print("precision is: ", precision)
            print("positive_predictive value is: ", positive_predictive_value)
            print("negative predictive value is: ", negative_predictive_value)
 #           print("number of hidden layers is: ", hidden_layer_num)
  #          model1.build_nn(num_hidden_layers=hidden_layer_num, num_hidden_neurons=267)

        # Vary number of hidden neurons cases:
        # 1) 1 hidden layer, 200 hidden neurons
        # 2) 1 hidden layer, 267 hidden neurons
        # 3) 1 hidden layer, 300 hidden neurons
        # 4) 1 hidden layer, 400 hidden neurons
        # 5) 1 hidden layer, 500 hidden neurons
  #      num_hidden_neurons_list = [200, 267, 300, 400, 500]
   #     for hidden_neuron_num in num_hidden_neurons_list:
   #         print("hidden neuron number is: ", hidden_neuron_num)
    #        model1.build_nn(num_hidden_layers=1, num_hidden_neurons=hidden_neuron_num)

        # Vary activation function
        # 1) 1 hidden layer, 267 hidden neurons, rectified linear unit activation function (relu)
        # 2) 1 hidden layer, 267 hidden neurons, scaled exponential linear activation function (selu)-self normalizing
        # 3) 1 hidden layer, 267 hidden neurons, hyperbolic tangent activation function (tanh)
        # 4) 1 hidden layer, 267 hidden neurons, sigmoid activation function (sigmoid)
        # 5) 1 hidden layer, 267 hidden neurons, linear activation function (linear)
        # activation_type_list = ["relu", "selu", "tanh", "sigmoid", "linear"]
        # for activation_type in activation_type_list:
        #     print("activation type is: ", activation_type)
        #     model1.build_nn(num_hidden_layers=1, num_hidden_neurons=267, activation=activation_type)
        #     model1.print_validation_loss_accuracy()
        #     model1.predict_unknown_datasets()
        #     model1.convert_classification_labels()
        #     TP, TN, FP, FN = model1.get_confusion_matrix(model1.y_test, model1.predictions)
        #     sensitivity, precision, positive_predictive_value, negative_predictive_value = model1.calculate_metrics(TP, TN, FP, FN)
        #     print("sensitivity is: ", sensitivity)
        #     print("precision is: ", precision)
        #     print("positive_predictive value is: ", positive_predictive_value)
        #     print("negative predictive value is: ", negative_predictive_value)

    ##################################
    ### Multiclass classification ####
    ##################################


if __name__ == '__main__':
    main()


