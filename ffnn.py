#!/usr/bin/env python
import keras.utils
from keras.models import Sequential
from keras.layers import Dense, Activation
import numpy as np
#from string import maketrans

class NN():
    def __init__(self):
        self.read_dict = {}

    def build_FF(self):
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

        # 4. Linear activation function model (ie. no activation function)
        model4 = Sequential()
        model4.add(Dense(,activation='linear'))
        model4.compile()

        # 5. Scaled Exponential Linear Unit model(Klambauer et al., 2017)
        model5 = Sequential()
        model5.add(Dense(,activation='selu'))
        model5.compile()

        # Creates 10 patient row vectors, each with 5 features  
        #data = np.random.random((10,1))
        my_list = ["ACT","CTA"]
        data = np.array(my_list)
        print(data)
        # Makes 10 row vectors, each a dimension of 1, with a value in between 0 to 2, not including 2
        labels = np.random.randint(4, size=(10, 1)) 
        print(labels)
        one_hot_labels = keras.utils.to_categorical(labels,num_classes=64)
        #print(one_hot_labels)

    def one_hot(self,string_to_covert):
        # Convert ATCG to 1234
        # Access Python 3's static str maketrans function
        numeric_string = string_to_convert.translate(str.maketrans("ACTG","1234"))
        new_list = list(map(int,[digit for digit in numeric_string]))
        
def main():
    instance = NN()
    instance.build_FF()

if __name__ == "__main__":
    main()
