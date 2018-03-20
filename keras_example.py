#!/usr/bin/env python
import numpy as np
from keras.models import Sequential
from keras.layers import Dense,Activation,Dropout

np.random.seed(1)

x_train = np.random.random((1000,20)) # number patients = 1000, number of features = 20
y_train = np.random.randint(2,size=(1000,1))
x_test = np.random.random((100,20))
y_test = np.random.randint(2,size=(100,1))

model = Sequential()
# Note input_dim (number of features) only needs to be specified once for the input layer and not for the additional layers
model.add(Dense(units=64,activation='relu',input_dim=20)) # units = output_dim's, input_dim = number of features/dimensionality of each patient vector
model.add(Dropout(0.5))
model.add(Dense(units=64,activation='relu')) # units = number of neurons in the layer
model.add(Dropout(0.5))
model.add(Dense(units=1,activation='sigmoid')) # 1 output node for binary classification
model.compile(optimizer='rmsprop',loss='binary_crossentropy',metrics=['accuracy'])

model.fit(x_train,y_train,epochs=20,batch_size=128,verbose=1,shuffle=True,validation_data=(x_test,y_test)) # shuffles data, print verbose info, and gives validation set
score = model.evaluate(x_test,y_test,batch_size=128)

print("Test score: " + str(score[0]))
print("Test accuracy: " + str(score[1]))
