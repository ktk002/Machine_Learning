#!/usr/bin/env python
# coding: utf8
import numpy as np
from keras.utils import np_utils
seqs = ["ACGTCCA","CGGATTG"]

CHARS = 'ACGT'
CHARS_COUNT = len(CHARS)

maxlen = max(map(len, seqs))
#res = np.zeros((len(seqs), CHARS_COUNT * maxlen), dtype=np.uint8)
#res = np.zeros((len(seqs), CHARS_COUNT * maxlen))
res = None
#for si, seq in enumerate(seqs):
#    seqlen = len(seq)
#    arr = np.chararray((seqlen,), buffer=seq)
#    for ii, char in enumerate(CHARS):
#        res[si][ii*seqlen:(ii+1)*seqlen][arr == char] = 1

#print(res)

# How to horizontally stack encoded numpy arrays with np_utils.to_categorical
for index,sequence in enumerate(seqs):
    new_string = list(map(int,[digit for digit in sequence.translate(str.maketrans("ACTG","0123"))]))
    one_hot = np.transpose(np_utils.to_categorical(new_string))
    for i,array in enumerate(one_hot):
        if i == 0:
            intermediate = array
        else:
            intermediate = np.hstack((intermediate,array))
    if index == 0:
        res = intermediate
    else:
        res = np.vstack((res,intermediate))
print(res)
