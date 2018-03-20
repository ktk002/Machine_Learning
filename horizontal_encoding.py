
import numpy as np

seqs = ["ACGTCCA","CGGATTG"]

CHARS = 'ACGT'
CHARS_COUNT = len(CHARS)

maxlen = max(map(len, seqs))
#res = np.zeros((len(seqs), CHARS_COUNT * maxlen), dtype=np.uint8)
res = np.zeros((len(seqs), CHARS_COUNT * maxlen))

for si, seq in enumerate(seqs):
    seqlen = len(seq)
    arr = np.chararray((seqlen,), buffer=seq)
    for ii, char in enumerate(CHARS):
        res[si][ii*seqlen:(ii+1)*seqlen][arr == char] = 1

print(res)
