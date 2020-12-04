"""
    Input:
        original sequence
        length of oligos
        group size
        division point
        restriction enzyme

    Output:
        position of oligos' ends

"""
# import
import re
from Bio.Seq import Seq
# Input
with open("/home/geneland/PycharmProjects/OligoDesign/oriSeq", "r") as f:
    oriSeq = Seq(f.read())
groupSize = 5
divPoint = 3
len_stem = 14
len_resRec = 6
len_oligo = 70
len_sticky =2
len_oriSeqInOligo = len_oligo - len_stem - len_resRec
startPoint = 0
resRec = {'BtsI': Seq('GCAGTG'), 'BsrDI': Seq('GCAATG')}
wrapMode = {"wrap5": "AGGTCTCT", "wrap3": "AGAGACC"}

oligoPos = {}
for i in range(groupSize):
    if i + 1 == 1:
        oligoPos["oligo" + str(i + 1)] = startPoint + len_oriSeqInOligo - len(wrapMode["wrap5"])
    elif i + 1 > 1 and i + 1 < groupSize:
        oligoPos["oligo" + str(i + 1)] = oligoPos["oligo" + str(i)] + len_oriSeqInOligo
    elif i + 1 == groupSize:
        oligoPos["oligo" + str(i + 1)] = oligoPos["oligo" + str(i)] + len_oriSeqInOligo - len(wrapMode["wrap3"])

# Output
print oriSeq
print oligoPos
print resRec['BsrDI'], resRec['BsrDI'][::-1], resRec['BsrDI'].reverse_complement(), resRec['BsrDI'].reverse_complement()[::-1]
print resRec['BtsI'], resRec['BtsI'].reverse_complement()