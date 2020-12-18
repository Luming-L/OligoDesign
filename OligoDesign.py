# -*- coding: UTF-8 -*-
"""
    The method (invented by LJW) can synthesize nearly any DNA sequence, including sequences with high GC content,
    low GC content and repeated regions. The synthesizing system contains an initial double-strand sequence, a group of
    oligos, enzymes (DNA ligase, DNA polymerase and restriction enzyme) and buffer.

    This script aims to design oligos. The output oligos have hairpin structure and restriction site.

    Design:
    - The composition of a oligo: part of the targeted sequence (main body of each oligo),
    restriction site and reverse complement sequence (2-10 bases overhang in 3' end)

    - main body length range (oligo)
    - main body of oligos overlap with each other. the length of the overlap is n.
    - no restriction site/avoid site in the main body

    - length and position of reverse complement sequence
    - unique overhang of each oligo in a group
    - deltaG



    6. A/B/C/D...片段之间有唯一重叠末端，可以连接

    ##### 拆成200bp
    ##### 对于200bp，添加wrap5和wrap3，拆分成60-75bp的片段（递归）
    ##### 每拆分出一组片段，判断是否符合要求，不符合要求继续递归
    ##### 要求包括：不包含酶切位点，唯一粘性末端
    ##### 考虑division point，构造hairpin结构

    Input:
        original sequence
        flanking sequence
        length of sub sequence
        length of oligos (max,min)

    Variables:
        strictly monotonically increasing sequences of indices: i1, i1', i2, i2', ..., ix, ix'
            ik - i(k-1)' = 2
            i1 = 1, ix = len(subSeq)
            1 <= k <= x
        group size
        division point
        restriction enzyme (i,)

    Output:
        position of oligos' ends
"""

# import
import re
from Bio.Seq import Seq


class OligoGroup():
    """a group of oligos for the synthesis of a sub sequence"""

    def __init__(self, subSeq, wrap5='AGGTCTCT', wrap3='AGAGACC'):
        """Initialize oligoGroup attributes."""
        self.subSeq = Seq(subSeq)
        self.wrap5 = Seq(wrap5)
        self.wrap3 = Seq(wrap3)
        self.wrap_seq = wrap5 + subSeq + wrap3

        # Other parameters
        self.groupSize = None
        self.divisionPoint = None
        self.len_oligo_min = 30
        self.len_stem = 14
        self.restrictionEnzyme = ""
        self.oligoGroup = {}  # {"oligo1":"ATCG",...}

    def restriction(self, wrap_seq):
        """determine the restriction enzyme, and get positions of restriction site"""
        restrictionSites = {'BtsI': Seq('GCAGTG'), 'BsrDI': Seq('GCAATG')}
        BtsI_num = len(
            list(re.finditer(r'(?=GCAGTG|GTGACG|CACTGC|CGTCAC)',
                             str(oriSeq))))  # BtsI, oriSeq必须是string，pattern最好是raw string
        BsrDI_num = len(list(re.finditer(r'(?=GCAATG|GTAACG|CATTGC|CGTTAC)', str(oriSeq))))  # BsrDI
        print BtsI_num
        print BsrDI_num

        for i in range(len(oligoGroups_starts) - 1):
            print "i", i
            # the first round
            if not BtsI_num and not BsrDI_num:  # no restriction sites
                print "no restriction sites"
            else:
                for site in re.finditer(r'(?=GCAGTG|GTGACG|CACTGC|CGTCAC)',
                                        str(oriSeq[oligoGroups_starts[i]:oligoGroups_starts[i + 1]])):
                    BtsI_sites_starts.append(site.start())
                for site in re.finditer(r'(?=GCAATG|GTAACG|CATTGC|CGTTAC)',
                                        str(oriSeq[oligoGroups_starts[i]:oligoGroups_starts[i + 1]])):
                    BsrDI_sites_starts.append(site.start())
                if BtsI_num >= BsrDI_num:
                    restrictionSites_starts = BtsI_sites_starts
                else:
                    restrictionSites_starts = BsrDI_sites_starts
                print "restrictionSites_starts", restrictionSites_starts
                print
        return


# Input
with open("/home/geneland/PycharmProjects/OligoDesign/oriSeq", "r") as f:
    oriSeq = Seq(f.read())
len_oriSeq = len(oriSeq)


len_restrictionSite = 6  # restriction recognition sites (BtsI or BsrDI)
len_sticky = 2
len_oriSeqInOligo = len_oligo - len_stem - len_restrictionSite

BtsI_sites_starts = []
BsrDI_sites_starts = []
restrictionSites_starts = []
oligoGroups_starts = []
oligoPositions = {}

startPoint = 0

if len_oriSeq <= 200:
    groups = 1
else:
    groups = int(len_oriSeq / 200)
for i in range(groups):
    oligoGroups_starts.append(i * 200)
oligoGroups_starts.append(len_oriSeq + 1)

print "oligoGroups_starts", len(oligoGroups_starts), oligoGroups_starts

# a possible way get oligos
for i in range(groupSize):
    if i + 1 == 1:
        oligoPositions["oligo" + str(i + 1)] = startPoint + len_oriSeqInOligo - len(wrapMode["wrap5"])
    elif i + 1 > 1 and i + 1 < groupSize:
        oligoPositions["oligo" + str(i + 1)] = oligoPositions["oligo" + str(i)] + len_oriSeqInOligo
    elif i + 1 == groupSize:
        oligoPositions["oligo" + str(i + 1)] = oligoPositions["oligo" + str(i)] + len_oriSeqInOligo - len(
            wrapMode["wrap3"])

# search the original sequence for BtsI and BsrDI restriction sites
# seq_BtsI = [str(restrictionSites["BtsI"]), str(restrictionSites["BtsI"][::-1]),
#             str(restrictionSites['BtsI'].reverse_complement()),
#             str(restrictionSites['BtsI'].reverse_complement()[::-1])]
# seq_BsrDI = [str(restrictionSites["BsrDI"]), str(restrictionSites["BsrDI"][::-1]),
#              str(restrictionSites['BsrDI'].reverse_complement()),
#              str(restrictionSites['BsrDI'].reverse_complement()[::-1])]
# print seq_BtsI
# print seq_BsrDI



# Output
print oriSeq
print oligoPositions

