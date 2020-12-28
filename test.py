# -*- coding: UTF-8 -*-


# import
import re
from Bio.Seq import Seq


REases = {"BtsI": 'GCAGTG', "BsrDI": 'GCAATG'}



def f1(seq_range=[0, 8], min=2, max=3, tem=[], result=[]):
    if seq_range[1] + 1 - seq_range[0] >= min and seq_range[1] + 1 - seq_range[0] <= max:
        result.append(tem + [seq_range])
    elif seq_range[1] + 1 - seq_range[0] == 0:
        result.append(tem)
    for i in range(min, max + 1, 1):
        for j in range(-min, -max - 1, -1):
            oligos = [[seq_range[0], seq_range[0] + i - 1], [seq_range[1] + j + 1, seq_range[1]]]
            remain = [seq_range[0] + i, seq_range[1] + j]
            if remain[1] + 1 - remain[0] >= min or remain[1] + 1 - remain[0] == 0:
                f1(seq_range=remain, min=min, max=max, tem=tem + oligos, result=result)
    return result


def list2dict(input=[]):
    output = []
    for i in range(len(input)):
        tem = {}
        for j in range(len(input[i])):
            tem[input[i][j][0]] = input[i][j]
        for index in range(1, len(tem.keys()) + 1):
            tem[index] = tem.pop(sorted(tem.keys())[index - 1])
        output.append(tem)
    return output


def determine_REase(sequence, group, default_REase="BtsI", REases=REases):
    """determine the restriction enzyme"""

    overlap = 1
    for REase in REases:
        for i in range(1, len(group) + 1):
            if i == len(group):
                sub = sequence[group[i][0], group[i][1]]
            else:
                sub = sequence[group[i][0], group[i][1] + overlap]
            if 存在酶切位点:
                break
        enzyme = REase
        break

    return enzyme
    # BtsI_num = len(
    #     list(re.finditer(r'(?=GCAGTG|GTGACG|CACTGC|CGTCAC)',
    #                      str(oriSeq))))  # BtsI, oriSeq必须是string，pattern最好是raw string
    # BsrDI_num = len(list(re.finditer(r'(?=GCAATG|GTAACG|CATTGC|CGTTAC)', str(oriSeq))))  # BsrDI
    # print BtsI_num
    # print BsrDI_num
    #
    # for i in range(len(oligoGroups_starts) - 1):
    #     print "i", i
    #     # the first round
    #     if not BtsI_num and not BsrDI_num:  # no restriction sites
    #         print "no restriction sites"
    #     else:
    #         for site in re.finditer(r'(?=GCAGTG|GTGACG|CACTGC|CGTCAC)',
    #                                 str(oriSeq[oligoGroups_starts[i]:oligoGroups_starts[i + 1]])):
    #             BtsI_sites_starts.append(site.start())
    #         for site in re.finditer(r'(?=GCAATG|GTAACG|CATTGC|CGTTAC)',
    #                                 str(oriSeq[oligoGroups_starts[i]:oligoGroups_starts[i + 1]])):
    #             BsrDI_sites_starts.append(site.start())
    #         if BtsI_num >= BsrDI_num:
    #             restrictionSites_starts = BtsI_sites_starts
    #         else:
    #             restrictionSites_starts = BsrDI_sites_starts
    #         print "restrictionSites_starts", restrictionSites_starts
    #         print
    # return


if __name__ == '__main__':
    x = "AGGTCTCTCTCTCTCTGGTACCAATCTAGAGGATCCCTCGAGGATTATGTGGAAAAAAAGCACCGACTCGGTGCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTCAACTTGCTATGCTGTTTCCAGCATAGCTCTGAAACTGAGGCAGGCGGGGATGAAGTGCCACGGATCATCTGCACAACTCTTTTAAATCAGCTTTGATCTATGTGGATAGCCGAGGTAGAGACC"
    results1 = f1(seq_range=[0, 8], min=2, max=3, tem=[], result=[])
    results2 = list2dict(input=results1)

    oligoGroup = results2[0]

    print oligoGroup

# def f(sequence="abc"*3,  tem="", min=2, max=3):
#     if (len(sequence) >= min and len(sequence) <= max) or len(sequence) == 0:
#         print tem + "(" + sequence + ")"
#     for i in range(min, max+1, 1):
#         for j in range(-min, -max-1, -1):
#             oligo1 = sequence[:i]
#             oligo2 = sequence[j:]
#             remain = sequence[i:j]
#             x = tem + "(" + oligo1 + ")" + "(" + oligo2 + ")"
#             if len(remain) >= min or i + abs(j) == len(sequence):
#                 f(sequence=remain, tem=x, min=min, max=max)
#
# def f2(sequence="abc"*3,  tem=[], min=2, max=3):
#     if len(sequence) >= min and len(sequence) <= max:
#         print tem + [sequence]
#     elif len(sequence) == 0:
#         print tem
#
#     for i in range(min, max+1, 1):
#         for j in range(-min, -max-1, -1):
#             oligos = [sequence[:i], sequence[j:]]
#             remain = sequence[i:j]
#             if len(remain) >= min or i + abs(j) == len(sequence):
#                 f2(sequence=remain, tem=tem + oligos, min=min, max=max)
