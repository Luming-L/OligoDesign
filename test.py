# -*- coding: UTF-8 -*-
def f1(seq_range=[0, 8], min=2, max=3, tem=[], result=[]):
    if seq_range[1] + 1 - seq_range[0] >= min and seq_range[1] + 1 - seq_range[0] <= max:
        result.append(tem + [seq_range])
    elif seq_range[1] + 1 - seq_range[0] == 0:
        result.append(tem)
    for i in range(min, max+1, 1):
        for j in range(-min, -max-1, -1):
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

if __name__ == '__main__':
    x = "AGGTCTCTCTCTCTCTGGTACCAATCTAGAGGATCCCTCGAGGATTATGTGGAAAAAAAGCACCGACTCGGTGCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTCAACTTGCTATGCTGTTTCCAGCATAGCTCTGAAACTGAGGCAGGCGGGGATGAAGTGCCACGGATCATCTGCACAACTCTTTTAAATCAGCTTTGATCTATGTGGATAGCCGAGGTAGAGACC"
    results1 = f1(seq_range=[0, 8], min=2, max=3, tem=[], result=[])
    results2 = list2dict(input=results1)

    for i in range(len(results2)):
        print "Group", i + 1, results2[i]




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

