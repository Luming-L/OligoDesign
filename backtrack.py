# coding=utf-8


# version 1.1
# recursion, print results as strings
def break_a_sequence_v1_1(sequence, tem, minimum, maximum):
    if minimum <= len(sequence) <= maximum:
        print tem + "(" + sequence + ")"
    if len(sequence) == 0:
        print tem
    for i in range(minimum, maximum + 1, 1):
        for j in range(-minimum, -maximum - 1, -1):
            subSeq1 = sequence[:i]
            subSeq2 = sequence[j:]
            remain = sequence[i:j]
            x = tem + "(" + subSeq1 + ")" + "(" + subSeq2 + ")"
            if len(remain) >= minimum or i + abs(j) == len(sequence):
                break_a_sequence_v1_1(sequence=remain, tem=x, minimum=minimum, maximum=maximum)


break_a_sequence_v1_1(sequence="abc" * 3, tem="", minimum=2, maximum=3)


# version 1.2
# recursion, print results as lists
def break_a_sequence_v1_2(sequence, tem, minimum, maximum):
    if minimum <= len(sequence) <= maximum:
        print tem + [sequence]
    elif len(sequence) == 0:
        print tem

    for i in range(minimum, maximum + 1, 1):
        for j in range(-minimum, -maximum - 1, -1):
            subSeqs = [sequence[:i], sequence[j:]]
            remain = sequence[i:j]
            if len(remain) >= minimum or i + abs(j) == len(sequence):
                break_a_sequence_v1_2(sequence=remain, tem=tem + subSeqs, minimum=minimum, maximum=maximum)


break_a_sequence_v1_2(sequence="abc" * 3, tem=[], minimum=2, maximum=3)


# version 2.0
# recursion, return all results as a list
def break_a_sequence_v2(seq_range, minimum, maximum, tem, results):
    remain_length = seq_range[1] + 1 - seq_range[0]
    if minimum <= remain_length <= maximum:
        results.append(tem + [seq_range])
    elif remain_length == 0:
        results.append(tem)
    for i in range(minimum, maximum + 1, 1):
        for j in range(-minimum, -maximum - 1, -1):
            subSeqs = [[seq_range[0], seq_range[0] + i - 1],
                       [seq_range[1] + j + 1, seq_range[1]]]
            remain = [seq_range[0] + i, seq_range[1] + j]
            remain_length = remain[1] + 1 - remain[0]
            if remain_length >= minimum or remain_length == 0:
                break_a_sequence_v2(seq_range=remain, minimum=minimum, maximum=maximum, tem=tem + subSeqs,
                                    results=results)
    return results


break_a_sequence_v2(seq_range=[0, 8], minimum=2, maximum=3, tem=[], results=[])


# version 3.0
# recursion, backtrack, return when one result has been found, print the result as list
def backtrack(seq_range, minimum, maximum, tem):
    remain_length = seq_range[1] + 1 - seq_range[0]
    if minimum <= remain_length <= maximum or remain_length == 0:  # meet the termination condition
        if minimum <= remain_length <= maximum:
            print tem + [seq_range]
        elif remain_length == 0:
            print tem
        return True
    for i in range(minimum, maximum + 1, 1):
        for j in range(-minimum, -maximum - 1, -1):
            if not isValid(seq_range, minimum, i, j):
                continue
            subs = [[seq_range[0], seq_range[0] + i - 1],
                    [seq_range[1] + j + 1, seq_range[1]]]
            seq_range = [seq_range[0] + i, seq_range[1] + j]  # do
            if backtrack(seq_range, minimum, maximum, tem=tem + subs):
                return True  # a result has been found
            seq_range = [seq_range[0] - i, seq_range[1] - j]  # undo
    return False


def isValid(seq_range, minimum, i, j):
    remain = [seq_range[0] + i, seq_range[1] + j]
    remain_length = remain[1] + 1 - remain[0]
    if remain_length >= minimum or remain_length == 0:
        return True
    else:
        return False


backtrack(seq_range=[0, 8], minimum=2, maximum=3, tem=[])
