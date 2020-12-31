# -*- coding: UTF-8 -*-
__metaclass__ = type


class PrimaryFragment:
    """primary fragment, can be synthesized by a group of oligos"""

    def __init__(self, sequence):
        self.sequence = sequence.replace(' ', '')

        # all possibilities of breaking the sequence given the length range of subSequences
        # all possibilities are stored in a list, each possibility is a dict
        self.subSequences_groups = None  # list

        # a feasible oligoGroup
        self.oligoGroup = None  # OligoGroup class

    def generate_subSequences(self, minimum, maximum):
        """
            Given a sequence, minimum and maximum length of subSequences,
            store all possibilities of breaking the sequence in self.subSequences_groups

            :param minimum: int, minimum length of subSequences
            :param maximum: int
        """

        def break_a_sequence(seq_range, tem, result):
            """
            break a sequence

            :param seq_range: list, [first character index, last character index] of the sequence
            :param tem: list, store subSequences in each partition round
            :param result: list, store possibilities of breaking the sequence

            :rtype result: list, store all possibilities of breaking the sequence
            """

            # the length of the remaining sequence
            remain_length = seq_range[1] + 1 - seq_range[0]

            # the sequence has been broken completely
            if minimum <= remain_length <= maximum:
                result.append(tem + [seq_range])
            elif remain_length == 0:
                result.append(tem)

            # the sequence has not been broken completely
            else:
                for i in range(minimum, maximum + 1, 1):  # left partition, partition point iterates from min to max
                    for j in range(-minimum, -maximum - 1,
                                   -1):  # right partition, partition point iterates from min to max

                        # given a remaining sequence (seq_range), a left (i) and a right partition point (j),
                        # generate 2 subSequences from left and right, respectively
                        # the subSequence is represented by [first character index, last character index]
                        subSequences = [[seq_range[0], seq_range[0] + i - 1],
                                        [seq_range[1] + j + 1, seq_range[1]]]

                        remain = [seq_range[0] + i,
                                  seq_range[1] + j]  # the start and end index of the remaining sequence
                        remain_length = remain[1] + 1 - remain[0]  # the length of the remain

                        if remain_length >= minimum or remain_length == 0:  # start the next round of partition
                            break_a_sequence(seq_range=remain, tem=tem + subSequences, results=result)
            return result

        def format_results(results):
            """
            format results, i.e. use a dict to store each possibility rather than a list

            :param results: list, store all possibilities of breaking the sequence, each possibility is a list

            :rtype output: list, store all possibilities of breaking the sequence, each possibility is a dict
            """
            output = []
            for group in range(len(results)):  # iterate all possibilities of breaking the sequence
                tem = {}  # use a dict to store a possibility

                for subSeq in range(
                        len(results[group])):  # iterate all subSequences of a possibility, store them in a dict
                    # in the temporary dict "tem", keys are start indices of subSequences and values are subSequences
                    tem[results[group][subSeq][0]] = results[group][subSeq]

                for new_key in range(1, len(tem) + 1):  # iterate to change keys of "tem"
                    # use the order of subSequence positions as new keys
                    # sort old keys, pop items according to the order, change to new keys
                    tem[new_key] = tem.pop(sorted(tem.keys())[new_key - 1])

                output.append(tem)  # append each possibility to the output
            return output








if __name__ == '__main__':
    x = "AGGTCTCTCTCTCTCTGGTACCAATCTAGAGGATCCCTCGAGGATTATGTGGAAAAAAAGCACCGACTCGGTGCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTCAACTTGCTATGCTGTTTCCAGCATAGCTCTGAAACTGAGGCAGGCGGGGATGAAGTGCCACGGATCATCTGCACAACTCTTTTAAATCAGCTTTGATCTATGTGGATAGCCGAGGTAGAGACC   "

    y = "abcabcabc"

    primaryFragment_1 = PrimaryFragment(x)
    primaryFragment_1.generate_subSequences(minimum=10, maximum=60)
    print len(primaryFragment_1.subSequences_groups)

    primaryFragment_1.add_oligoGroup({})
    print primaryFragment_1.oligoGroup
