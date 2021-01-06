# coding=utf-8
__metaclass__ = type

from Bio.Seq import Seq


class PrimaryFragment:
    """primary fragment, can be synthesized by a group of oligos"""

    def __init__(self, original_sequence, wrap, iREase, vector):
        self.original_sequence = original_sequence.replace(" ", "")
        self.wrap = wrap
        self.sequence = self.wrap.wrap5 + self.original_sequence + self.wrap.wrap3
        self.iREase = iREase
        self.vector = vector
        self.division_point = None
        self.rc_fragments_length = None
        self.subSeq_group = None
        self.oligo_group = None

    def change_division_point(self, new_division_point):
        self.division_point = int(new_division_point)

    def generate_oligos(self, minimum, maximum):

        def backtrack(seq_range, tem):
            remain_length = seq_range[1] + 1 - seq_range[0]
            if minimum <= remain_length <= maximum or remain_length == 0:
                if minimum <= remain_length <= maximum:
                    unformatted_subSeq_group = tem + [seq_range]
                elif remain_length == 0:
                    unformatted_subSeq_group = tem
                if hairpins_are_valid(unformatted_subSeq_group):
                    return True
            for i in range(minimum, maximum + 1, 1):
                for j in range(-minimum, -maximum - 1, -1):
                    if not isValid(seq_range, i, j):
                        continue
                    subs = [[seq_range[0], seq_range[0] + i - 1],
                            [seq_range[1] + j + 1, seq_range[1]]]
                    seq_range = [seq_range[0] + i, seq_range[1] + j]  # do
                    if backtrack(seq_range, tem=tem + subs):
                        return True
                    seq_range = [seq_range[0] - i, seq_range[1] - j]  # undo
            return False

        def isValid(seq_range, i, j):
            remain = [seq_range[0] + i, seq_range[1] + j]
            remain_length = remain[1] + 1 - remain[0]
            if remain_length >= minimum or remain_length == 0:
                return True
            else:
                return False

        def hairpins_are_valid(unformatted_subSeq_group):
            def format_the_result(a_list):
                a_dict = {}
                for subSeq in range(len(a_list)):
                    a_dict[a_list[subSeq][0]] = a_list[subSeq]
                for new_key in range(1, len(a_dict) + 1):
                    a_dict[new_key] = a_dict.pop(sorted(a_dict.keys())[new_key - 1])
                return a_dict

            def construct_hairpins(sequence, vector, subSeq_group, division_point, iREase):
                hairpins = {}

                # make subSeqs overlapped
                overlapped_subSeq_group = {len(subSeq_group): subSeq_group[len(subSeq_group)]}
                for key in range(1, len(subSeq_group)):
                    overlapped_subSeq_group[key] = [subSeq_group[key][0], subSeq_group[key][1] + iREase.cleave_location]

                # first half of subSeqs are in negative chain, second in positive
                if division_point is None:  # calculate division point
                    division_point = int(round(float(len(subSeq_group)) / 2))
                for key in range(1, len(overlapped_subSeq_group) + 1):  # extract sequences in positive chain
                    overlapped_subSeq_group[key] = sequence[overlapped_subSeq_group[key][0]:
                                                            overlapped_subSeq_group[key][1] + 1]
                for key in range(1, division_point + 1):  # half to negative, i.e.reverse complement
                    overlapped_subSeq_group[key] = str(Seq(overlapped_subSeq_group[key]).reverse_complement())

                # calculate the length of rc fragments
                rc_fragments_length = 14

                # construct hairpins
                for key in range(1, len(overlapped_subSeq_group) + 1):
                    if key == 1 or key == len(overlapped_subSeq_group):
                        tem_frag = overlapped_subSeq_group[key][
                                   -vector.sticky_end_length - rc_fragments_length:-vector.sticky_end_length]
                    else:
                        tem_frag = overlapped_subSeq_group[key][
                                   -iREase.cleave_location - rc_fragments_length:-iREase.cleave_location]

                    rc_fragment = str(Seq(tem_frag).reverse_complement())
                    hairpins[key] = rc_fragment + iREase.recognition_site + overlapped_subSeq_group[key]

                return hairpins

            def final_check(hairpins):
                # unique sticky ends
                all_sticky_ends = []
                for key in range(1, len(hairpins) + 1):
                    if key == 1 or key == len(hairpins):
                        all_sticky_ends.append(hairpins[key][-self.vector.sticky_end_length:])
                    else:
                        all_sticky_ends.append(hairpins[key][-self.iREase.cleave_location:])
                print all_sticky_ends
                if len(set(all_sticky_ends)) == len(all_sticky_ends):
                    print self.subSeq_group
                    print self.oligo_group
                    return True
                else:
                    self.subSeq_group = None
                    self.oligo_group = None
                    return False

            self.subSeq_group = format_the_result(unformatted_subSeq_group)

            self.oligo_group = construct_hairpins(sequence=self.sequence,
                                                  vector=self.vector,
                                                  subSeq_group=self.subSeq_group,
                                                  division_point=self.division_point,
                                                  iREase=self.iREase
                                                  )

            return final_check(self.oligo_group)

        backtrack(seq_range=[0, len(self.sequence) - 1], tem=[])

    def print_info(self):
        return

