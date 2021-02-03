from Bio.Seq import Seq
from iREase import IREase
from vector import Vector
from wrap import Wrap


class PrimaryFragment:
    """
    primary fragment, can be synthesized by a group of oligos
    length: 100 bases - 600 bases
    """

    def __init__(self, original_sequence, wrap, iREase, vector):
        self.original_sequence = original_sequence.replace(" ", "")
        self.wrap = wrap
        self.iREase = iREase
        self.vector = vector
        self.division_point = None
        self.rc_fragments_length = None  # the reverse complement fragment to form hairpin structure
        self.subSeq_group = None
        self.oligo_group = None

    def generate_oligos(self, minimum, maximum, size):
        """ generate oligos """
        sequence = self.wrap.wrap5 + self.original_sequence + self.wrap.wrap3

        def backtrack(seq_range, tem):
            """main sub function in generate_oligos"""
            remain_length = seq_range[1] + 1 - seq_range[0]
            if minimum <= remain_length <= maximum or remain_length == 0:  # means all subSequences meet the requirement
                unformatted_subSeq_group = None
                if minimum <= remain_length <= maximum:
                    unformatted_subSeq_group = tem + [seq_range]
                elif remain_length == 0:
                    unformatted_subSeq_group = tem
                print "Sequence breaking has been finished!"
                if hairpins_are_valid(unformatted_subSeq_group):
                    return True
            for i in range(minimum, maximum + 1, 1):
                for j in range(-minimum, -maximum - 1, -1):
                    print "i = " + str(i), "j = " + str(j)
                    if not isValid(seq_range, i, j):
                        continue
                    subs = [[seq_range[0], seq_range[0] + i - 1],
                            [seq_range[1] + j + 1, seq_range[1]]]
                    seq_range = [seq_range[0] + i, seq_range[1] + j]  # cut
                    if backtrack(seq_range, tem=tem + subs):
                        return True
                    seq_range = [seq_range[0] - i, seq_range[1] - j]  # cancel
            return False

        def isValid(seq_range, i, j):
            """ judge whether the remaining sequence length is not less than than minimum """
            remain = [seq_range[0] + i, seq_range[1] + j]
            remain_length = remain[1] + 1 - remain[0]
            if remain_length >= minimum or remain_length == 0:
                return True
            else:
                print "The remaining sequence is " + str(remain_length) + " bases. Not valid, go to next round."
                return False

        def hairpins_are_valid(unformatted_subSeq_group):
            """ construct hairpins and judge whether they are all valid """

            def format_the_result(a_list):
                """ use a dict to represent a group of subSequences """
                a_dict = {}
                for subSeq in range(len(a_list)):
                    a_dict[a_list[subSeq][0]] = a_list[subSeq]
                for new_key in range(1, len(a_dict) + 1):
                    a_dict[new_key] = a_dict.pop(sorted(a_dict.keys())[new_key - 1])
                for num in range(1, len(a_dict) + 1):
                    print "subSeq " + str(num) + ":", str(a_dict[num][1] - a_dict[num][0] + 1) + " bases", a_dict[num]
                return a_dict

            def construct_hairpins(seq, vector, subSeq_group, iREase):
                """ construct hairpins """
                print "construct hairpins"
                hairpins = {}

                # make subSeqs overlapped
                overlapped_subSeq_group = {len(subSeq_group): subSeq_group[len(subSeq_group)]}
                for key in range(1, len(subSeq_group)):
                    overlapped_subSeq_group[key] = [subSeq_group[key][0], subSeq_group[key][1] + iREase.cleave_location]

                print "make subSeqs overlapped"
                for num in range(1, len(overlapped_subSeq_group) + 1):
                    print "subSeq " + str(num) + ":", \
                        str(overlapped_subSeq_group[num][1] - overlapped_subSeq_group[num][0] + 1) + " bases", \
                        overlapped_subSeq_group[num]

                # first half of subSeqs are in negative chain, second in positive
                print "first half of subSeqs are in negative chain, second in positive"
                division_point = int(round(float(len(subSeq_group)) / 2))  # calculate division point
                self.division_point = division_point
                print "The division point is " + str(self.division_point) + "."
                for key in range(1, len(overlapped_subSeq_group) + 1):  # extract sequences in positive chain
                    overlapped_subSeq_group[key] = seq[
                                                   overlapped_subSeq_group[key][0]:overlapped_subSeq_group[key][1] + 1]
                for num in range(1, len(overlapped_subSeq_group) + 1):
                    print "subSeq " + str(num) + ":", \
                        str(len(overlapped_subSeq_group[num])) + " bases", \
                        overlapped_subSeq_group[num]
                print ""
                for key in range(1, division_point + 1):  # half to negative, i.e.reverse complement
                    overlapped_subSeq_group[key] = str(Seq(overlapped_subSeq_group[key]).reverse_complement())

                for num in range(1, len(overlapped_subSeq_group) + 1):
                    print "subSeq " + str(num) + ":", \
                        str(len(overlapped_subSeq_group[num])) + " bases", \
                        overlapped_subSeq_group[num]

                # calculate the length of rc fragments
                rc_fragments_length = 14

                # construct hairpins
                print "hairpins are"
                for key in range(1, len(overlapped_subSeq_group) + 1):
                    if key == 1 or key == len(overlapped_subSeq_group):
                        tem_frag = overlapped_subSeq_group[key][
                                   -vector.sticky_end_length - rc_fragments_length:-vector.sticky_end_length]
                    else:
                        tem_frag = overlapped_subSeq_group[key][
                                   -iREase.cleave_location - rc_fragments_length:-iREase.cleave_location]

                    rc_fragment = str(Seq(tem_frag).reverse_complement())
                    hairpins[key] = rc_fragment + iREase.recognition_site + overlapped_subSeq_group[key]
                    print "oligo " + str(key) + ":", \
                        str(len(hairpins[key])) + " bases", \
                        rc_fragment + " " + iREase.recognition_site + " " + overlapped_subSeq_group[key]

                return hairpins

            def final_check(hairpins):
                """judge whether hairpins are all valid"""
                # unique sticky ends
                all_sticky_ends = []
                for key in range(1, len(hairpins) + 1):
                    if key == 1 or key == len(hairpins):
                        all_sticky_ends.append(hairpins[key][-self.vector.sticky_end_length:])
                    else:
                        all_sticky_ends.append(hairpins[key][-self.iREase.cleave_location:])
                print all_sticky_ends
                if len(set(all_sticky_ends)) == len(all_sticky_ends):
                    return True
                else:
                    self.division_point = None
                    self.rc_fragments_length = None
                    self.subSeq_group = None
                    self.oligo_group = None
                    return False

            # certain group size
            if len(unformatted_subSeq_group) != size:
                return False

            self.subSeq_group = format_the_result(unformatted_subSeq_group)

            self.oligo_group = construct_hairpins(seq=sequence,
                                                  vector=self.vector,
                                                  subSeq_group=self.subSeq_group,
                                                  iREase=self.iREase)

            return final_check(self.oligo_group)

        backtrack(seq_range=[0, len(sequence) - 1], tem=[])


if __name__ == '__main__':
    # the original sequence of primary fragment
    pf_seq = "CTCTGGTACCAATCTAGAGGATCCCTCGAGGATTATGTGGAAAAAAAGCACCGACTCGGTGCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTCAACTTGCTATGCTGTTTCCAGCATAGCTCTGAAACTGAGGCAGGCGGGGATGAAGTGCCACGGATCATCTGCACAACTCTTTTAAATCAGCTTTGATCTATGTGGATAGCCGAGGTGGTACTAATACTAGTCTTTGTTGTCGTCCAATTGCGTA"

    # wrap
    BsaI = Wrap(name="BsaI", wrap5="AGGTCTCT", wrap3="AGAGACC", recognition_site="GGTCTC", cleave_location=4)

    # iREase
    BtsI = IREase(name="BtsI", recognition_site="GCAGTG", cleave_location=2)

    # vector
    vector1 = Vector(name="vector1", sticky_end_length=3, sticky_end1="GGT", sticky_end2="AGG")

    # create a PF object
    primaryFragment = PrimaryFragment(original_sequence=pf_seq, wrap=BsaI, iREase=BtsI, vector=vector1)

    # generate oligos for a PF

    primaryFragment.generate_oligos(minimum=50, maximum=56, size=5)
