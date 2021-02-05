from Bio.Seq import Seq
from iREase import IREase
from vector import Vector
from wrap import Wrap
from config import configs


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
        size = size

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

                def clean():
                    self.division_point = None
                    self.rc_fragments_length = None
                    self.subSeq_group = None
                    self.oligo_group = None

                # all oligos only have one recognition site
                for hairpin in hairpins.values():
                    if self.iREase.find_recognition_site(hairpin) != 1:
                        print "More than one recognition site in the hairpin."
                        clean()
                        return False

                # sticky ends of hairpins should be unique
                all_sticky_ends_of_hairpins = []
                for key in range(1, len(hairpins) + 1):
                    if key == 1 or key == len(hairpins):  # the first and the last hairpins
                        all_sticky_ends_of_hairpins.append(hairpins[key][-self.vector.sticky_end_length:])
                    else:
                        all_sticky_ends_of_hairpins.append(hairpins[key][-self.iREase.cleave_location:])
                print "all_sticky_ends_of_hairpins", all_sticky_ends_of_hairpins
                if len(set(all_sticky_ends_of_hairpins)) != len(all_sticky_ends_of_hairpins):
                    clean()
                    return False

                # sticky ends of intermediate products should be unique
                all_sticky_ends_of_inter_products = []
                for key in range(1, len(hairpins) + 1):
                    seoip_start = hairpins[key].find(self.iREase.recognition_site) + len(self.iREase.recognition_site)
                    # sticky ends of intermediate products
                    all_sticky_ends_of_inter_products.append(
                        str(Seq(
                            hairpins[key][seoip_start:seoip_start + self.iREase.cleave_location]).reverse_complement()))
                all_sticky_ends_of_inter_products.append(
                    str(Seq(hairpins[1][-self.vector.sticky_end_length:]).reverse_complement()))
                all_sticky_ends_of_inter_products.append(
                    str(Seq(hairpins[len(hairpins)][-self.vector.sticky_end_length:]).reverse_complement()))
                print "all_sticky_ends_of_inter_products", all_sticky_ends_of_inter_products
                if len(set(all_sticky_ends_of_inter_products)) != len(all_sticky_ends_of_inter_products):
                    clean()
                    return False

                # all intermediate products should not be ligated
                right_sticky_ends_of_inter_products = []
                left_sticky_ends_of_inter_products = []
                for key in range(1, self.division_point + 1):
                    rseoip_start = hairpins[key].find(self.iREase.recognition_site) + len(self.iREase.recognition_site)
                    right_sticky_ends_of_inter_products.append(str(Seq(
                        hairpins[key][rseoip_start:rseoip_start + self.iREase.cleave_location]).reverse_complement()))
                for key in range(self.division_point + 1, len(hairpins) + 1):
                    lseoip_start = hairpins[key].find(self.iREase.recognition_site) + len(self.iREase.recognition_site)
                    left_sticky_ends_of_inter_products.append(str(Seq(
                        hairpins[key][lseoip_start:lseoip_start + self.iREase.cleave_location]).reverse_complement()))
                right_sticky_ends_of_inter_products.append(
                    str(Seq(hairpins[1][-self.vector.sticky_end_length:]).reverse_complement()))
                left_sticky_ends_of_inter_products.append(
                    str(Seq(hairpins[len(hairpins)][-self.vector.sticky_end_length:]).reverse_complement()))
                print "right_sticky_ends_of_inter_products", right_sticky_ends_of_inter_products
                print "left_sticky_ends_of_inter_products", left_sticky_ends_of_inter_products

                for rseoip in right_sticky_ends_of_inter_products:
                    if str(Seq(rseoip).reverse_complement()) in left_sticky_ends_of_inter_products:
                        if not right_sticky_ends_of_inter_products.index(
                                rseoip) == self.division_point - 1 and left_sticky_ends_of_inter_products.index(
                            Seq(rseoip).reverse_complement()) == 0:
                            clean()
                            return False

                return True

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
    # configurations
    primaryFragment_length_range = configs["primaryFragment_length_range"]
    subSequence_length_range = configs["subSequence_length_range"]
    reverse_complementary_length_range = configs["reverse_complementary_length_range"]
    oligo_group_size_range = configs["oligo_group_size_range"]

    # the original sequence of primary fragment
    pf_seq = "AGGTCTCTAGGTGGTACTAATACTAGTCTTTGTTGTCGTCCAATTGCGTAATGGGCCGGCCCATACTGCAATACATGTCCTGAAAGGCTTCATGGCCCACTACGAAATGCTTTTCTCCTACAGTTTATCTTACTTCTTCACATCACGTGGTTTCCGACGTACCCAGTGTTCCCGGCTTCCAGCATTTGCTGGTAGCACCAGTAGAAGACGCCTGTCTTGTGCTATGGAGAGACC"

    # wrap
    BsaI = Wrap(name="BsaI", wrap5="AGGTCTCT", wrap3="AGAGACC", recognition_site="GGTCTC", cleave_location=4)

    # iREase
    BtsI = IREase(name="BtsI", recognition_site="GCAGTG", cleave_location=2)

    # vector
    vector1 = Vector(name="vector1", sticky_end_length=3, sticky_end1="GGT", sticky_end2="AGG")

    # create a PF object
    primaryFragment = PrimaryFragment(original_sequence=pf_seq, wrap=BsaI, iREase=BtsI, vector=vector1)

    # generate oligos for a PF
    while True:
        try:
            class DesignError(Exception):
                pass


            sequence_length1 = len(
                primaryFragment.wrap.wrap5 + primaryFragment.original_sequence + primaryFragment.wrap.wrap3)
            print "The sequence is " + str(sequence_length1) + \
                  " bases, please input the subSequence group size, the minimum and maximum length of subSequences."
            print "Three parameters above should meet these requirements: "
            print "1. " + \
                  str(oligo_group_size_range[0]) + " <= size <= " + str(oligo_group_size_range[1])
            print "2. " + \
                  str(subSequence_length_range[0]) + " <= minimum <= sequence_length/size <= maximum <= " + str(
                subSequence_length_range[1])

            size_input = input("\nsize: ")
            if type(size_input) == int and oligo_group_size_range[0] <= size_input <= oligo_group_size_range[1]:
                if subSequence_length_range[0] <= sequence_length1 / size_input <= subSequence_length_range[1]:
                    print "sequence_length/size is " + str(sequence_length1 / size_input) + "."
                else:
                    raise DesignError
            else:
                raise DesignError

            minimum_input = input("minimum: ")
            if type(minimum_input) == int and subSequence_length_range[
                0] <= minimum_input <= sequence_length1 / size_input:
                pass
            else:
                raise DesignError

            maximum_input = input("maximum: ")
            if type(maximum_input) == int and sequence_length1 / size_input <= maximum_input <= \
                    subSequence_length_range[1]:
                pass
            else:
                raise DesignError

            if maximum_input - minimum_input > 10:
                while True:
                    selection = raw_input("maximum - minimum > 10, it may take a long time to break the sequence.\n"
                                          "Do you want to continue? (yes/no)\n")
                    if selection == "yes":
                        break
                    if selection == "no":
                        raise DesignError

            primaryFragment.generate_oligos(minimum=minimum_input, maximum=maximum_input, size=size_input)
            if primaryFragment.oligo_group is None:
                raise DesignError

        except (DesignError, NameError, SyntaxError):
            print "Invalid number, please try again!\n\n\n"
        else:
            break

    #

    #
    #
    #
    #     except (DesignError, NameError, SyntaxError, ValueError):
    #
