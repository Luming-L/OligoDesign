import logging
from Bio.Seq import Seq
from iREase import IREase
from vector import Vector
from wrap import Wrap
from config import configs


class PrimaryFragment:
    """
    primary fragment, can be synthesized by a group of oligos
    """

    def __init__(self, original_sequence, wrap, iREase, vector):
        self.original_sequence = original_sequence
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
                if minimum <= remain_length <= maximum:  # store the result
                    unformatted_subSeq_group = tem + [seq_range]
                elif remain_length == 0:  # store the result
                    unformatted_subSeq_group = tem
                if len(unformatted_subSeq_group) == size:  # certain group size
                    if hairpins_are_valid(unformatted_subSeq_group):
                        logging.info("Got a group of valid hairpins!\n")
                        return True
            for i in range(minimum, maximum + 1, 1):
                for j in range(-minimum, -maximum - 1, -1):
                    if not isValid(seq_range, i, j):
                        continue
                    subs = [[seq_range[0], seq_range[0] + i - 1],
                            [seq_range[1] + j + 1, seq_range[1]]]  # get two subSequences
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
                return False

        def hairpins_are_valid(unformatted_subSeq_group):
            """ construct hairpins and judge whether they are all valid """

            def format_the_result(a_list):
                """ use a dict to represent a group of subSequences """
                a_dict = {}
                for i in range(len(a_list)):
                    a_dict[a_list[i][0]] = a_list[i]
                for i in range(1, len(a_dict) + 1):
                    a_dict[i] = a_dict.pop(sorted(a_dict.keys())[i - 1])
                for i in range(1, len(a_dict) + 1):
                    logging.info("subSeq " + str(i) + ": " + str(a_dict[i][1] - a_dict[i][0] + 1) + " bases " + str(a_dict[i]))
                return a_dict

            def construct_hairpins(seq, vector, subSeq_group, iREase):
                """ construct hairpins """
                logging.info("- overlapping subSequences, and then getting the sequence of PF for subSequences...")
                hairpins = {}

                # make subSeqs overlapped
                overlapped_subSeq_group = {len(subSeq_group): subSeq_group[len(subSeq_group)]}
                for i in range(1, len(subSeq_group)):
                    overlapped_subSeq_group[i] = [subSeq_group[i][0], subSeq_group[i][1] + iREase.cleave_location]
                for i in range(1, len(overlapped_subSeq_group) + 1):  # extract sequences in positive chain
                    overlapped_subSeq_group[i] = seq[overlapped_subSeq_group[i][0]:overlapped_subSeq_group[i][1] + 1]
                for i in range(1, len(overlapped_subSeq_group) + 1):
                    logging.info("subSeq " + str(i) + ": " + str(len(overlapped_subSeq_group[i])) + " bases " + str(overlapped_subSeq_group[i]))

                # first half of subSequences are in negative chain, second in positive
                division_point = int(round(float(len(subSeq_group)) / 2))  # calculate division point
                self.division_point = division_point
                logging.info("- getting the reverse complement sequence of PF for the first half of subSequences...")
                for i in range(1, division_point + 1):  # half to negative, i.e.reverse complement
                    overlapped_subSeq_group[i] = str(Seq(overlapped_subSeq_group[i]).reverse_complement())

                for i in range(1, len(overlapped_subSeq_group) + 1):
                    logging.info("subSeq " + str(i) + ": " + str(len(overlapped_subSeq_group[i])) + " bases " + overlapped_subSeq_group[i])

                # calculate the length of reverse complement fragments
                rc_fragments_length = 14

                # construct hairpins: ligate reverse complement fragments, reconstruction sites and subSequences
                logging.info("- constructing hairpins...")
                for i in range(1, len(overlapped_subSeq_group) + 1):
                    if i == 1 or i == len(overlapped_subSeq_group):  # rc fragments of first and last subSequences
                        tem_frag = overlapped_subSeq_group[i][
                                   -vector.sticky_end_length - rc_fragments_length:-vector.sticky_end_length]
                    else:  # rc fragments of other subSequences
                        tem_frag = overlapped_subSeq_group[i][
                                   -iREase.cleave_location - rc_fragments_length:-iREase.cleave_location]

                    rc_fragment = str(Seq(tem_frag).reverse_complement())
                    hairpins[i] = rc_fragment + iREase.recognition_site + overlapped_subSeq_group[i]
                    logging.info("oligo " + str(i) + ": " + str(len(hairpins[i])) + " bases " + rc_fragment + " " + iREase.recognition_site + " " + overlapped_subSeq_group[i])

                return hairpins

            def final_check(hairpins):
                """judge whether hairpins are all valid"""

                def clean():
                    """ clean when hairpins are not valid"""
                    self.division_point = None
                    self.rc_fragments_length = None
                    self.subSeq_group = None
                    self.oligo_group = None

                logging.info("- checking hairpins...")
                # all oligos only have one recognition site
                for hairpin in hairpins.values():
                    if self.iREase.count_recognition_site(hairpin) != 1:
                        logging.info("More than one recognition site in the hairpin. Failed!")
                        clean()
                        return False

                # sticky ends of hairpins should be unique
                all_sticky_ends_of_hairpins = []
                for i in range(1, len(hairpins) + 1):
                    if i == 1 or i == len(hairpins):  # the first and the last hairpins
                        all_sticky_ends_of_hairpins.append(hairpins[i][-self.vector.sticky_end_length:])
                    else:  # other hairpins
                        all_sticky_ends_of_hairpins.append(hairpins[i][-self.iREase.cleave_location:])
                logging.info("all_sticky_ends_of_hairpins: " + str(all_sticky_ends_of_hairpins))
                if len(set(all_sticky_ends_of_hairpins)) != len(all_sticky_ends_of_hairpins):
                    logging.info("sticky ends of hairpins are not unique. Failed!")
                    clean()
                    return False

                # sticky ends of intermediate products should be unique
                all_sticky_ends_of_inter_products = []
                for i in range(1, len(hairpins) + 1):
                    seoip_start = hairpins[i].find(self.iREase.recognition_site) + len(self.iREase.recognition_site)
                    all_sticky_ends_of_inter_products.append(  # sticky ends of intermediate products
                        str(Seq(
                            hairpins[i][seoip_start:seoip_start + self.iREase.cleave_location]).reverse_complement()))
                all_sticky_ends_of_inter_products.append(  # sticky ends of the vector
                    str(Seq(hairpins[1][-self.vector.sticky_end_length:]).reverse_complement()))
                all_sticky_ends_of_inter_products.append(
                    str(Seq(hairpins[len(hairpins)][-self.vector.sticky_end_length:]).reverse_complement()))
                logging.info("all_sticky_ends_of_inter_products: " + str(all_sticky_ends_of_inter_products))
                if len(set(all_sticky_ends_of_inter_products)) != len(all_sticky_ends_of_inter_products):
                    logging.info("sticky ends of intermediate products are not unique. Failed!")
                    clean()
                    return False

                # all intermediate products should not be ligated
                right_sticky_ends_of_inter_products = all_sticky_ends_of_inter_products[:self.division_point]
                left_sticky_ends_of_inter_products = all_sticky_ends_of_inter_products[self.division_point:-2]
                right_sticky_ends_of_inter_products.append(all_sticky_ends_of_inter_products[-2])
                left_sticky_ends_of_inter_products.append(all_sticky_ends_of_inter_products[-1])
                logging.info("right_sticky_ends_of_inter_products: " + str(right_sticky_ends_of_inter_products))
                logging.info("left_sticky_ends_of_inter_products: " + str(left_sticky_ends_of_inter_products))

                for rseoip in right_sticky_ends_of_inter_products:
                    if str(Seq(rseoip).reverse_complement()) in left_sticky_ends_of_inter_products:
                        if not right_sticky_ends_of_inter_products.index(rseoip) == self.division_point - 1 and \
                                left_sticky_ends_of_inter_products.index(str(Seq(rseoip).reverse_complement())) == 0:
                            logging.info("intermediate products can be ligated. Failed!")
                            clean()
                            return False

                return True

            self.subSeq_group = format_the_result(unformatted_subSeq_group)

            self.oligo_group = construct_hairpins(seq=sequence,
                                                  vector=self.vector,
                                                  subSeq_group=self.subSeq_group,
                                                  iREase=self.iREase)

            return final_check(self.oligo_group)

        logging.info("- adding wrap sequences to the PF (" + str(len(self.original_sequence)) + " bases),  and then breaking the new sequence...")
        backtrack(seq_range=[0, len(sequence) - 1], tem=[])


if __name__ == '__main__':
    # the original sequence of primary fragment
    pf_seq = "AGGTCTCTTAGCTCTGAAACCTTCGAGCTCCACCGCTCCATGCCACGGATCATCTGCACAACTCTTTTAAATCAGCTTTGATCTATGTGGATAGCCGAGGTGGTACTAATACTAGTCTTTGTTGTCGTCCAATTGCGTAATGGGCCGGCCCATACTGCAATACATGTCCTGAAAGGCTTCATGGCCCACTACGAAATGCTAGAGAAGAGACC"

    # wrap
    BsaI = Wrap(name="BsaI", wrap5="AGGTCTCT", wrap3="AGAGACC", recognition_site="GGTCTC", cleave_location=4)

    # iREase
    BtsI = IREase(name="BtsI", recognition_site="GCAGTG", cleave_location=2)

    # vector
    vector1 = Vector(name="vector1", sticky_end_length=3, sticky_end1="GGT", sticky_end2="AGG")

    # create a PF object
    primaryFragment1 = PrimaryFragment(original_sequence=pf_seq, wrap=BsaI, iREase=BtsI, vector=vector1)

    # generate oligos for a PF
    primaryFragment1.generate_oligos(minimum=40, maximum=50, size=5)  # consider oligo_group_size_range, subSequence_length_range
