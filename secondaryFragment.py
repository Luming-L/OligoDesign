import logging

from config import configs
from iREase import IREase
from primaryFragment import PrimaryFragment
from vector import Vector
from wrap import Wrap


class SecondaryFragment:
    """
    secondary fragment, can be synthesized by a group of primary fragments
    """

    def __init__(self, original_sequence, vectors, wraps, iREases):
        self.original_sequence = original_sequence
        self.vectors = vectors
        self.wraps = wraps
        self.iREases = iREases
        self.primaryFragments_generator = []
        self.primaryFragments = {}

    def create_primaryFragments_generator(self, minimum, maximum):

        def generator(seq_range, container):
            remain_length = seq_range[1] + 1 - seq_range[0]
            if remain_length == 0 or minimum <= remain_length <= maximum:  # the last subSeq
                unformatted_subSeq_group = None
                if remain_length == 0:
                    unformatted_subSeq_group = container
                elif minimum <= remain_length <= maximum:
                    unformatted_subSeq_group = container + [seq_range]
                valid_primaryFragment_group = get_valid_primaryFragments(unformatted_subSeq_group)
                if valid_primaryFragment_group is not None:
                    yield valid_primaryFragment_group
            else:
                for i in range(minimum, maximum + 1, 1):
                    for j in range(-minimum, -maximum - 1, -1):
                        if not isValid(seq_range, i, j):
                            continue
                        subs = [[seq_range[0], seq_range[0] + i - 1], [seq_range[1] + j + 1, seq_range[1]]]
                        for result in generator([seq_range[0] + i, seq_range[1] + j], container=container + subs):
                            yield result

        def isValid(seq_range, i, j):
            """ judge whether the remaining sequence length is not less than than minimum """
            remain = [seq_range[0] + i, seq_range[1] + j]
            remain_length = remain[1] + 1 - remain[0]
            if remain_length >= minimum or remain_length == 0:
                return True
            else:
                return False

        def get_valid_primaryFragments(unformatted_subSeq_group):
            """ get valid primary fragments, each valid pf has a vector, a wrap and a iREase """

            def format_the_result(a_list):
                """ use a dict to represent a group of subSequences """
                a_dict = {}
                for i in range(len(a_list)):
                    a_dict[a_list[i][0]] = a_list[i]
                for i in range(1, len(a_dict) + 1):
                    a_dict[i] = a_dict.pop(sorted(a_dict.keys())[i - 1])
                return a_dict

            def determine_wrap(seq):
                """ determine the wrap of a primary fragment """
                wrap_usable = None
                for a_wrap in self.wraps:
                    if a_wrap.count_recognition_site(seq) == 0:
                        wrap_usable = a_wrap
                        break
                return wrap_usable

            def determine_iREase(seq):
                """ determine the iREase of a primary fragment """
                iREase_usable = None
                for a_iREase in self.iREases:
                    if a_iREase.count_recognition_site(seq) == 0:
                        iREase_usable = a_iREase
                        break
                return iREase_usable

            subSeq_group = format_the_result(unformatted_subSeq_group)
            primaryFragment_group = {}

            # determine the wrap and the iREase for each PF sequence
            # initiate each PF object
            for index in range(1, len(subSeq_group) + 1):
                pf_seq = self.original_sequence[subSeq_group[index][0]:subSeq_group[index][1] + 1]  # sequence of PF
                wrap_of_pf = determine_wrap(pf_seq)  # wrap
                if wrap_of_pf is not None:
                    iREase_of_pf = determine_iREase(wrap_of_pf.wrap5 + pf_seq + wrap_of_pf.wrap3)  # iREase
                    if iREase_of_pf is not None:
                        primaryFragment_group[index] = PrimaryFragment(original_sequence=pf_seq,
                                                                       wrap=wrap_of_pf,
                                                                       iREase=iREase_of_pf,
                                                                       vector=self.vectors[0])
                    else:
                        return None
                else:
                    return None

            return primaryFragment_group

        self.primaryFragments_generator = generator(seq_range=[0, len(self.original_sequence) - 1], container=[])
        self.primaryFragments = {}

    def next_primaryFragments(self):
        self.primaryFragments = next(self.primaryFragments_generator)


if __name__ == '__main__':
    # the original sequence of secondary fragment
    sf_seq = "CTCTGGTACCAATCTAGAGGATCCCTCGAGGATTATGTGGAAAAAAAGCACCGACTCGGTGCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTCAACTTGCTATGCTGTTTCCAGCATAGCTCTGAAACTGAGGCAGGCGGGGATGAAGTGCCACGGATCATCTGCACAACTCTTTTAAATCAGCTTTGATCTATGTGGATAGCCGAGGTGGTACTAATACTAGTCTTTGTTGTCGTCCAATTGCGTAATGGGCCGGCCCATACTGCAATACATGTCCTGAAAGGCTTCATGGCCCACTACGAAATGCTTTTCTCCTACAGTTTATCTTACTTCTTCACATCACGTGGTTTCCGACGTACCCAGTGTTCCCGGCTTCCAGCATTTGCTGGTAGCACCAGTAGAAGACGCCTGTCTTGTGCTATGGTCCCTGACTGCACATCTGATTCCTCCAAGATCCATGCATGCCTGATAACTTTAAGTTGCTTCAGAAGAACTTTAAGTGATCTGTTCGTATGTTTAAAGATTCCTTAATCGTCGAGGATTATGTGGAAAAAAAGCACCGACTCGGTGCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTCAACTTGCTATGCTGTTTCCAGCATAGCTCTGAAACCTTCGAGCTCCACCGCTCCATGCCACGGATCATCTGCACAACTCTTTTAAATCAGCTTTGATCTATGTGGATAGCCGAGGTGGTACTAATACTAGTCTTTGTTGTCGTCCAATTGCGTAATGGGCCGGCCCATACTGCAATACATGTCCTGAAAGGCTTCATGGCCCACTACGAAATGCT"

    # configurations
    primaryFragment_length_range = configs["primaryFragment_length_range"]

    all_vectors = []
    all_wraps = []
    all_iREases = []

    # all vectors
    for vector in configs["vectors"]:  # use the first one by default
        all_vectors.append(Vector(name=vector["name"],
                                  sticky_end_length=vector["sticky_end_length"],
                                  sticky_end1=vector["sticky_end1"],
                                  sticky_end2=vector["sticky_end2"]))

    # all wraps
    for wrap in configs["wraps"]:
        all_wraps.append(Wrap(name=wrap["name"],
                              wrap5=wrap["wrap5"],
                              wrap3=wrap["wrap3"],
                              recognition_site=wrap["recognition_site"],
                              cleave_location=wrap["cleave_location"]))

    # all iREases
    for iREase in configs["iREases"]:
        all_iREases.append(IREase(name=iREase["name"],
                                  recognition_site=iREase["recognition_site"],
                                  cleave_location=iREase["cleave_location"]))

    secondaryFragment1 = SecondaryFragment(original_sequence=sf_seq, vectors=all_vectors, wraps=all_wraps,
                                           iREases=all_iREases)

    secondaryFragment1.create_primaryFragments_generator(minimum=205, maximum=208)

    secondaryFragment1.next_primaryFragments()

    for pf in secondaryFragment1.primaryFragments.values():
        print pf.original_sequence

    secondaryFragment1.next_primaryFragments()

    for pf in secondaryFragment1.primaryFragments.values():
        print pf.original_sequence
