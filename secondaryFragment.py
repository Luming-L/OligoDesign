from primaryFragment import PrimaryFragment


class SecondaryFragment:
    """
    secondary fragment, can be synthesized by a group of primary fragments
    length: 200 nt - 1500 nt
    """

    def __init__(self, original_sequence, iREase_set, wrap_set, vector):
        self.original_sequence = original_sequence.replace(" ", "")
        self.sequence = self.original_sequence
        self.iREase_set = iREase_set
        self.wrap_set = wrap_set
        self.vector = vector
        self.subSeq_group = {}
        self.iREase_group = {}
        self.wrap_group = {}
        self.primaryFragment_group = {}

    def generate_primaryFragments(self, minimum, maximum):

        def backtrack(seq_range, tem):
            remain_length = seq_range[1] + 1 - seq_range[0]
            if minimum <= remain_length <= maximum or remain_length == 0:  # meet the termination condition
                if minimum <= remain_length <= maximum:
                    unformatted_subSeq_group = tem + [seq_range]
                elif remain_length == 0:
                    unformatted_subSeq_group = tem
                print unformatted_subSeq_group
                if primaryFragments_are_valid(unformatted_subSeq_group):
                    return True
            for i in range(minimum, maximum + 1, 1):
                for j in range(-minimum, -maximum - 1, -1):
                    if not isValid(seq_range, i, j):
                        continue
                    subs = [[seq_range[0], seq_range[0] + i - 1],
                            [seq_range[1] + j + 1, seq_range[1]]]
                    seq_range = [seq_range[0] + i, seq_range[1] + j]  # do
                    if backtrack(seq_range, tem=tem + subs):
                        return True  # a result has been found
                    seq_range = [seq_range[0] - i, seq_range[1] - j]  # undo
            return False

        def isValid(seq_range, i, j):
            remain = [seq_range[0] + i, seq_range[1] + j]
            remain_length = remain[1] + 1 - remain[0]
            if remain_length >= minimum or remain_length == 0:
                return True
            else:
                return False

        def primaryFragments_are_valid(unformatted_subSeq_group):

            def format_the_result(a_list):
                a_dict = {}
                for subSeq in range(len(a_list)):
                    a_dict[a_list[subSeq][0]] = a_list[subSeq]
                for new_key in range(1, len(a_dict) + 1):
                    a_dict[new_key] = a_dict.pop(sorted(a_dict.keys())[new_key - 1])
                return a_dict

            def determine_iREase(primaryFragment):
                iREase_usable = None
                for iREase in self.iREase_set:
                    if not iREase.find_recognition_site(primaryFragment):
                        iREase_usable = iREase
                        break
                return iREase_usable

            def determine_wrap(primaryFragment):
                wrap_usable = None
                for wrap in self.wrap_set:
                    if not wrap.find_recognition_site(primaryFragment):
                        wrap_usable = wrap
                        break
                return wrap_usable

            def check_and_instantiate_primaryFragments():
                if None not in self.iREase_group.values() and None not in self.wrap_group.values():
                    for index in range(1, len(self.subSeq_group) + 1):
                        container2 = self.sequence[self.subSeq_group[key][0]:self.subSeq_group[key][1] + 1]
                        self.primaryFragment_group[index] = PrimaryFragment(original_sequence=container2,
                                                                            wrap=self.wrap_group[index],
                                                                            iREase=self.iREase_group[index],
                                                                            vector=self.vector)
                    return True
                else:
                    print "2"
                    # self.subSeq_group = {}
                    # self.iREase_group = {}
                    # self.wrap_group = {}
                    return False

            self.subSeq_group = format_the_result(unformatted_subSeq_group)
            print self.subSeq_group

            for key in range(1, len(self.subSeq_group) + 1):
                container1 = self.sequence[self.subSeq_group[key][0]:self.subSeq_group[key][1] + 1]
                self.iREase_group[key] = determine_iREase(container1)
                self.wrap_group[key] = determine_wrap(container1)
            return check_and_instantiate_primaryFragments()
        backtrack([0, len(self.sequence) - 1], tem=[])
