# -*- coding: UTF-8 -*-
__metaclass__ = type


class PrimaryFragment:
    """primary fragment, can be synthesized by a group of oligos"""

    def __init__(self, sequence, subSequences_groups, oligoGroup):
        self.sequence = sequence
        # all possibilities of breaking the sequence given the length range of sub sequences
        # all possibilities are stored in a list, each possibility is a dict
        self.subSequences_groups = subSequences_groups
        self.oligoGroup = oligoGroup

    def generate_subSequences(self, minimum, maximum):
        def break_a_string():
            return

        def format_results():
            return

        self.subSequences_groups = results