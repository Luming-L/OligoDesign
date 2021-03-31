configs = {"oligo_length_range": [1, 70],  # the sequence can be synthesized by machine directly
           "primaryFragment_length_range": [71, 300],  # subSequence of secondary fragment
           "secondaryFragment_length_range": [301, 1500],

           "default_reverse_complement_fragment_length": 14,  # "reverse_complementary_length_range": [6, 30]

           "primaryFragment_num_range": [2, 10],
           "oligo_num_range": [2, 10],

           "vectors": [{"name": "vector1", "sticky_end_length": 3, "sticky_end1": "GGT", "sticky_end2": "AGG"}],

           "wraps": [{"name": "BsaI", "wrap5": "AGGTCTCT", "wrap3": "AGAGACC", "recognition_site": "GGTCTC",
                      "cleave_location": 4},
                     {"name": "BsmBI", "wrap5": "AGGCGTCTCT", "wrap3": "AGAGACGACC", "recognition_site": "CGTCTC",
                      "cleave_location": 4}],

           "iREases": [{"name": "BtsI", "recognition_site": "GCAGTG", "cleave_location": 2},
                       {"name": "BsrDI", "recognition_site": "GCAATG", "cleave_location": 2}]

           }
