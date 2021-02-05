configs = {"secondaryFragment_length_range": [100, 1500],
           "primaryFragment_length_range": [50, 300],
           "subSequence_length_range": [15, 50],
           "reverse_complementary_length_range": [6, 30],
           "oligo_group_size_range": [2, 8],

           "vectors": [{"name": "vector1", "sticky_end_length": 3, "sticky_end1": "GGT", "sticky_end2": "AGG"}],

           "wraps": [{"name": "BsaI", "wrap5": "AGGTCTCT", "wrap3": "AGAGACC", "recognition_site": "GGTCTC",
                      "cleave_location": 4},
                     {"name": "BsmBI", "wrap5": "AGGCGTCTCT", "wrap3": "AGAGACGACC", "recognition_site": "CGTCTC",
                      "cleave_location": 4}],

           "iREases": [{"name": "BtsI", "recognition_site": "GCAGTG", "cleave_location": 2},
                       {"name": "BsrDI", "recognition_site": "GCAATG", "cleave_location": 2}]

           }
