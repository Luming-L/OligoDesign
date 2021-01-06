# coding=utf-8
__metaclass__ = type

from Bio.Seq import Seq


class IREase:
    """iREase: restriction enzyme"""

    def __init__(self, name, recognition_site, cleave_location):
        """

        :rtype: object
        """
        self.name = name
        self.recognition_site = recognition_site
        self.cleave_location = cleave_location  # cleave at the downstream 2bp of the recognition site

    def find_recognition_site(self, sequence):
        """judge whether the sequence has sites that will be cleaved by the restriction enzyme"""
        found = False

        sites = [self.recognition_site,  # the recognition site
                 self.recognition_site[::-1],  # reverse sequence of the site
                 str(Seq(self.recognition_site).reverse_complement()),  # reverse complement sequence of the site
                 str(Seq(self.recognition_site).reverse_complement())[::-1]]  # complement sequence of the site

        for site in sites:
            if sequence.find(site) != -1:
                found = True
                break

        return found


if __name__ == '__main__':
    BtsI = IREase("BtsI", "GCAGTG", 2)

    test_seq1 = "GCAGTGTTT"
    test_seq2 = "CACTGCTTT"
    test_seq3 = "GTGACGTTT"
    test_seq4 = "CGTCACTTT"
    test_seq5 = "CGTTTTTTT"

    print "find BtsI recognition site in test seq1", BtsI.find_recognition_site(test_seq1)
    print "find BtsI recognition site in test seq2", BtsI.find_recognition_site(test_seq2)
    print "find BtsI recognition site in test seq3", BtsI.find_recognition_site(test_seq3)
    print "find BtsI recognition site in test seq4", BtsI.find_recognition_site(test_seq4)
    print "find BtsI recognition site in test seq5", BtsI.find_recognition_site(test_seq5)
