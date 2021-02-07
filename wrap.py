import re
from Bio.Seq import Seq


class Wrap:
    """ Wrap, with wrap5 and wrap3, is used to ligate a primary fragment (PF) to a vector and ligate several PFs.

        wrap5 and wrap3 are short sequences containing recognition sites.
        The first several bases of wrap5 and the last several bases of wrap3 are complementary to sticky ends of the vector.

        wrap5 and wrap3 will be added to the head and the tail of a PF.
        With them, a PF can be ligated to a vector for amplification.
        After digestion, several PFs can be ligated by their sticky ends."""

    def __init__(self, name, wrap5, wrap3, recognition_site, cleave_location):
        self.name = name
        self.wrap5 = wrap5
        self.wrap3 = wrap3
        self.recognition_site = recognition_site
        self.cleave_location = cleave_location  # the length of the sticky end after digestion

    def count_recognition_site(self, sequence):
        """ how many sites in the sequence that will be cleaved by the restriction enzyme """

        # the sequence and the reverse complement sequence of the site
        regex = "(?=" + self.recognition_site + "|" + str(Seq(self.recognition_site).reverse_complement()) + ")"

        # find all sites, including sites overlapped
        sites = list(re.finditer(regex, sequence))

        # the number of sites found
        site_num = len(sites)

        return site_num


if __name__ == '__main__':
    sequence1 = "ATCGGGTCTCATCGGAGACGGTCTCGAGACCGAGACG"

    BsaI = Wrap(name="BsaI", wrap5="AGGTCTCT", wrap3="AGAGACC", recognition_site="GGTCTC", cleave_location=4)
    BsmBI = Wrap(name="BsmBI", wrap5="AGGCGTCTCT", wrap3="AGAGACGACC", recognition_site="CGTCTC", cleave_location=4)

    print "BsaI wrap:" + \
          "\nwrap5 " + BsaI.wrap5 + \
          "\nwrap3 " + BsaI.wrap3 + \
          "\nrecognition site " + BsaI.recognition_site + " " + str(Seq(BsaI.recognition_site).reverse_complement()) + "\n"

    print "BsmBI wrap:" + \
          "\nwrap5 " + BsmBI.wrap5 + \
          "\nwrap3 " + BsmBI.wrap3 + \
          "\nrecognition site " + BsmBI.recognition_site + " " + str(Seq(BsmBI.recognition_site).reverse_complement()) + "\n"

    print "Does the sequence " + "'" + str(sequence1) + "'" + " contain recognition site?"
    print "BsaI: ", BsaI.count_recognition_site(sequence1)
    print "BsmBI: ", BsmBI.count_recognition_site(sequence1)

