import re
from Bio.Seq import Seq


class IREase:
    """ iREase, restriction enzyme, is used to ligate oligos. """

    def __init__(self, name, recognition_site, cleave_location):
        self.name = name
        self.recognition_site = recognition_site
        self.cleave_location = cleave_location  # the length of the sticky end after digestion

    def find_recognition_site(self, sequence):
        """ how many sites in the sequence that will be cleaved by the restriction enzyme """

        # the sequence and the reverse complement sequence of the site
        regex = "(?=" + self.recognition_site + "|" + str(Seq(self.recognition_site).reverse_complement()) + ")"

        # find all sites, including sites overlapped
        sites = list(re.finditer(regex, sequence))

        # the number of sites found
        site_num = len(sites)

        return site_num


if __name__ == '__main__':
    sequence1 = "CACTGCGCAGTGCAGTGCAATGCATTGC"

    BtsI = IREase(name="BtsI", recognition_site="GCAGTG", cleave_location=2)
    BsrDI = IREase(name="BsrDI", recognition_site="GCAATG", cleave_location=2)

    print "iREase BtsI:" + \
          "\nrecognition site " + BtsI.recognition_site + " " + str(Seq(BtsI.recognition_site).reverse_complement()) + "\n"

    print "iREase BsrDI:" + \
          "\nrecognition site " + BsrDI.recognition_site + " " + str(Seq(BsrDI.recognition_site).reverse_complement()) + "\n"

    print "Does the sequence " + "'" + str(sequence1) + "'" + " contain recognition site?"
    print "BtsI: ", BtsI.find_recognition_site(sequence1)
    print "BsrDI: ", BsrDI.find_recognition_site(sequence1)
