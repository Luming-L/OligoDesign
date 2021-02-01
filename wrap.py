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

    def find_recognition_site(self, sequence):
        """judge whether the sequence has sites that will be cleaved by the restriction enzyme"""
        found = False

        sites = [self.recognition_site,  # the recognition site
                 str(Seq(self.recognition_site).reverse_complement())]  # reverse complement sequence of the site

        for site in sites:
            if sequence.find(site) != -1:
                found = True
                break

        return found


if __name__ == '__main__':
    sequence = "ATCGGGTCTCATCGGAGACG"

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

    print "Does the sequence " + "'" + str(sequence) + "'" + " contain recognition site?"
    print "BsaI: ", BsaI.find_recognition_site(sequence)
    print "BsmBI: ", BsmBI.find_recognition_site(sequence)
