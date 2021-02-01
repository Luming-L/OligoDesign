from Bio.Seq import Seq


class IREase:
    """ iREase, restriction enzyme, is used to ligate oligos. """

    def __init__(self, name, recognition_site, cleave_location):
        self.name = name
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
    sequence = "GCAGTGTTTGTGACG"

    BtsI = IREase(name="BtsI", recognition_site="GCAGTG", cleave_location=2)
    BsrDI = IREase(name="BsrDI", recognition_site="GCAATG", cleave_location=2)

    print "iREase BtsI:" + \
          "\nrecognition site " + BtsI.recognition_site + " " + str(Seq(BtsI.recognition_site).reverse_complement()) + "\n"

    print "iREase BsrDI:" + \
          "\nrecognition site " + BsrDI.recognition_site + " " + str(Seq(BsrDI.recognition_site).reverse_complement()) + "\n"

    print "Does the sequence " + "'" + str(sequence) + "'" + " contain recognition site?"
    print "BtsI: ", BtsI.find_recognition_site(sequence)
    print "BsrDI: ", BsrDI.find_recognition_site(sequence)
