# coding=utf-8
__metaclass__ = type

from Bio.Seq import Seq


class Wrap:
    """wrap"""

    def __init__(self, name, wrap5, wrap3, recognition_site, cleave_location):

        self.name = name

        self.wrap5 = wrap5

        self.wrap3 = wrap3

        self.recognition_site = recognition_site

        self.cleave_location = cleave_location

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



