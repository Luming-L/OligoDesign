# coding=utf-8
__metaclass__ = type


class Wrap:
    """wrap"""

    def __init__(self, name, wrap5, wrap3, cleave_location):
        self.name = name

        self.wrap5 = wrap5

        self.wrap3 = wrap3

        self.cleave_location = cleave_location


if __name__ == '__main__':
    BsaI = Wrap(name="BsaI", wrap5="AGGTCTCT", wrap3="AGAGACC", cleave_location=4)
    print BsaI.wrap5
    print BsaI.wrap3
    print BsaI.cleave_location
