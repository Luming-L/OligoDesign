class Vector:
    """ The primary fragment will be ligated to the vector.
        Vector is the initial double-strand sequence in outside-to-inside mode"""

    def __init__(self, name, sticky_end_length, sticky_end1, sticky_end2):
        self.name = name
        self.sticky_end_length = sticky_end_length
        self.sticky_end1 = sticky_end1
        self.sticky_end2 = sticky_end2


if __name__ == '__main__':
    vector1 = Vector(name="vector1", sticky_end_length=3, sticky_end1="GGT", sticky_end2="AGG")
    print "The length of each sticky end in " + vector1.name + " is " + str(vector1.sticky_end_length) + "."
    print "The sticky ends of " + vector1.name + " are " + vector1.sticky_end1 + " and " + vector1.sticky_end2 + "."
