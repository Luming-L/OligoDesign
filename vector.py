class Vector:
    """ The primary fragment will be ligated to the vector.
        Vector is the initial double-strand sequence in outside-to-inside mode"""

    def __init__(self, name, sticky_end_length):
        self.name = name
        self.sticky_end_length = sticky_end_length


if __name__ == '__main__':
    vector1 = Vector(name="vector1", sticky_end_length=3)
    print "The length of each sticky end in " + vector1.name + " is " + str(vector1.sticky_end_length) + "."
