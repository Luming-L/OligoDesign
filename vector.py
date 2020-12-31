# -*- coding: UTF-8 -*-
__metaclass__ = type


class Vector:
    """Vector, the initial double-strand sequence"""

    def __init__(self, name, sticky_end_length):
        self.name = name
        self.sticky_end_length = sticky_end_length


if __name__ == '__main__':
    vector1 = Vector("vector1", 3)
    print vector1.name
    print vector1.sticky_end_length
