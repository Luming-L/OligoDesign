def oligoDesign(sequence="ATCG", tem="", min=60, max=75, len_sticky=2):
    """
    Given a DNA sequence, generate a group of oligos
    :param sequence: the DNA sequence
    :param tem: oligos that has been generated
    :param min: the max length of oligos
    :param max: the minimum length of oligos
    :param len_sticky: the length of sticky ends
    :return: oligoGroup: a group of oligos
    """
    if sequence == "":
        oligoGroup = tem
        print(oligoGroup)
    for i in range(1, len(sequence) + 1):
        oligo = sequence[:i]
        remain = sequence[i:]
        x = tem + "(" + oligo + ")"
        if len(oligo) >= min and len(oligo) <= max:
            oligoDesign(sequence=remain, tem=x, min=min, max=max, len_sticky=len_sticky)



if __name__ == "__main__":
    # input string
    oligoDesign(sequence="abcde", tem="", min=2, max=3, len_sticky=2)
