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

"""
Input:
"abcabcabc"
Output:
(ab)(ca)(bc)(abc)
(ab)(ca)(bca)(bc)
(ab)(cab)(ca)(bc)
(abc)(ab)(ca)(bc)
(abc)(abc)(abc)
"""

def f(sequence="abc"*3, tem="", min=2, max=3):
    if (len(sequence) >= min and len(sequence) <= max) or len(sequence) == 0:
        print tem + "(" + sequence + ")"
    for i in range(min, max+1, 1):
        for j in range(-min, -max-1, -1):
            oligo1 = sequence[:i]
            oligo2 = sequence[j:]
            remain = sequence[i:j]
            x = tem + "(" + oligo1 + ")" + "(" + oligo2 + ")"
            if len(remain) >= min or i + abs(j) == len(sequence):
                f(sequence=remain, tem=x, min=min, max=max)

"AGGTCTCTCTCTCTCTGGTACCAATCTAGAGGATCCCTCGAGGATTATGTGGAAAAAAAGCACCGACTCGGTGCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTCAACTTGCTATGCTGTTTCCAGCATAGCTCTGAAACTGAGGCAGGCGGGGATGAAGTGCCACGGATCATCTGCACAACTCTTTTAAATCAGCTTTGATCTATGTGGATAGCCGAGGTAGAGACC"