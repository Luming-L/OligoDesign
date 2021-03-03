from iREase import IREase
from primaryFragment import PrimaryFragment
from secondaryFragment import SecondaryFragment
from wrap import Wrap
from vector import Vector
from config import configs


# some functions
def generate_primaryFragments_for_secondaryFragment(secondaryFragment, group_size_range, subSequence_length_range):
    """ generate a group of primaryFragments for secondaryFragment """
    for group_size in range(group_size_range[0], group_size_range[1] + 1):
        subSequence_length = int(
            len(secondaryFragment.original_sequence) / group_size)  # average subSeq length given certain group size
        step = 1
        while subSequence_length_range[0] + step <= subSequence_length <= subSequence_length_range[
            1] - step and step < 11:
            secondaryFragment.create_primaryFragments_generator(minimum=subSequence_length - step,
                                                                maximum=subSequence_length + step,
                                                                num=group_size)
            if secondaryFragment.primaryFragments is None:
                step = step + 1
            else:
                return


def generate_oligos_for_primaryFragment(primaryFragment, group_size_range, subSequence_length_range):
    """ generate a group of oligos for primaryFragment """
    for group_size in range(group_size_range[0], group_size_range[1] + 1):
        subSequence_length = int(len(primaryFragment.original_sequence) / group_size)
        step = 1
        while subSequence_length_range[0] + step <= subSequence_length <= subSequence_length_range[
            1] - step and step < 11:
            primaryFragment.generate_oligos(minimum=subSequence_length - step,
                                            maximum=subSequence_length + step,
                                            num=group_size)
            if primaryFragment.oligos is None:
                step = step + 1
            else:
                return


# configurations
oligo_length_range = configs["oligo_length_range"]
primaryFragment_length_range = configs["primaryFragment_length_range"]
secondaryFragment_length_range = configs["secondaryFragment_length_range"]
default_reverse_complement_fragment_length = configs["default_reverse_complement_fragment_length"]
primaryFragment_group_size_range = configs["primaryFragment_group_size_range"]
oligo_group_size_range = configs["oligo_group_size_range"]

all_vectors = []
all_wraps = []
all_iREases = []

# all vectors
for vector in configs["vectors"]:  # use the first one by default
    all_vectors.append(Vector(name=vector["name"],
                              sticky_end_length=vector["sticky_end_length"],
                              sticky_end1=vector["sticky_end1"],
                              sticky_end2=vector["sticky_end2"]))

# all wraps
for wrap in configs["wraps"]:
    all_wraps.append(Wrap(name=wrap["name"],
                          wrap5=wrap["wrap5"],
                          wrap3=wrap["wrap3"],
                          recognition_site=wrap["recognition_site"],
                          cleave_location=wrap["cleave_location"]))

# all iREases
for iREase in configs["iREases"]:
    all_iREases.append(IREase(name=iREase["name"],
                              recognition_site=iREase["recognition_site"],
                              cleave_location=iREase["cleave_location"]))

# input sequence
input_seq = raw_input("\nPlease input a sequence: \n")

# initialize a secondaryFragment object
sf = SecondaryFragment(original_sequence=input_seq, vectors=all_vectors, wraps=all_wraps, iREases=all_iREases)

# generate Primary Fragments
if oligo_length_range[0] <= len(input_seq) <= oligo_length_range[1]:
    print "The input sequence can be synthesized by machine directly."
elif primaryFragment_length_range[0] <= len(input_seq) <= primaryFragment_length_range[1]:
    generate_primaryFragments_for_secondaryFragment(sf, [1, 1], primaryFragment_length_range)
elif secondaryFragment_length_range[0] <= len(input_seq) <= secondaryFragment_length_range[1]:
    generate_primaryFragments_for_secondaryFragment(sf, primaryFragment_group_size_range, primaryFragment_length_range)
elif len(input_seq) > secondaryFragment_length_range[1]:
    print "The input sequence is too long."

# generate oligos
while True:
    try:
        sf.next_primaryFragments()
        for pf in sf.primaryFragments.values():
            generate_oligos_for_primaryFragment(pf, oligo_group_size_range, [15, 50])
            if pf.oligos is None:
                raise Exception
    except Exception as e:
        if type(e) == StopIteration:
            print "The oligos of the sequence should be designed manually."
            break
    else:
        break

# results
print "\nResults:"
for i in sf.primaryFragments.keys():
    print "\nprimary fragment " + str(i)
    print "length: " + str(len(sf.primaryFragments[i].original_sequence))
    print "wrap: " + sf.primaryFragments[i].wrap.name
    print "iREase: " + sf.primaryFragments[i].iREase.name
    print "vector: " + sf.primaryFragments[i].vector.name
    for j in sf.primaryFragments[i].oligos:
        print "oligo " + str(j) + ": " + sf.primaryFragments[i].oligos[j]

