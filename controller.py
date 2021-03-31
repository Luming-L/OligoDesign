from iREase import IREase
from primaryFragment import PrimaryFragment
from secondaryFragment import SecondaryFragment
from wrap import Wrap
from vector import Vector
from config import configs


# some functions
def generate_primaryFragments(secondaryFragment, subSeq_num_range, subSeq_length_range):
    """
        Given a secondaryFragment object, range of subSequence number and range of subSequence length,
        generate a group of primaryFragments.
    """
    for subSeq_num in range(subSeq_num_range[0], subSeq_num_range[1] + 1):
        # calculate the average subSeq length given certain number of subSeq
        subSeq_length = int(len(secondaryFragment.original_sequence) / subSeq_num)
        # calculate minimum and maximum subSeq length when breaking secondaryFragment sequence
        # minimum = subSeq_length - step, maximum=subSeq_length + step
        step = 1
        # requirements: minimum => subSeq_num_range[0], maximum <= subSeq_num_range[1], maximum - minimum <= 10
        while subSeq_length_range[0] + step <= subSeq_length <= subSeq_length_range[1] - step and step < 11:
            secondaryFragment.create_primaryFragments_generator(minimum=subSeq_length - step,
                                                                maximum=subSeq_length + step,
                                                                num=subSeq_num)
            try:
                # get PFs
                secondaryFragment.next_primaryFragments()
            except StopIteration:
                step = step + 1
            else:
                return


def generate_oligos(primaryFragment, subSeq_num_range, subSeq_length_range):
    """
        Given a primaryFragment object, range of subSequence number and range of subSequence length,
        generate a group of oligos.
    """
    for subSeq_num in range(subSeq_num_range[0], subSeq_num_range[1] + 1):
        subSeq_length = int(len(primaryFragment.original_sequence) / subSeq_num)
        step = 1
        while subSeq_length_range[0] + step <= subSeq_length <= subSeq_length_range[1] - step and step < 11:
            primaryFragment.generate_oligos(minimum=subSeq_length - step,
                                            maximum=subSeq_length + step,
                                            num=subSeq_num)
            if primaryFragment.oligos is None:
                step = step + 1
            else:
                return


def print_oligos(primaryFragment):
    iREaseA = primaryFragment.iREase
    vectorA = primaryFragment.vector
    for index in primaryFragment.oligos.keys():
        site1 = primaryFragment.oligos[index].find(iREaseA.recognition_site)
        if index == 1 or index == len(primaryFragment.oligos):
            site2 = -vectorA.sticky_end_length
        else:
            site2 = -iREaseA.cleave_location
        print 'oligo ' + str(index) + ': ' + \
              '\033[0;33m' + primaryFragment.oligos[index][:site1] + '\033[0m' + \
              '\033[0;32m' + primaryFragment.oligos[index][site1:site1 + len(iREaseA.recognition_site)] + '\033[0m' + \
              primaryFragment.oligos[index][site1 + len(iREaseA.recognition_site):site2] + \
              '\033[0;31m' + primaryFragment.oligos[index][site2:] + '\033[0m'


# configurations
oligo_length_range = configs["oligo_length_range"]
primaryFragment_length_range = configs["primaryFragment_length_range"]
secondaryFragment_length_range = configs["secondaryFragment_length_range"]
default_reverse_complement_fragment_length = configs["default_reverse_complement_fragment_length"]
primaryFragment_num_range = configs["primaryFragment_num_range"]
oligo_num_range = configs["oligo_num_range"]

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
    generate_primaryFragments(sf, [1, 1], primaryFragment_length_range)
    if sf.primaryFragments is None:
        print "The oligos of the sequence should be designed manually."
elif secondaryFragment_length_range[0] <= len(input_seq) <= secondaryFragment_length_range[1]:
    generate_primaryFragments(sf, primaryFragment_num_range, primaryFragment_length_range)
    if sf.primaryFragments is None:
        print "The oligos of the sequence should be designed manually."
elif len(input_seq) > secondaryFragment_length_range[1]:
    print "The input sequence is too long."

# generate oligos
while sf.primaryFragments is not None:
    try:
        for pf in sf.primaryFragments.values():
            generate_oligos(pf, oligo_num_range, [15, 50])  # generate oligos for each PF
            if pf.oligos is None:
                sf.next_primaryFragments() # if cannot generate oligos, get a new PF group
                raise Exception  # iterate PFs to generate oligos for each PF again
    except Exception as e:
        if type(e) == StopIteration:
            print "The oligos of the sequence should be designed manually."
            break
    else:
        break

# results
if sf.primaryFragments is not None:
    print "\nResults:"
    for i in sf.primaryFragments.keys():
        print "\nprimary fragment " + str(i)
        print "length: " + str(len(sf.primaryFragments[i].original_sequence))
        print "wrap5: " + sf.primaryFragments[i].wrap.wrap5 + "; wrap3: " + sf.primaryFragments[i].wrap.wrap3
        print "iREase: " + sf.primaryFragments[i].iREase.name
        print "vector sticky end1: " + sf.primaryFragments[i].vector.sticky_end1 + \
              "; vector sticky end2: " + sf.primaryFragments[i].vector.sticky_end2
        print_oligos(sf.primaryFragments[i])



