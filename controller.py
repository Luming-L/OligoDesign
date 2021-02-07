from iREase import IREase
from primaryFragment import PrimaryFragment
from secondaryFragment import SecondaryFragment
from wrap import Wrap
from vector import Vector
from config import configs


# custom Error class
class DesignError(Exception):
    pass


# some functions
def generate_oligos_from_primaryFragment(primaryFragment):
    """ generate oligos from primaryFragment """

    while True:
        try:
            sequence_length = len(
                primaryFragment.wrap.wrap5 + primaryFragment.original_sequence + primaryFragment.wrap.wrap3)

            # prompt
            print "The sequence is " + str(sequence_length) + \
                  " bases, please input the subSequence group size, the minimum and maximum length of subSequences."
            print "Three parameters above should meet these requirements: "
            print "1. " + \
                  str(oligo_group_size_range[0]) + " <= size <= " + str(oligo_group_size_range[1])
            print "2. " + \
                  str(subSequence_length_range[0]) + " <= minimum <= sequence_length/size <= maximum <= " + str(
                subSequence_length_range[1])

            # input size
            size_input = input("\nsize: ")
            if type(size_input) == int and oligo_group_size_range[0] <= size_input <= oligo_group_size_range[1]:
                if subSequence_length_range[0] <= sequence_length / size_input <= subSequence_length_range[1]:
                    print "sequence_length/size is " + str(sequence_length / size_input) + "."
                else:
                    raise DesignError
            else:
                raise DesignError

            # input minimum
            minimum_input = input("minimum: ")
            if type(minimum_input) == int and subSequence_length_range[
                0] <= minimum_input <= sequence_length / size_input:
                pass
            else:
                raise DesignError

            # input maximum
            maximum_input = input("maximum: ")
            if type(maximum_input) == int and sequence_length / size_input <= maximum_input <= subSequence_length_range[
                1]:
                pass
            else:
                raise DesignError

            # prompt
            if maximum_input - minimum_input > 10:
                while True:
                    selection = raw_input("maximum - minimum > 10, it may take a long time to break the sequence.\n"
                                          "Do you want to continue? (yes/no)\n")
                    if selection == "yes":
                        break
                    if selection == "no":
                        raise DesignError

            primaryFragment.generate_oligos(minimum=minimum_input, maximum=maximum_input, size=size_input)

            if primaryFragment.oligo_group is None:
                raise DesignError

        except (DesignError, NameError, SyntaxError) as e:
            print e
            print "Please try again!\n\n\n"

        else:
            break


# configurations
secondaryFragment_length_range = configs["secondaryFragment_length_range"]
primaryFragment_length_range = configs["primaryFragment_length_range"]
subSequence_length_range = configs["subSequence_length_range"]
reverse_complementary_length_range = configs["reverse_complementary_length_range"]
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

with open("input.txt", "r") as f:
    sf_seq = f.read()


# for pf in range(1, len(secondaryFragment_1.primaryFragment_group) + 1):
#     secondaryFragment_1.primaryFragment_group[pf].generate_oligos(37, 45)
#     print len(secondaryFragment_1.primaryFragment_group[pf].original_sequence)
# for pf in range(1, len(secondaryFragment_1.primaryFragment_group) + 1):
#     print "primary fragment " + str(pf)
#     # for index in sorted(secondaryFragment_1.primaryFragment_group[pf].oligo_group.keys()):
#     #     print "oligo " + str(index) + ":" + secondaryFragment_1.primaryFragment_group[pf].oligo_group[index]
#     print secondaryFragment_1.primaryFragment_group[pf].oligo_group


# primaryFragment_1 = PrimaryFragment(original_sequence=x, wrap=BsaI, iREase=BtsI, vector=vector1)
# print len(primaryFragment_1.sequence)
# primaryFragment_1.generate_oligos(minimum=45, maximum=55)


