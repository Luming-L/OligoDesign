from iREase import IREase
from primaryFragment import PrimaryFragment
from secondaryFragment import SecondaryFragment
from wrap import Wrap
from vector import Vector
from config import configs
import logging


# custom Error class
class DesignError(Exception):
    pass


# some functions
def generate_primaryFragments_for_secondaryFragment(secondaryFragment, group_size_range, subSequence_length_range):
    """ generate a group of primaryFragments for secondaryFragment """
    for group_size in range(group_size_range[0], group_size_range[1] + 1):
        subSequence_length = int(len(secondaryFragment.original_sequence) / group_size)
        step = 1
        while subSequence_length_range[0] + step <= subSequence_length <= subSequence_length_range[1] - step:
            secondaryFragment.generate_primaryFragments(minimum=subSequence_length - step,
                                                        maximum=subSequence_length + step, size=group_size)
            if secondaryFragment.primaryFragment_group is None:
                step = step + 1
            else:
                return True
    return False


def generate_oligos_for_primaryFragment(primaryFragment, group_size_range, subSequence_length_range):
    """ generate a group of oligos for primaryFragment """
    for group_size in range(group_size_range[0], group_size_range[1] + 1):
        subSequence_length = int(len(primaryFragment.original_sequence) / group_size)
        step = 1
        while subSequence_length_range[0] + step <= subSequence_length <= subSequence_length_range[1] - step:
            primaryFragment.generate_oligos(minimum=subSequence_length - step, maximum=subSequence_length + step,
                                            size=group_size)
            if primaryFragment.oligo_group is None:
                step = step + 1
            else:
                return True
    return False


logging.basicConfig(level=logging.INFO, filename='oligoDesign.log')
logging.info('STARTING PROGRAM')

# configurations
oligo_length_range = configs["oligo_length_range"]
subSequence_of_primaryFragment_length_range = configs["subSequence_of_primaryFragment_length_range"]
primaryFragment_length_range = configs["primaryFragment_length_range"]
secondaryFragment_length_range = configs["secondaryFragment_length_range"]
reverse_complementary_length_range = configs["reverse_complementary_length_range"]
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

sf_seq = raw_input("\nPlease input a sequence: \n")

sf = SecondaryFragment(original_sequence=sf_seq, vectors=all_vectors, wraps=all_wraps, iREases=all_iREases)

logging.info(
    '\n############################## Step 1: Generate Primary Fragments for a Secondary Fragment ##############################')
generate_primaryFragments_for_secondaryFragment(sf, primaryFragment_group_size_range, primaryFragment_length_range)

logging.info(
    '\n############################## Step 2: Generate oligos for each Primary Fragment ##############################')
for pf in sf.primaryFragment_group.values():
    generate_oligos_for_primaryFragment(pf, oligo_group_size_range, subSequence_of_primaryFragment_length_range)

for pf in sf.primaryFragment_group.values():
    print pf.oligo_group

logging.info('\nENDING PROGRAM')



