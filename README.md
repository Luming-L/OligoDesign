# Design oligos
This script aims to design oligos for gene synthesis. The gene synthesis method can be found in patent CN104212791A. You can find the latest version [here](https://github.com/Luming-L/OligoDesign).
## Requirements
- python 2.7
- biopython
- re
## Usage
- revise parameters in config.py
- `python OligoDesign/controller.py`
## Definitions
oligo: a fragment can be synthesized by machine directly. Each oligo used in this method contains a reverse complement fragment, a restriction site, and a subSequence of the PF. By the reverse complement fragment, the oligo can form hairpin structure.

primary fragment (PF): a fragment can be synthesized by a group of oligos. A PF has its vector, wrap and iREase.

secondary fragment (SF): a fragment can be synthesized by a group of primary fragments.

vector: the PF will be ligated to the vector.

wrap: including wrap5 and wrap3, is used to ligate a primary fragment (PF) to a vector and ligate several PFs. wrap5 and wrap3 are short sequences containing recognition sites; And the first several bases of wrap5, and the last several bases of wrap3 are complementary to sticky ends of the vector. wrap5 and wrap3 will be added to the head and the tail of a PF. With them, a PF can be ligated to a vector for an amplification. After digestion, several PFs can be ligated by their sticky ends.
 
iREase: restriction enzyme, is used to ligate oligos.
## Workflow
- input sequence
- generate PFs
- generate oligos
- output oligos
### generate PFs
- 1.break the sequence of SF
- 2.determine wrap for each subSequence
- 3.determine iREase for each subSequence
- 4.initiate each PF object with its subSequence, vector, wrap and iREase
### generate oligos
- 1.add wrap sequences to the PF
- 2.break the new sequence
- 3.construct valid hairpins (oligos)
- 3.1 make subSequences overlapped
- 3.2 the first half of subSequences comes from the reverse complement sequence of PF (negative chain), the second comes from the sequence of PF (positive chain)
- 3.3 calculate the length of reverse complement fragment for each subSequence
- 3.4 ligate reverse complement fragments, reconstruction sites and subSequences
- 3.5 judge whether hairpins are valid
- 3.5 1 each only has one recognition site
- 3.5.2 sticky ends of hairpins should be unique
- 3.5.3 sticky ends of intermediate products should be unique
- 3.5.4 all intermediate products should not be ligated (the final product can be ligated)
## Pseudocode
The core process in both "generate PFs" and "generate oligos" is breaking sequence. We use backtrack algorithm to do it. We also create a generator to store all possibilities of PFs. If oligos of a PF cannot be generated, we will `next(generator)` to get another group of PFs and generate oligos again.

Ref: [break a string](https://www.geeksforgeeks.org/print-ways-break-string-bracket-form/)
### create PF generator
```
def create_primaryFragments_generator(seq, subSeq length range):
    generator(remain_seq, subSeq length range)
def generator(remain_seq, subSeq length range):
    if sequence breaking was finished:
        get_valid_primaryFragments(subSeqs)
        if PFs are valid:
            yield PFs
            
    for i in subSeq length range:
        for j in subSeq length range:
            if not isValid(seq, i, j):
                continue
            cut i bases from the head of the sequence
            cut j bases from the tail of the sequence
            make the remaining sequence to be the new sequence
            for result in generator(remain_seq, subSeq length range):
                yield result
def isValid(seq, i, j):
    judge whether the remaining sequence length is not less than minimum subSequence length
    
def get_valid_primaryFragments(subSeqs): 
    determine the wrap and the iREase for each PF sequence
```
### generate oligos
```
def generate_oligos(seq, subSeq length range):
    backtrack(remain_seq, subSeq length range)
    
def backtrack(remain_seq, subSeq length range):
    if sequence breaking was finished:
        if hairpins_are_valid(subSeqs):
            return True
            
    for i in subSeq length range:
        for j in subSeq length range:
            if not isValid(seq, i, j):
                continue
            cut i bases from the head of the sequence
            cut j bases from the tail of the sequence
            make the remaining sequence to be the new sequence
            if backtrack(remain_seq, subSeq length range):
                return True
            cancel cutting 
            
    return False

def isValid(seq, i, j):
    judge whether the remaining sequence length is not less than minimum subSequence length
    
def hairpins_are_valid(subSeqs): 
    construct_hairpins
    judge whether all hairpins are valid
```
## To be continued
deltaG, secondary structure 7-12 1.5kj
