# Design oligos
This script aims to design oligos for gene synthesis. The gene synthesis method can be found in patent CN104212791A. The implementation is based on object-oriented programming.
## Definition
oligo: a fragment can be synthesized by machine directly. Each oligo used in this method contains a reverse complement fragment, a reconstruction site, and a subSequence of the PF. By the reverse complement fragment, the oligo can form hairpin structure.

primary fragment (PF): a fragment can be synthesized by a group of oligos. A PF has its vector, wrap and iREase.

secondary fragment (SF): a fragment can be synthesized by a group of primary fragments.
## Workflow
Input a SF, output oligos.
### generate PFs for a SF
- 1.break the sequence of SF
- 2.determine wrap for each subSequence
- 3.determine iREase for each subSequence
- 4.initiate each PF object with its subSequence, vector, wrap and iREase
### generate oligos for a PF
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
The core process in both "generate PFs for a SF" and "generate oligos for a PF" is breaking sequence. We use backtrack algorithm to do it.

Ref: [break a string](https://www.geeksforgeeks.org/print-ways-break-string-bracket-form/)
### generate PFs
```
def generate_primaryFragments(seq, subSeq length range):
    backtrack(seq)
    
def backtrack(seq):
    if sequence breaking was finished:
        if primaryFragments_are_valid(subSeqs):
            return True
            
    for i in subSeq length range:
        for j in subSeq length range:
            if not isValid(seq, i, j):
                continue
            cut i bases from the head of the sequence
            cut j bases from the tail of the sequence
            make the remaining sequence to be the new sequence
            if backtrack(seq):
                return True
            cancel cutting 
            
    return False

def isValid(seq, i, j):
    judge whether the remaining sequence length is not less than minimum subSequence length
    
def primaryFragments_are_valid(subSeqs): 
    determine the wrap and the iREase for each subSeq
    if each PF has its wrap and iREase, initiate all PF objects
```
### generate oligos
```
def generate_oligos(seq, subSeq length range):
    backtrack(seq)
    
def backtrack(seq):
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
            if backtrack(seq):
                return True
            cancel cutting 
            
    return False

def isValid(seq, i, j):
    judge whether the remaining sequence length is not less than minimum subSequence length
    
def hairpins_are_valid(subSeqs): 
    construct_hairpins
    judge whether all hairpins are valid
```
## Object-oriented Programming
# To be continued
deltaG, secondary structure 7-12 1.5kj



    


