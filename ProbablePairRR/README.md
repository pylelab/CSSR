# ProbablePairRR #
``ProbablePairRR`` provides a front-end to the algorithm used by
``CSSR --fast 0`` to calculate the base pairing probabilities
from the sequence using thermodyanmics parameters.

## Installation ##
```bash
cd ..
make exe/ProbablePairRR
```

## Usage ##
```bash
exe/ProbablePairRR input.fasta output.rr --sequence
```
input.fasta should include the RNA sequence with characters A, U, C, and G.
output.rr has three columns: 
nucleotide 1, nucleotide 2, probability of pairing.
