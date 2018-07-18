# HiddenMarkovModel
Predicts protein states using Viterbi Algorithm
This algorithm maximizes the probability of observing the state (transmembrane or soluble) path recursively using dynamic programming.  
It constructs a viterbi path matrix that calculates the highest probability along a single path through the Hidden Markov Model (trellis).
Start and transition probabilities are calculated from the state sequences file, Probabilies for amino acids in the transmembrane state
are calculated from the transmembrane file, and soluble probabilies from the soluble sequence file. 
