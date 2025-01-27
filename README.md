# Prophage-SOS-dependency-Predictor(PSP) 
PSP is a novel bioinformatics tool to predict prophage induction modes by analyzing the heterology index (HI) of LexA protein binding to target DNA, classifying prophages into SOS-dependent (SdPs) and SOS-independent (SiPs).
## Dependencies
* PSP is a Python script that relies on:
```Biopython
DIAMOND
MEME
Python3
scikit-learn
```
## Input files
PSP needs four files as inputs,i.e.,
* -hf: host genome in fasta format
* -vf: a single viral genome in fasta format
* -motif: a motif file provided by psp as ```COBRA_joining_status.txt``` 


#
## Attention

