# Prophage-SOS-dependency-Predictor(PSP) 
PSP is a novel bioinformatics tool to predict prophage induction modes by analyzing the heterology index (HI) of LexA protein binding to target DNA, classifying prophages into SOS-dependent (SdPs) and SOS-independent (SiPs).

## Dependencies
* PSP is a Python script that relies on:
```Biopython
DIAMOND
MEME
Python3
scikit-learn
Prodigal
```

## Installation
(1) conda (recommended)
```
conda create -n psp
conda activate psp
conda install zhangmujie::psp
```
usage: ```psp -h```

(2) git (install dependencies mentioned above first)
```
git clone https://github.com/mujiezhang/PSP.git
cd PSP
python psp.py -h
```

## Input files
PSP needs two files as inputs,i.e.,
* ```-hf```: a host genome in fasta format
* ```-vf```: a single viral genome in fasta format 

other parameters
* ```-wd```: woking path to save result files

## How to run
The users can only specify the required parameters:
* install through git
```
python psp.py -hf host-genome.fasta -vf virus-genome.fasta -wd output_dir
```
* install through conda
```
psp -hf host-genome.fasta -vf virus-genome.fasta -wd output_dir
```

for example (install through git):
```
python psp.py -hf test/E.coli-HS.fasta -vf test/phiECO1.fasta -wd test/
```

Running this example with one core takes approximately two minutes. And result file is stored in file ```*_prediction.tsv```

## Attention
* PSP is designed for complete host and corresponding complete virus for that host. Using incomplete genome as input may influence the prediction accuracy.

## Citation
''''''
