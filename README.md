<div align="center">

# Prophage-SOS-dependency-Predictor(PSOSP)

<img src="https://github.com/user-attachments/assets/e8b7d80a-5a7b-4630-bef1-125c971b0bf2" alt="psp-logo" width="600" />


</div>
PSOSP (<b>P</b>rophage <b>SOS</b> dependency <b>P</b>redictor) is a novel bioinformatics tool to predict prophage induction modes by analyzing the heterology index (HI) of LexA protein binding to target DNA, classifying prophages into SOS-dependent (SdPs) and SOS-independent (SiPs).
[toc]
# Table of contents
-[PSOSP](# Prophage-SOS-dependency-Predictor(PSOSP))

## Webserver
我们提供PSP的在线网站供用户进行快速预测噬菌体的sos 调控类型：https://vee-lab.sjtu.edu.cn/PSP/

## Dependencies
* PSOSP is a Python script that relies on:
```Biopython
DIAMOND=2.1.11
MEME=5.5.5
Python3
scikit-learn=1.6.1
Prodigal=2.6.3
biopython=1.85
checkv=1.0.3
```

## Installation
(1) conda (recommended)
```
conda config --set channel_priority flexible
conda create -n psosp
conda activate psosp
conda install zhangmujie::psosp
```
usage: ```psosp -h```

(2) git (install dependencies mentioned above first)
```
git clone https://github.com/mujiezhang/PSOSP.git
cd PSOSP
python psosp.py -h
```

## Input files
PSOSP needs two files as inputs,i.e.,
* ```-hf```: a host genome in fasta format
* ```-vf```: a single viral genome in fasta format 

other parameters
* ```-wd```: woking path to save result files

## How to run
The users can only specify the required parameters:
* install through git
```
python psosp.py -hf host-genome.fasta -vf virus-genome.fasta -wd output_dir
```
* install through conda
```
psosp -hf host-genome.fasta -vf virus-genome.fasta -wd output_dir
```

for example (install through git):
```
python psosp.py -hf test/E.coli-HS.fasta -vf test/phiECO1.fasta -wd test/
```

Running this example with one core takes approximately two minutes. And result file is stored in file ```*_prediction.tsv```

## Attention
* PSOSP is designed for complete host and corresponding complete virus for that host. Using incomplete genome as input may influence the prediction accuracy.

## Citation
''''''
