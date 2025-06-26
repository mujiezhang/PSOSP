
# PSOSP: Prophage-SOS-dependency-Predictor
PSOSP (<b>P</b>rophage <b>SOS</b> dependency <b>P</b>redictor) is a novel bioinformatics tool to predict prophage induction modes by analyzing the heterology index (HI) of LexA protein binding to target DNA, classifying prophages into SOS-dependent (SdPs) and SOS-independent (SiPs).

![PSOSP-LOGO](https://github.com/user-attachments/assets/3f02c2e6-5596-4293-970f-45e60139e6a9)



## Table of contents
<!-- TOC -->
- [PSOSP: Prophage-SOS-dependency-Predictor](#psosp-prophage-sos-dependency-predictor)
- [Webserver](#webserver)
- [Background](#background)
- [Workflow of PSOSP](#workflow-of-psosp)
- [Experimental validation](#experimental-validation)
- [Input requirements](#input-requirements)
- [Dependencies](#dependencies)
- [Installation](#installation)
- [Input files](#input-files)
- [How to run](#how-to-run)
- [Outputs](#outputs)
- [Citation](#citation)

<!-- /TOC -->

## Webserver
We provide an online platform (PSOSP) for rapid prediction of bacteriophage induction modes: **https://vee-lab.sjtu.edu.cn/PSOSP/**. There you can upload your host and virus genomes and get the prediction results.

## Background
Temperate phages integrate into the bacterial host genome as prophages. Under normal conditions, the LexA protein binds to the SOS box within the prophage, repressing the expression of phage-related genes and maintaining the lysogenic state. Upon external stimuli (such as exposure to DNA-damaging agents), the RecA protein is activated, leading the self-cleavage of LexA and its dissociation from the SOS box. This relieves the prophage repression, triggering the temperate phage to enter the lytic cycle and thereby facilitating its proliferation.

![psosp-theory](https://github.com/user-attachments/assets/654a77e1-dbb6-44bb-9719-0fe4fca7519c)


## Workflow of PSOSP
- LexA & Canonical SOS Box (CBS) Identification : Scanning the host genome to identify LexA protein and canonical SOS boxes (CSBs) located upstream of the lexA gene
  
- Heterology Index (HI) Calculation: Identifying potential SOS boxes (PSBs) across bacterial genomes, calculating the Heterology Index (HI) for each PSB and establishing classification thresholds (HI<sub>c1</sub> and HI<sub>c2</sub>) via Mean Shift clustering results
  
- PSB scan in prophage: Scanning PSBs within prophage promoter regions and determining of the minimum HI (HI<sub>min</sub>)
  
- Prophage categoriation : Evaluating the ability of LexA binding to prophage promoter regions by comparing HImin with thresholds
   - HI<sub>min</sub> ≤ HI<sub>c1</sub> → **SdP** (SOS-dependent Prophage)
   - HI<sub>min</sub> ≥ HI<sub>c2</sub> → **SiP** (SOS-independent Prophage)
   - HI<sub>c1</sub> < HI<sub>min</sub> < HI<sub>c2</sub> → **SuP** (SOS-uncertain Prophage)


![PSOSP_workflow](https://github.com/user-attachments/assets/c3a0334f-ce4f-4533-960a-ae8c19b71514)



## Experimental validation
We have validated PSOSP's accuracy using 14 experimentally confirmed bacteriophages spanning 10 viral families (including 2 Peduoviridae, 3 Inoviridae, and 9 distinct novel families), with their hosts covering 7 bacterial genera (Salmonella, Escherichia, Vibrio, Pseudomonas, Serratia_J, Hafnia, and Shewanella) across 3 bacterial orders (Enterobacterales, Enterobacterales_A, and Pseudomonadales). Significantly, all PSOSP predictions for these bacteriophages showed complete consistency with experimental evidence, demonstrating the tool's versatility and reliability across broad taxonomic ranges.

![experiment_validation](https://github.com/user-attachments/assets/f39bc3c6-a18b-4bf4-9459-d14176d76289)


## Input requirements
**For host**:
- Host Taxon Suitability: PSOSP is primarily suitable for Gammaproteobacteria. If your host belongs to another taxon, PSOSP is unlikely to produce meaningful results.

- Genome Quality: **We advise using host genomes with a completeness score above 90%,** since low-quality genomes may lose the LexA protein and lead to poor results. You can assess the completeness of your genome assembly using [**CheckM2**](https://github.com/chklovski/CheckM2).

- Multi-Contig Genomes: If the host genome consists of multiple contigs, ensure the input host genome file contains all contigs (i.e., provide the genome assembly as a single multi-contig file).

**For prophage**:

- Genome Quality: PSOSP utilizes CheckV for quality assessment. **Predictions for viruses with a CheckV-estimated completeness >90% are relatively reliable.** If you are certain your viral genome is complete, you may disregard the CheckV results in the output file.

- Multiple Inputs: The input viral genome file can contain sequences for multiple viruses (prophages).

- Host Association: The input viruses must be prophages integrated within the specific input host genome. Predicting associations for mismatched virus-host pairs is meaningless.
  

## Dependencies
* PSOSP is a Python script that relies on:
```
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
* ```-faa```: host protein sequences in fasta format (optional)

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

for example (install through conda):
```
psosp -hf test/E.coli-HS.fasta -vf test/phiECO1.fasta -wd test/
```


## Outputs
In this example, the results of PSOSP's analysis will be written to the `test/result` directory, which will look like this:
```
test/result
├── virus_wp2-phage-sp1-sp2-sp3_checkv
├── host_wp2.fna.faa_lexa_blast.tsv
├── host_wp2.fna_whole_genome_HI.tsv
├── host_wp2_prodigal.faa
╰── virus_wp2-phage-sp1-sp2-sp3_prediction.tsv
```

1. `virus_wp2-phage-sp1-sp2-sp3_checkv`: this directory contain results of [**CheckV**](https://bitbucket.org/berkeleylab/checkv/src/master/)
2. `host_wp2.fna.faa_lexa_blast.tsv`: blast result of LexA protein
3. `host_wp2.fna_whole_genome_HI.tsv`: HI clusters of all potential SOS box in host genome using MeanShift method.
4. `host_wp2_prodigal.faa`: protein sequences produced by prodigal
5. `virus_wp2-phage-sp1-sp2-sp3_prediction.tsv`: **prediction results of PSOSP**

A detailed overview of **`virus_wp2-phage-sp1-sp2-sp3_prediction.tsv`**:

|                **host**                |                **virus**                |                **prediction_result**                |                **prediction_quality**                |                **completeness**                |                **contamination**                |                **viral-HI(min)**                |                **box-seq**                |                **box-seq_start_pos**                |                **box-seq_strand**                |                **confidence_window_lower**                |                **confidence_window_upper**                |                **blast_status**                |                **fimo_status**                |
|:--------------------------------------:|:----------------------------------------:|:------------------------------------------------:|:------------------------------------------------:|:--------------------------------------------:|:--------------------------------------------:|:--------------------------------------------:|:----------------------------------------:|:------------------------------------------------:|:--------------------------------------------:|:----------------------------------------------------:|:----------------------------------------------------:|:--------------------------------------------:|:------------------------------------------:|
|               host_wp2.fna              |                    sp1                   |          SiP (SOS-independent Prophage)          |                      High                         |                    100.0                     |                     0.0                      |                   14.2017                    |         CACTGTATTATTATACCACA              |                     29652                      |                      -                       |                      11.8522                       |                      13.3752                       |                   Blast_OK                   |                   Fimo_OK                  |
|               host_wp2.fna              |                    sp2                   |           SdP (SOS-dependent Prophage)           |                      High                         |                    100.0                     |                     0.0                      |                   11.7239                    |         TTATGTATGTATATTCAGCA              |                     17757                      |                      -                       |                      11.8522                       |                      13.3752                       |                   Blast_OK                   |                   Fimo_OK                  |
|               host_wp2.fna              |                    sp3                   |          SiP (SOS-independent Prophage)          |                      High                         |                    99.62                     |                     0.0                      |                   15.4715                    |         CACTGTATAAAAAAACATAC              |                     1225                       |                      -                       |                      11.8522                       |                      13.3752                       |                   Blast_OK                   |                   Fimo_OK                  |

- `host`: Input host filename
- `virus`: Viral identifier in input FASTA
- `prediction_result`: prediction induction mode of PSOSP. **SiPs**: SOS-independent Prophage; **SuPs**: SOS-uncertain Prophage; **SdPs**: SOS-dependent Prophage
- `prediction quality`: **`High`** for viral completeness between 90%-100%; **`Medium`** for viral completeness between 50%-90%; **'Low'** for viral completeness lower than 50%
- `completeness`: Estimated viral completeness (CheckV)
- `contamination`: Estimated viral contamination (CheckV)
- `viral-HI(min)`: Minimal HI in viral genome
- `box-seq`: sequence of potential sos box with minimal HI
- `box-seq_start_pos`: box start position in viral genome
- `box-seq_strand`: + (forward) or - (reverse)
- `confidence_window_lower`: Threshold HI<sub>c1</sub> (HI<sub>min</sub> ≤ HI<sub>c1</sub> → `SdP`) 
- `confidence_window_upper`: Threshold HI<sub>c2</sub> (HI<sub>min</sub> ≥ HI<sub>c2</sub> → `SiP`) 
- `blast_status`: 'Blast_OK' (LexA homologs found) or fails if absent
- `fimo_status`: 'Fimo_OK' (SOS-box detected upstream of LexA) or fails if absent

## Citation
''''''
