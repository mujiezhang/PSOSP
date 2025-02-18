import sys
import os
import re
import math
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from sklearn.cluster import MeanShift
import numpy as np
import argparse


# Set up argument parser
parser = argparse.ArgumentParser(description="PSP (Prophage SOS Predictor)")
parser.add_argument('-hf', '--host_fasta', required=True, help="Host genome fasta file")
parser.add_argument('-vf', '--virus_fasta', required=True, help="Virus genome fasta file")
parser.add_argument('-motif', '--motif_file', required=True, help="19 motifs meme file")
parser.add_argument('-lexa', '--lexa_db', required=True, help="LexA diamond database")
parser.add_argument('-wd', '--working_dir', required=True, help="Output directory path")


def run_prodigal(input_fna: str, output_faa: str, mode: str = "meta"):
    """Execute Prodigal gene prediction"""
    cmd = [
        'prodigal',
        '-i', input_fna,
        '-a', output_faa,
        '-p', mode,
        '-q'
    ]
    subprocess.run(cmd, check=True)

def parse_prodigal_fasta(fasta_file):
    genes = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        description = record.description.replace(' ','').split('#')
        # Parse the description to extract contig_id, start, end, strand information
        contig_id = '_'.join(description[0].split("_")[:-1])
        start = int(description[1])
        end = int(description[2])
        strand = description[3]
        genes.append((contig_id, start, end, strand, record.id, record.seq))
    return genes
def group_genes(genes):
    clusters = []
    current_cluster = [genes[0]]
    for i in range(1, len(genes)):
        prev_gene = genes[i-1]
        curr_gene = genes[i]

        if curr_gene[0] == prev_gene[0] and curr_gene[3] == prev_gene[3] and curr_gene[1] - prev_gene[2] < 20:
            current_cluster.append(curr_gene)
        else:
            clusters.append(current_cluster)
            current_cluster = [curr_gene]
    clusters.append(current_cluster)
    return clusters

def extract_flanking_sequence(contig_seq, gene_start, gene_end, strand, upstream=300, downstream=50):
    if strand == '1':
        start = max(0, gene_start - upstream)
        end = min(len(contig_seq), gene_start + downstream)
        flanking_seq = contig_seq[start:end]
    else:
        start = max(0, gene_end - downstream)
        end = min(len(contig_seq), gene_end + upstream)
        flanking_seq = contig_seq[start:end].reverse_complement()
    return flanking_seq

def fetch_regulon(blastfile):
    b={}
    with open(blastfile,'r') as f:
        for j in f:
            j=j.strip('\n').split('\t')
            gene=j[1].split('|')[0]
            if float(j[4])>=30:
                if gene not in b.keys():
                    b[gene]=[j]
                else:
                    b[gene].append(j)
    return b

def run_host(fna1, host_faa, blast_out, intergenic_region, regulon_region):

    # Load genome sequences
    genome_records = SeqIO.to_dict(SeqIO.parse(fna1, "fasta"))

    # Parse prodigal output
    genes = parse_prodigal_fasta(host_faa)

    # Group genes into clusters
    gene_clusters = group_genes(genes)

    regulon=fetch_regulon(blast_out)

    gene2cluster={}
    cluster_seq={}
    with open(intergenic_region,'w') as inter:
        for cluster in gene_clusters:

            if cluster[0][3] == '1':
                target_gene = cluster[0]
            else:
                target_gene = cluster[-1]

            contig_id = target_gene[0]

            cluster_gene=''
            for i in cluster:
                gene2cluster[i[4]]=target_gene[4]
                cluster_gene+=(i[4]+'-')

            gene_start = target_gene[1]
            gene_end = target_gene[2]
            strand = target_gene[3]
            contig_seq = genome_records[contig_id].seq
            flanking_seq = extract_flanking_sequence(contig_seq, gene_start, gene_end, strand)

            cluster_seq[target_gene[4]]=flanking_seq
            inter.write(f">{target_gene[4]} {cluster_gene} {gene_start}-{gene_end}\n{flanking_seq}\n")

    with open(regulon_region,'w') as reg:
        head=[]
        for i in regulon.keys():
            #head=[]
            for j in regulon[i]:
                protein=j[0]
                target=gene2cluster[protein]
                seq=cluster_seq[target]
                if target not in head:
                    head.append(target)
                    reg.write(f">{protein}-{i} {target}\n{seq}\n")
                else:
                    print(i,regulon[i])

def run_virus(fna2, virus_faa, virus_intergenic_region):
    genome_records = SeqIO.to_dict(SeqIO.parse(fna2, "fasta"))
    genes = parse_prodigal_fasta(virus_faa)
    gene_clusters = group_genes(genes)
    gene2cluster={}
    cluster_seq={}
    with open(virus_intergenic_region,'w') as inter:
        for cluster in gene_clusters:
            if cluster[0][3] == '1':
                target_gene = cluster[0]
            else:
                target_gene = cluster[-1]
            contig_id = target_gene[0]
            cluster_gene=''
            for i in cluster:
                gene2cluster[i[4]]=target_gene[4]
                cluster_gene+=(i[4]+'-')
            gene_start = target_gene[1]
            gene_end = target_gene[2]
            strand = target_gene[3]
            contig_seq = genome_records[contig_id].seq
            flanking_seq = extract_flanking_sequence(contig_seq, gene_start, gene_end, strand)
            cluster_seq[target_gene[4]]=flanking_seq
            inter.write(f">{target_gene[4]} {cluster_gene} {gene_start}-{gene_end}\n{flanking_seq}\n")

base_index={'A':0,'C':2,'G':3,'T':1}
position_matrix={0: [9, 25, 2, 2], 1: [24, 1, 5, 8], 2: [0, 1, 37, 0], 3: [0, 38, 0, 0], 4: [0, 0, 0, 38], 5: [1, 30, 0, 7], 6: [26, 6, 1, 5], 7: [3, 30, 3, 2], 8: [21, 5, 0, 12], 9: [9, 24, 2, 3], 10: [24, 9, 3, 2], 11: [5, 21, 12, 0], 12: [30, 3, 2, 3], 13: [6, 26, 5, 1], 14: [30, 1, 7, 0], 15: [0, 0, 38, 0], 16: [38, 0, 0, 0], 17: [1, 0, 0, 37], 18: [1, 24, 8, 5], 19: [25, 9, 2, 2]}
def calculate_hi_index(sequence):
    hi=0
    for i in range(0,len(sequence)):
        natural=position_matrix[i][base_index[sequence[i]]]+0.5
        consens=max(position_matrix[i])+0.5
        ratio=consens/natural
        hi+=math.log(ratio)
    return hi
    
def analyze_sos_boxes(intergenic_region):
    fna={}
    with open(intergenic_region,'r') as f:
        for i in f:
            i=i.strip('\n')
            if i.startswith('>'):
                key=i
                fna[key]=''
            else:
                fna[key]+=i
    pattern_box={}
    for i in fna.keys():
        seq=fna[i]
        for j in range(2,len(seq)-18):
            if seq[j:j+3]=='CTG' and seq[j+13:j+16]=='CAG':
                pattern_box[i+' '+str(j-2)+'-'+str(j+17)]=seq[j-2:j+18]
    return pattern_box

def cluster_and_analyze_hi(pattern_box, fna1, output_dir):
    a3=[]
    hi=[]
    label=[]
    regu=[]
    for i in pattern_box.keys():
        key=i
        if  re.match("^[ATGC]*$", pattern_box[i]):
                hh=calculate_hi_index(pattern_box[i])
                regu.append(key+'\t'+pattern_box[i])
                a3.append(hh)
    a3, regu = zip(*sorted(zip(a3,regu)))
    x3 = np.array(a3).reshape(-1,1)
    ms3 = MeanShift(max_iter=1000).fit(x3)
    ms_label3=ms3.labels_.flatten().tolist()
    print('MeanShift done')

    host_hi=output_dir+'/'+fna1.split('/')[-1]+'_whole_genome_HI.tsv'
    hi_dic={}
    with open(host_hi,'w') as ff:
        ff.write('HI'+'\t'+'Labels'+'\t'+'protein'+'\t'+'box'+'\n')
        for i in range(0,len(a3)):
            ff.write(str('%.4f'% a3[i])+'\t'+str(ms_label3[i])+'\t'+regu[i]+'\n')
            if ms_label3[i] not in hi_dic.keys():
                hi_dic[ms_label3[i]]=[float('%.4f'% a3[i])]
            else:
                hi_dic[ms_label3[i]].append(float('%.4f'% a3[i]))
    return a3, ms_label3, regu, hi_dic

def predict_phage_type(virus_intergenic_region, fna2, a3, ms_label3, regu, hi_dic, output_dir):
    a4={}
    virus_hi=output_dir+'/'+fna2.split('/')[-1]+'_prediction.tsv'
    hi_final=100
    box_final=''
    with open(virus_intergenic_region,'r') as fii:
        for i in fii:
            i=i.strip('\n')
            if i.startswith('>'):
                key=i.strip('>')
                key2='_'.join(i.strip('>').split()[0].split('_')[:-1])
                if key2 not in a4.keys():
                    a4[key2]=[100,'']
            else:
                seq=i
                for j in range(2,len(seq)-18):
                    seq1=seq[j-2:j+18]
                    if  re.match("^[ATGC]*$", seq1):
                        seq1index=calculate_hi_index(seq1)
                        if seq1index<a4[key2][0]:
                            flag1=1
                            a4[key2][0]=seq1index
                            hi_final=seq1index
                            box_final=seq1
                            a4[key2][1]=key2+'\t'+str(j-2)+'-'+str(j+17)+'\t'+str('%.2f' %seq1index)+'\t'+seq1+'\n'
                    else:
                        pass

    PEAK_EXTENSION_FACTOR = 0.05
    HI_AVERAGE = 15.182220079133675

    # Extract peak boundaries
    keys = list(hi_dic.keys())
    first_peak = hi_dic[keys[0]]
    second_peak = hi_dic[keys[1]]

    # Validate peak sizes
    a3_midpoint = len(a3) // 2
    if len(hi_dic) == 2:
        if len(first_peak) >= a3_midpoint:
            sys.exit('Primary peak exceeds 50% sample threshold')
        
        # Calculate critical thresholds
        interpeak_center = (first_peak[-1] + second_peak[0]) / 2
        confidence_window = (
            interpeak_center - (interpeak_center - first_peak[0]) * PEAK_EXTENSION_FACTOR,
            interpeak_center + (second_peak[-1] - interpeak_center) * PEAK_EXTENSION_FACTOR
        )

    elif len(hi_dic) == 3:
        third_peak = hi_dic[keys[2]]
        if len(first_peak) >= (len(a3) // 3):
            sys.exit('Primary peak exceeds tertiary threshold')
        
        # Determine optimal cluster boundary
        primary_boundary = (first_peak[-1] + second_peak[0]) / 2
        secondary_boundary = (second_peak[-1] + third_peak[0]) / 2
        
        interpeak_center = primary_boundary if abs(HI_AVERAGE - primary_boundary) <= abs(secondary_boundary - HI_AVERAGE) else secondary_boundary
        confidence_window = (
            interpeak_center - (interpeak_center - first_peak[0]) * PEAK_EXTENSION_FACTOR,
            interpeak_center + (third_peak[-1] - interpeak_center) * PEAK_EXTENSION_FACTOR
        )

    else:
        sys.exit(f"Unsupported cluster configuration: {len(hi_dic)} clusters detected")

    # Predictive classification
    phage_class = (
        'SdP (SOS-dependent Prophage)' if hi_final <= confidence_window[0] else
        'SiP (SOS-independent Prophage)' if hi_final >= confidence_window[1] else 
        'SuP (SOS-uncertain Prophage)'
    )

    # Output results
    with open(virus_hi, 'w') as output_file:
        output_file.write(
            f"HI(min)\tbox-seq\tprediction_result\n"
            f"{hi_final}\t{box_final}\t{phage_class}\n"
        )

def run_fimo(memefile, regulon_region, fimo_path):
    os.system('fimo --oc '+fimo_path +' '+memefile+' '+regulon_region)

    fimo_result={}

    fimo_file=fimo_path+'/fimo.tsv'
    with open(fimo_file,'r') as fi:
        for i in fi:
            if i.startswith('TACTGTA'):
                i=i.strip('\n').split('\t')
                if i[5]=='+':
                    name=i[2]+'||'+i[3]+'-'+i[4]
                    fimo_result[name]=i[-1]
    if len(fimo_result)==0:
        sys.exit("No canonical SOS box detected from LexA promoter region thus script stopped")

def main():
    args = parser.parse_args()

    # Assign arguments to variables
    fna1 = args.host_fasta
    fna2 = args.virus_fasta
    memefile = args.motif_file
    lexa_db = args.lexa_db
    output_dir = args.working_dir

    os.makedirs(output_dir, exist_ok=True)
    host_faa = os.path.join(output_dir, fna1.split('/')[-1]+".faa")
    print("Running Prodigal for host...")
    run_prodigal(fna1, host_faa)
    # DIAMOND analysis
    blast_out = os.path.join(output_dir, fna1.split('/')[-1]+".faa_lexa_blast.tsv")
    print("Running DIAMOND for LexA...")

    subprocess.run([
        'diamond', 'blastp',
        '--db', lexa_db,
        '--query', host_faa,
        '--out', blast_out,
        '--threads', '64',  
        '--evalue', '1e-10',      
        '--max-target-seqs', '1', 
        '--id', '30',             
        '--outfmt', '6',          
        'qseqid', 'sseqid', 'qlen', 'slen', 'qcovhsp', 'scovhsp',
        'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend',
        'sstart', 'send', 'evalue', 'bitscore'
    ], check=True)
    if os.path.getsize(blast_out) == 0:
        sys.exit("No LexA homologs detected - analysis terminated")

    intergenic_region=output_dir+'/'+fna1.split('/')[-1]+'_intergenic_region.fna'
    regulon_region=output_dir+'/'+fna1.split('/')[-1]+'_regulon_region.fna'

    run_host(fna1, host_faa, blast_out, intergenic_region, regulon_region)

    virus_faa = os.path.join(output_dir, fna2.split('/')[-1]+".faa")
    print("Running Prodigal for virus...")
    run_prodigal(fna2, virus_faa)
    virus_intergenic_region=output_dir+'/'+fna2.split('/')[-1]+'_intergenic_region.fna'
    run_virus(fna2, virus_faa, virus_intergenic_region)

    fimo_path=output_dir+'/'+fna1.split('/')[-1]+'_fimo'
    run_fimo(memefile, regulon_region, fimo_path)

    pattern_box = analyze_sos_boxes(intergenic_region)
    a3, ms_label3, regu, hi_dic = cluster_and_analyze_hi(pattern_box, fna1, output_dir)
    predict_phage_type(virus_intergenic_region, fna2, a3, ms_label3, regu, hi_dic, output_dir)

if __name__ == "__main__":
    main()
