import sys
import os
import re
import math
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from sklearn.cluster import MeanShift
import numpy as np

print('import done')
fna1=sys.argv[1]#host fna
fna2=sys.argv[2]#virus faa
memefile=sys.argv[3]#19 motifs meme file
lexa_db=sys.argv[4]# lexa diamond database
output_dir=sys.argv[5]#working dir path


MOTIF_SEQUENCES = [
    "TACTGTATGTTTATACAGTA", "TACTGTGTATATATACAGTA", "TACTGTATATACATACAGCA",
    "TACTGTATATAAAAACAGTA", "TACTGTATAAATAAACAGTT", "TACTGTATATAAAACCAGTT",
    "TACTGTATGAGCATACAGTA", "TACTGTACATCCATACAGTA", "AACTGTTTTTTTATCCAGTA",
    "TGCTGTATATACTCACAGCA", "TACTGTATATTCATTCAGGT", "AACTGTATATACACCCAGGG",
    "TCCTGTTAATCCATACAGCA", "ATCTGTATATATACCCAGCT", "TATTGGCTGTTTATACAGTA",
    "CGCTGGATATCTATCCAGCA", "TACTGTACACAATAACAGTA", "GACTGTATAAAACCACAGCC",
    "ACCTGTAGGATCGTACAGGT"
]

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
def parse_prodigal_output(faa_file: str):
    """Parse Prodigal gene predictions"""
    genes = []
    for record in SeqIO.parse(faa_file, "fasta"):
        desc = record.description.replace(' ', '').split('#')
        contig = '_'.join(desc[0].split("_")[:-1])
        start, end, strand = int(desc[1]), int(desc[2]), desc[3]
        genes.append((contig, start, end, strand, record.id, str(record.seq)))
    return genes

def cluster_genes(genes: list, max_spacing: int = 20):
    """Cluster adjacent genes"""
    if not genes:
        return []
    
    clusters = []
    current = [genes[0]]
    
    for gene in genes[1:]:
        last = current[-1]
        if (gene[0] == last[0] and 
            gene[3] == last[3] and 
            (gene[1] - last[2]) < max_spacing):
            current.append(gene)
        else:
            clusters.append(current)
            current = [gene]
    
    clusters.append(current)
    return clusters

def extract_flanking_sequence(genome: dict, gene: tuple, upstream=300, downstream=50):
    """Extract intergenic sequence"""
    contig_seq = str(genome[gene[0]].seq)
    if gene[3] == '1':
        start = max(0, gene[1] - upstream)
        end = min(len(contig_seq), gene[1] + downstream)
        return contig_seq[start:end]
    else:
        start = max(0, gene[2] - downstream)
        end = min(len(contig_seq), gene[2] + upstream)
        return str(Seq(contig_seq[start:end]).reverse_complement()),str(start),str(end)

os.makedirs(output_dir, exist_ok=True)
host_faa = os.path.join(output_dir, fna1.split('/')[-1]+".faa")
print("Running Prodigal for host...")
#run_prodigal(fna1, host_faa)
# DIAMOND analysis
blast_out = os.path.join(output_dir, fna1.split('/')[-1]+".faa_lexa_blast.tsv")
print("Running DIAMOND for LexA...")

subprocess.run([
    'diamond', 'blastp',
    '--db', lexa_db,
    '--query', host_faa,
    '--out', blast_out,
    '--threads', '64',        # 拆分成两个元素
    '--evalue', '1e-10',      # 拆分成两个元素
    '--max-target-seqs', '1', # 拆分成两个元素
    '--id', '30',             # 拆分成两个元素
    '--outfmt', '6',          # '6'作为单独参数，字段列表可能需要调整
    'qseqid', 'sseqid', 'qlen', 'slen', 'qcovhsp', 'scovhsp',
    'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend',
    'sstart', 'send', 'evalue', 'bitscore'
], check=True)
if os.path.getsize(blast_out) == 0:
    sys.exit("No LexA homologs detected - analysis terminated")
    # Extract intergenic regions
print("Extracting host intergenic regions...")

genes = parse_prodigal_output(host_faa)
clusters = cluster_genes(genes)
genome = SeqIO.to_dict(SeqIO.parse(fna1, "fasta"))
    
intergenic_fna = os.path.join(output_dir, fna1.split('/')[-1]+"_intergenic_region.fna")
with open(intergenic_fna, 'w') as fh:
    for cluster in clusters:
        target = cluster[0] if cluster[0][3] == '1' else cluster[-1]
        seq = extract_flanking_sequence(genome, target)
        fh.write(f">{target[4]}\n{seq}\n")

