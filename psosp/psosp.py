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
parser = argparse.ArgumentParser(description="PSOSP (Prophage SOS-dependency Predictor)")
subparsers = parser.add_subparsers(dest='command', help='Available commands')

# Add test command
test_parser = subparsers.add_parser('test', help='Run test with sample data')

# Add normal prediction command (make it default)
predict_parser = subparsers.add_parser('predict', help='Run prophage prediction (default)')
predict_parser.add_argument('-hf', '--host_fasta', required=True, help="Host genome fasta file")
predict_parser.add_argument('-vf', '--virus_fasta',required=True, help="Virus genome fasta file (can contain multiple viruses)")
predict_parser.add_argument('-wd', '--working_dir', required=True, help="Output directory path")
predict_parser.add_argument('-hfaa', '--host_faa', required=False, help="Host faa file path (optional, will run prodigal if not provided)")
predict_parser.add_argument('-db', '--checkv_db', required=False, help="checkv reference database path (optional)")

# For backward compatibility, also add arguments to main parser
parser.add_argument('-hf', '--host_fasta', required=False, help="Host genome fasta file")
parser.add_argument('-vf', '--virus_fasta',required=False, help="Virus genome fasta file (can contain multiple viruses)")
parser.add_argument('-wd', '--working_dir', required=False, help="Output directory path")
parser.add_argument('-hfaa', '--host_faa', required=False, help="Host faa file path (optional, will run prodigal if not provided)")
parser.add_argument('-db', '--checkv_db', required=False, help="checkv reference database path (optional)")
parser.add_argument('-v', '--version', action='version', version='PSOSP v1.1')
def run_checkv(input_fna: str, output_dir: str, db_dir: str ):
    """Execute CheckV analysis"""
    cmd = [
        'checkv', 'end_to_end',
        input_fna,
        output_dir,
        '-d', db_dir,
        '--remove_tmp',
        '--restart',
        '-t', '4',
        '--quiet'
    ]
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(f"CheckV completed: {output_dir}")
    except subprocess.CalledProcessError as e:
        print(f"CheckV command failed: {' '.join(cmd)}")
        print(f"Return code: {e.returncode}")
        print(f"stdout: {e.stdout}")
        print(f"stderr: {e.stderr}")
        raise

def run_prodigal(input_fna: str, output_faa: str, mode: str = "meta"):
    """Execute Prodigal gene prediction"""
    cmd = [
        'prodigal',
        '-i', input_fna,
        '-a', output_faa,
        '-p', mode,
        '-q'
    ]
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(f"Prodigal completed: {output_faa}")
    except subprocess.CalledProcessError as e:
        print(f"Prodigal command failed: {' '.join(cmd)}")
        print(f"Return code: {e.returncode}")
        print(f"stdout: {e.stdout}")
        print(f"stderr: {e.stderr}")
        raise

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

def run_virus_multi(fna2, virus_faa, virus_intergenic_regions_dir):
    """Modified function to support processing multiple viruses"""
    genome_records = SeqIO.to_dict(SeqIO.parse(fna2, "fasta"))
    genes = parse_prodigal_fasta(virus_faa)
    
    # Group genes by virus ID
    genes_by_virus = {}
    for gene in genes:
        virus_id = gene[0]  # contig_id is the virus ID
        if virus_id not in genes_by_virus:
            genes_by_virus[virus_id] = []
        genes_by_virus[virus_id].append(gene)
    
    virus_intergenic_files = {}
    
    # Create intergenic region files for each virus
    for virus_id, virus_genes in genes_by_virus.items():
        gene_clusters = group_genes(virus_genes)
        gene2cluster = {}
        cluster_seq = {}
        
        virus_intergenic_file = os.path.join(virus_intergenic_regions_dir, f"{virus_id}_intergenic_region.fna")
        virus_intergenic_files[virus_id] = virus_intergenic_file
        
        with open(virus_intergenic_file, 'w') as inter:
            for cluster in gene_clusters:
                if cluster[0][3] == '1':
                    target_gene = cluster[0]
                else:
                    target_gene = cluster[-1]
                
                contig_id = target_gene[0]
                cluster_gene = ''
                for i in cluster:
                    gene2cluster[i[4]] = target_gene[4]
                    cluster_gene += (i[4] + '-')
                
                gene_start = target_gene[1]
                gene_end = target_gene[2]
                strand = target_gene[3]
                contig_seq = genome_records[contig_id].seq
                flanking_seq = extract_flanking_sequence(contig_seq, gene_start, gene_end, strand)
                cluster_seq[target_gene[4]] = flanking_seq
                inter.write(f">{target_gene[4]} {cluster_gene} {gene_start}-{gene_end} strand={strand}\n{flanking_seq}\n")
    
    return virus_intergenic_files

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
    print('[6/6]MeanShift done')

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

def predict_phage_type_multi(virus_intergenic_files, fna1, fna2, a3, ms_label3, regu, hi_dic, output_dir, has_blast_results, has_fimo_results, checkv_output_dir):
    """
    Modified prediction function that can handle multiple viruses and write results to the same file
    """
    host_basename = os.path.basename(fna1)
    virus_basename = os.path.basename(fna2)
    virus_name, _ = os.path.splitext(virus_basename)
    prediction_file = os.path.join(output_dir, f"{virus_name}_prediction.tsv")
    
    # Read CheckV quality assessment results (if available)
    checkv_quality = {}
    if checkv_output_dir:
        quality_summary_file = os.path.join(checkv_output_dir, 'quality_summary.tsv')
        if os.path.exists(quality_summary_file):
            with open(quality_summary_file, 'r') as f:
                header = f.readline().strip().split('\t')
                try:
                    completeness_idx = header.index('completeness')
                    contamination_idx = header.index('contamination')
                    contig_id_idx = header.index('contig_id')
                    
                    for line in f:
                        fields = line.strip().split('\t')
                        if len(fields) > max(completeness_idx, contamination_idx, contig_id_idx):
                            contig_id = fields[contig_id_idx]
                            completeness = fields[completeness_idx] if fields[completeness_idx] != 'NA' else 'NA'
                            contamination = fields[contamination_idx] if fields[contamination_idx] != 'NA' else 'NA'
                            checkv_quality[contig_id] = {
                                'completeness': completeness,
                                'contamination': contamination
                            }
                except ValueError:
                    print("Warning: Could not find completeness/contamination columns in CheckV output")
        else:
            print(f"Warning: CheckV quality summary file not found: {quality_summary_file}")
    # If no checkv_output_dir is provided, checkv_quality remains empty
    
    # Create or write prediction results file
    with open(prediction_file, 'w') as output_file:
        header = "host\tvirus\tprediction_result\tprediction_quality\tcompleteness\tcontamination\tviral-HI(min)\tbox-seq\tbox-seq_start_pos\tbox-seq_strand\tconfidence_window_lower\tconfidence_window_upper\tblast_status\tfimo_status\n"
        output_file.write(header)
        
        # Generate prediction results for each virus
        for virus_id, virus_intergenic_file in virus_intergenic_files.items():
            # 1. Calculate HI(min) and box-seq for each virus
            hi_final = 100
            box_final = ''
            box_start_pos = 'NA'
            box_contig = 'NA'
            box_strand = 'NA'
            
            if os.path.exists(virus_intergenic_file):
                with open(virus_intergenic_file, 'r') as fii:
                    current_key = ''
                    current_gene_info = {}
                    for line in fii:
                        line = line.strip('\n')
                        if line.startswith('>'):
                            current_key = line.strip('>')
                            # Parse gene information: format is "{gene_id} {cluster_gene} {gene_start}-{gene_end} strand={strand}"
                            parts = current_key.split()
                            if len(parts) >= 4:
                                gene_id = parts[0]
                                gene_coords = parts[2].split('-')
                                strand_info = parts[3] if parts[3].startswith('strand=') else 'strand=1'
                                strand = strand_info.split('=')[1] if '=' in strand_info else '1'
                                
                                if len(gene_coords) == 2:
                                    try:
                                        gene_start = int(gene_coords[0])
                                        gene_end = int(gene_coords[1])
                                        # Extract contig_id from gene_id (assuming format is contig_gene)
                                        contig_id = '_'.join(gene_id.split('_')[:-1])
                                        current_gene_info = {
                                            'gene_start': gene_start,
                                            'gene_end': gene_end,
                                            'strand': strand,
                                            'contig_id': contig_id,
                                            'gene_id': gene_id
                                        }
                                    except ValueError:
                                        current_gene_info = {}
                        else:
                            seq = line
                            if current_gene_info:
                                for j in range(2, len(seq)-18):
                                    seq1 = seq[j-2:j+18]
                                    if re.match("^[ATGC]*$", seq1):
                                        seq1index = calculate_hi_index(seq1)
                                        if seq1index < hi_final:
                                            hi_final = seq1index
                                            box_final = seq1
                                            
                                            # Calculate absolute position in virus genome
                                            gene_start = current_gene_info['gene_start']
                                            gene_end = current_gene_info['gene_end']
                                            strand = current_gene_info['strand']
                                            contig_id = current_gene_info['contig_id']
                                            
                                            # Calculate absolute position based on extract_flanking_sequence logic
                                            if strand == '1':  # Forward strand gene
                                                # flanking region: start = max(0, gene_start - 300), end = min(len, gene_start + 50)
                                                flanking_start = max(0, gene_start - 300)  # upstream=300
                                                box_abs_pos = flanking_start + (j-2)
                                            else:  # Reverse strand gene
                                                # flanking region: start = max(0, gene_end - 50), end = min(len, gene_end + 300)
                                                flanking_start = max(0, gene_end - 50)  # downstream=50
                                                # For reverse strand, sequence is reverse complemented, so special handling is needed
                                                # box position j-2 in reverse complement sequence needs to be converted to original position
                                                flanking_length = len(seq)
                                                box_abs_pos = flanking_start + flanking_length - (j-2) - 20  # 20 is box length
                                            
                                            box_start_pos = box_abs_pos
                                            box_contig = contig_id
                                            box_strand = '+' if strand == '1' else '-'
            
            # 2. Determine prediction status based on BLAST and FIMO results
            phage_class = 'NA'
            confidence_window_lower = 'NA'
            confidence_window_upper = 'NA'
            blast_status = 'Blast_OK' if has_blast_results else 'No_Blast_Results'
            fimo_status = 'NA' # Default for when BLAST fails

            if not has_blast_results:
                phage_class = 'No LexA homologs detected - cannot predict'
            elif not has_fimo_results:
                fimo_status = 'No_Fimo_Results'
                phage_class = 'No canonical SOS box detected - cannot predict'
            else:
                # Only proceed with full prediction when both BLAST and FIMO are successful
                fimo_status = 'Fimo_OK'
                PEAK_EXTENSION_FACTOR = 0.05
                HI_AVERAGE = 15.182220079133675
                
                if not hi_dic or len(hi_dic.keys()) < 2:
                    phage_class = "Prediction failed: Not enough host HI clusters."
                else:
                    keys = list(hi_dic.keys())
                    first_peak = hi_dic[keys[0]]
                    second_peak = hi_dic[keys[1]]
                    a3_midpoint = len(a3) // 2
                    
                    confidence_window = None

                    if len(hi_dic) == 2:
                        if len(first_peak) >= a3_midpoint:
                            phage_class = 'Primary peak exceeds 50% sample threshold - uncertain prediction'
                        else:
                            interpeak_center = (first_peak[-1] + second_peak[0]) / 2
                            confidence_window = (
                                interpeak_center - (interpeak_center - first_peak[0]) * PEAK_EXTENSION_FACTOR,
                                interpeak_center + (second_peak[-1] - interpeak_center) * PEAK_EXTENSION_FACTOR
                            )
                    elif len(hi_dic) == 3:
                        third_peak = hi_dic[keys[2]]
                        if len(first_peak) >= (len(a3) // 3):
                            phage_class = 'Primary peak exceeds tertiary threshold - uncertain prediction'
                        else:
                            primary_boundary = (first_peak[-1] + second_peak[0]) / 2
                            secondary_boundary = (second_peak[-1] + third_peak[0]) / 2
                            interpeak_center = primary_boundary if abs(HI_AVERAGE - primary_boundary) <= abs(secondary_boundary - HI_AVERAGE) else secondary_boundary
                            confidence_window = (
                                interpeak_center - (interpeak_center - first_peak[0]) * PEAK_EXTENSION_FACTOR,
                                interpeak_center + (third_peak[-1] - interpeak_center) * PEAK_EXTENSION_FACTOR
                            )
                    else:
                        phage_class = f"Unsupported cluster configuration: {len(hi_dic)} clusters detected"

                    if confidence_window:
                        confidence_window_lower = f'{confidence_window[0]:.4f}'
                        confidence_window_upper = f'{confidence_window[1]:.4f}'
                        phage_class = (
                            'SdP (SOS-dependent Prophage)' if hi_final <= confidence_window[0] else
                            'SiP (SOS-independent Prophage)' if hi_final >= confidence_window[1] else 
                            'SuP (SOS-uncertain Prophage)'
                        )
            
            # 3. Format output data
            hi_final_str = f'{hi_final:.4f}' if hi_final != 100 else 'NA'
            box_final_str = box_final if box_final else 'NA'
            box_start_pos_str = str(box_start_pos) if box_start_pos != 'NA' else 'NA'
            box_strand_str = box_strand if box_strand != 'NA' else 'NA'
            
            # Get CheckV quality information
            if checkv_output_dir and virus_id in checkv_quality:
                completeness = checkv_quality[virus_id]['completeness']
                contamination = checkv_quality[virus_id]['contamination']
            elif checkv_output_dir:
                # CheckV was run but no results for this virus
                completeness = 'NA'
                contamination = 'NA'
            else:
                # CheckV was not run (no --db provided)
                completeness = '-'
                contamination = '-'
            
            # Determine prediction quality based on completeness
            if completeness == 'NA' or completeness == '-':
                prediction_quality = 'low'
            else:
                try:
                    completeness_val = float(completeness)
                    if completeness_val >= 90:
                        prediction_quality = 'high'
                    elif completeness_val >= 50:
                        prediction_quality = 'medium'
                    else:
                        prediction_quality = 'low'
                except ValueError:
                    prediction_quality = 'low'
            
            line = f"{host_basename}\t{virus_id}\t{phage_class}\t{prediction_quality}\t{completeness}\t{contamination}\t{hi_final_str}\t{box_final_str}\t{box_start_pos_str}\t{box_strand_str}\t{confidence_window_lower}\t{confidence_window_upper}\t{blast_status}\t{fimo_status}\n"
            output_file.write(line)

def run_fimo(memefile, regulon_region, fimo_path):
    os.system('fimo --oc '+fimo_path +' '+memefile+' '+regulon_region + '> /dev/null 2>&1')

    fimo_result={}

    fimo_file=fimo_path+'/fimo.tsv'
    if os.path.exists(fimo_file):
        with open(fimo_file,'r') as fi:
            for i in fi:
                if i.startswith('TACTGTA'):
                    i=i.strip('\n').split('\t')
                    if i[5]=='+':
                        name=i[2]+'||'+i[3]+'-'+i[4]
                        fimo_result[name]=i[-1]
    
    # Return fimo results without exiting the program
    return fimo_result

def run_test():
    """Run test with sample data"""
    # This script will be located inside the 'psosp' package directory,
    # so __file__ will point there.
    script_path = os.path.abspath(__file__)
    script_dir = os.path.dirname(script_path)
    
    # Construct paths relative to the script's location
    test_virus = os.path.join(script_dir, "test/data/virus_wp2-phage-sp1-sp2-sp3.fna")
    test_host = os.path.join(script_dir, "test/data/host_wp2.fna")
    test_output = "temp"
    
    print("=" * 60)
    print("PSOSP Test Mode")
    print("=" * 60)
    print(f"Running test with:")
    print(f"  Virus file: {test_virus}")
    print(f"  Host file: {test_host}")
    print(f"  Output directory: {test_output}")
    print("=" * 60)
    
    # Check if test files exist
    if not os.path.exists(test_virus):
        print(f"Error: Test virus file not found: {test_virus}")
        return
    if not os.path.exists(test_host):
        print(f"Error: Test host file not found: {test_host}")
        return
    
    try:
        # Create test arguments
        import types
        test_args = types.SimpleNamespace()
        test_args.host_fasta = test_host
        test_args.virus_fasta = test_virus
        test_args.working_dir = test_output
        test_args.host_faa = None
        test_args.checkv_db = None
        
        # Run the prediction with test data
        run_prediction(test_args)
        
        print("\n" + "=" * 60)
        print("Test completed successfully!")
        print("=" * 60)
        
    except Exception as e:
        print(f"\nTest failed with error: {e}")
        import traceback
        traceback.print_exc()
    
    finally:
        # Clean up temp directory
        if os.path.exists(test_output):
            print(f"Cleaning up temporary directory: {test_output}")
            import shutil
            shutil.rmtree(test_output)
            print("Temporary files cleaned up.")

def run_prediction(args):
    """Run the main prediction pipeline"""
    # Assign arguments to variables
    fna1 = args.host_fasta
    fna2 = args.virus_fasta
    output_dir = args.working_dir
    host_faa_provided = args.host_faa
    checkv_db = args.checkv_db

    # The script_dir will correctly point to the 'psosp' package directory
    # where the 'files' directory is also located.
    script_path = os.path.abspath(__file__)
    script_dir = os.path.dirname(script_path)
    memefile = os.path.join(script_dir, "files/19-motifs-meme.txt")
    lexa_db = os.path.join(script_dir, "files/uniprot_swiss_prot_LexA.dmnd")

    # Check and create output directory
    if os.path.exists(output_dir):
        print(f"Using existing output directory: {output_dir}")
    else:
        print(f"Creating output directory: {output_dir}")
        os.makedirs(output_dir, exist_ok=True)
    
    # Create temporary directory for virus intergenic region files
    virus_intergenic_dir = os.path.join(output_dir, 'virus_intergenic_temp')
    os.makedirs(virus_intergenic_dir, exist_ok=True)
    
    # Check the number of sequences in virus file
    virus_records = list(SeqIO.parse(fna2, "fasta"))
    num_viruses = len(virus_records)
    print(f"Detected {num_viruses} virus sequences")
    
    # Run CheckV analysis (only if --db is provided)
    virus_basename = os.path.basename(fna2)
    virus_name, _ = os.path.splitext(virus_basename)
    checkv_output_dir = None
    if checkv_db:
        checkv_output_dir = os.path.join(output_dir, f"{virus_name}_checkv")
        if os.path.exists(os.path.join(checkv_output_dir, 'quality_summary.tsv')):
            print(f"[1/6]Found existing CheckV results: {checkv_output_dir}")
        else:
            print("[1/6]Running CheckV virus quality assessment...")
            os.makedirs(checkv_output_dir, exist_ok=True)
            run_checkv(fna2, checkv_output_dir, checkv_db)
    else:
        print("[1/6]Skipping CheckV analysis (no database provided)")
        checkv_output_dir = None
    
    # Process host faa file
    if host_faa_provided:
        host_faa = host_faa_provided
        print("[2/6]Using provided host faa file:", host_faa)
    else:
        host_basename = os.path.basename(fna1)
        host_name, _ = os.path.splitext(host_basename)
        host_faa = os.path.join(output_dir, f"{host_name}_prodigal.faa")
        
        if os.path.exists(host_faa):
            print(f"[2/6]Found existing host protein file: {host_faa}")
        else:
            print("[2/6]Running Prodigal for host...")
            run_prodigal(fna1, host_faa)
    
    # DIAMOND analysis
    blast_out = os.path.join(output_dir, fna1.split('/')[-1]+".faa_lexa_blast.tsv")
    print("[3/6]Running DIAMOND for LexA...")

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
    ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
    # Check blast results
    has_blast_results = os.path.getsize(blast_out) > 0
    
    # Run virus analysis
    virus_faa = os.path.join(output_dir, fna2.split('/')[-1]+".faa")
    print("[4/6]Running Prodigal for virus...")
    run_prodigal(fna2, virus_faa)
    
    # Use new multi-virus processing function
    virus_intergenic_files = run_virus_multi(fna2, virus_faa, virus_intergenic_dir)
    
    if not has_blast_results:
        print("No LexA homologs detected - will write this to prediction file")
        # Write prediction file without blast results
        predict_phage_type_multi(virus_intergenic_files, fna1, fna2, [], [], [], {}, output_dir, False, False, checkv_output_dir)
        
        # Clean up temporary files
        os.system(f'rm -rf {virus_intergenic_dir}')
        os.system('rm '+virus_faa)
        return

    intergenic_region=output_dir+'/'+fna1.split('/')[-1]+'_intergenic_region.fna'
    regulon_region=output_dir+'/'+fna1.split('/')[-1]+'_regulon_region.fna'

    run_host(fna1, host_faa, blast_out, intergenic_region, regulon_region)

    fimo_path=output_dir+'/'+fna1.split('/')[-1]+'_fimo'
    fimo_result = run_fimo(memefile, regulon_region, fimo_path)
    
    # Check fimo results
    has_fimo_results = len(fimo_result) > 0
    
    if not has_fimo_results:
        print("No canonical SOS box detected - will write this to prediction file")
        predict_phage_type_multi(virus_intergenic_files, fna1, fna2, [], [], [], {}, output_dir, True, False, checkv_output_dir)
    else:
        print("[5/6]Both BLAST and FIMO successful - proceeding with normal prediction")
        pattern_box = analyze_sos_boxes(intergenic_region)
        a3, ms_label3, regu, hi_dic = cluster_and_analyze_hi(pattern_box, fna1, output_dir)
        predict_phage_type_multi(virus_intergenic_files, fna1, fna2, a3, ms_label3, regu, hi_dic, output_dir, True, True, checkv_output_dir)
    
    # Clean up temporary files
    os.system('rm '+intergenic_region)
    os.system('rm '+regulon_region)
    os.system(f'rm -rf {virus_intergenic_dir}')
    os.system('rm -rf '+fimo_path)
    os.system('rm '+virus_faa)
    
    print("PSOSP completed")

def main():
    args = parser.parse_args()
    
    # Handle different commands
    if args.command == 'test':
        run_test()
    elif args.command == 'predict':
        run_prediction(args)
    else:
        # For backward compatibility - if no subcommand is used, treat as normal prediction
        # Check if required arguments are provided
        if not args.host_fasta or not args.virus_fasta or not args.working_dir:
            print("Error: Required arguments missing.")
            print("Usage:")
            print("  python psosp.py test                              # Run test with sample data")
            print("  python psosp.py predict -vf <virus> -hf <host> -wd <output>  # Run prediction")
            print("  python psosp.py -vf <virus> -hf <host> -wd <output>          # Run prediction (legacy)")
            parser.print_help()
            return
        run_prediction(args)

if __name__ == "__main__":
    main()
