import os, glob
import pandas as pd
from Bio import SeqIO

R1 = str(snakemake.input['R1'])
R2 = str(snakemake.input['R2'])
THREADS = str(snakemake.threads)
assemblies = snakemake.input['assemblies']
hisat2_dir = str(snakemake.params['out_dir'])
output_file = str(snakemake.output)

if not os.path.exists(hisat2_dir):
    os.mkdir(hisat2_dir)

# Index
indices = []
number_of_bases = {}

for assembly_path in assemblies:
    # Create mapping index
    base_name = os.path.splitext(os.path.basename(assembly_path))[0]
    indices.append(base_name)
    index_prefix = f'{hisat2_dir}/{base_name}'
    os.system(f'hisat2-build -q {assembly_path} {index_prefix}')

    # Count number of bases
    with open(assembly_path, 'r') as fastaReader:
        # original_fasta_sequences = SeqIO.to_dict(SeqIO.parse(fastaReader, 'fasta'))
        n_bases = 0
        for record in SeqIO.parse(fastaReader, 'fasta'):
            n_bases += len(record.seq)
    number_of_bases[base_name] = n_bases

mapping_infos = pd.DataFrame(index = ['mapping_rate', 'scaling_factor', 'mapping_rate_per_base'], columns = indices)
mapping_infos = mapping_infos.fillna(0.0)

dummy_sam_file = f'{hisat2_dir}/tmp.sam'
combined_metric = {}

for index in indices:
    
    metric_filename = f'{hisat2_dir}/{index}.txt'
    index_file = f'{hisat2_dir}/{index}'

    # Mapping
    if not os.path.exists(metric_filename):
        os.system(f'hisat2 -x {index_file} -1 {R1} -2 {R2} --threads {THREADS} -S {dummy_sam_file} &> {metric_filename}')

    # Read overall mapping rate
    with open(metric_filename) as reader:
        last_line = reader.readlines()[-1]
        percentage_mapped = float(last_line.split('%')[0])
        mapping_infos[index]['mapping_rate'] = percentage_mapped
        combined_metric[index] = percentage_mapped / number_of_bases[index]

    # Remove index
    for f in glob.glob(f"{hisat2_dir}/{index}.*.ht2"):
        os.remove(f)


# Mapping rate in combination with number of bases
scaling_factor = 0
while not any(x >= 0.1 for x in [a for a in combined_metric.values()]):
    for key, value in combined_metric.items():
        combined_metric[key] = value*10
    scaling_factor += 1

print(f'Scaling factor is {scaling_factor}.')

for index in indices:
    mapping_infos[index]['scaling_factor'] = scaling_factor
    mapping_infos[index]['mapping_rate_per_base'] = combined_metric[index]

mapping_infos.to_csv(output_file, sep = '\t', header = True, index = True)

if os.path.exists(dummy_sam_file):
    os.remove(dummy_sam_file)
