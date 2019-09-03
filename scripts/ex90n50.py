import os
from shutil import copyfile
import pandas as pd

inputFiles = snakemake.input['fasta']
leftReads = str(snakemake.input['R1'])
rightReads = str(snakemake.input['R2'])
base_output_folder = str(snakemake.output)
threads = str(snakemake.threads)

if not os.path.exists(base_output_folder):
    os.mkdir(base_output_folder)
import time
for inputFile in inputFiles:
    base_name = os.path.splitext(os.path.basename(inputFile))[0]
    output_dir = f'{base_output_folder}/{base_name}'
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

os.chdir(base_output_folder)

for inputFile in inputFiles:
    base_name = os.path.splitext(os.path.basename(inputFile))[0]
    input_fasta = f'{base_name}.fasta'

    os.chdir(base_name)

    copyfile(inputFile, input_fasta)
    os.system(f'align_and_estimate_abundance.pl --thread_count {threads} --transcripts {input_fasta} --est_method salmon '
              f'--trinity_mode --prep_reference')
    abundance_dir = f'salmon'
    os.system(f'align_and_estimate_abundance.pl --thread_count {threads} --transcripts {input_fasta} --seqType fq --left {leftReads} '
              f'--right {rightReads} --est_method salmon --trinity_mode --output_dir {abundance_dir}')
    os.system(f'abundance_estimates_to_matrix.pl --est_method salmon --gene_trans_map none --out_prefix salmon --name_sample_by_basedir {abundance_dir}/quant.sf')
    os.system(f'contig_ExN50_statistic.pl salmon.isoform.TPM.not_cross_norm {input_fasta} | tee ExN50.stats')

    os.chdir('..')

