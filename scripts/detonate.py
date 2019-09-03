# -*- coding: utf-8 -*-
import argparse
import os
from os.path import splitext, exists, abspath, basename
import sys
import logging
import subprocess

def basename_without_ext(path):
    """Accepts a path and returns the basename without extension.
    Arguments:
        path {[type]} -- path of a file
    Returns:
        [type] -- basename without extension of the given file
    """
    return splitext(basename(path))[0]


################
## DETONATE
################
## 1) run RSEM-EVAL on all assemblies, this calculates scores for each assembly based on the provided reads (so without reference)
##
## 2) after this score, we can additionally run REF-EVAL with a reference transcriptome
# command:  sh scripts/detonate.sh eco_SRR1173967.ds 10 /mnt/dessertlocal/projects/transcriptome_assembly/data/genomes/eco/eco.cdna.fa

# python3 detonate.py -t 40 -o ECO-TEST -L 99 -1 /mnt/prostlocal/lasse/snakemake_pipelines/multiAssembly/output/eco_ERR779269/fastp/eco_ERR779269_1.fastq -2 /mnt/prostlocal/lasse/snakemake_pipelines/multiAssembly/output/eco_ERR779269/fastp/eco_ERR779269_2.fastq -A /mnt/prostlocal/lasse/snakemake_pipelines/merger/output/ecoSample.PE/grouper.fasta --transcripts /mnt/prostlocal/lasse/transcripto/reference_genomes/eco_str_k_12_substr_mg1655.ASM584v2.cdna.all.fa
# python3 detonate.py -t 40 -o ECO-TEST -L 75 -1 /mnt/prostlocal/lasse/snakemake_pipelines/multiAssembly/output/eco_SRR1793813/fastp/eco_SRR1793813_1.fastq -2 /mnt/prostlocal/lasse/snakemake_pipelines/multiAssembly/output/eco_SRR1793813/fastp/eco_SRR1793813_2.fastq -A /mnt/prostlocal/lasse/snakemake_pipelines/merger/output/eco_SRR1793813.PE/grouper.fasta --transcripts /mnt/prostlocal/lasse/transcripto/reference_genomes/eco_str_k_12_substr_mg1655.ASM584v2.cdna.all.fa

# python3 detonate.py -t 20 -o ECO-TEST -L 75 -1 /mnt/prostlocal/lasse/snakemake_pipelines/multiAssembly/output/eco_SRR1793813/fastp/eco_SRR1793813_1.fastq -2 /mnt/prostlocal/lasse/snakemake_pipelines/multiAssembly/output/eco_SRR1793813/fastp/eco_SRR1793813_2.fastq -A /mnt/prostlocal/lasse/snakemake_pipelines/multiAssembly/output/eco_SRR1793813/PE_assemblies/final/rnaspades_solo.fasta --transcripts /mnt/prostlocal/lasse/transcripto/reference_genomes/eco_str_k_12_substr_mg1655.ASM584v2.cdna.all.fa
    
# eco_SRR1793813_ds_1.fastq / ANDERER HEADER / DS
# python3 detonate.py -t 20 -o test_mmu_ERR1041328_ds -L 99 -1 /mnt/prostlocal/lasse/transcripto/reads/test_mmu_ERR1041328_ds_1.fastq -2 /mnt/prostlocal/lasse/transcripto/reads/test_mmu_ERR1041328_ds_2.fastq -A /mnt/prostlocal/lasse/snakemake_pipelines/multiAssembly/output/test_mmu_ERR1041328_ds/PE_assemblies/rnaspades/transcripts.fasta --transcripts /mnt/prostlocal/lasse/transcripto/reference_genomes/mmu_transcripts_gencode_vM21.fa


# python3 detonate.py -t 20 -o ECO-TEST -L 145 -1 /mnt/prostlocal/lasse/snakemake_pipelines/multiAssembly/output/test_eco_ERR2511955_ds/fastp/test_eco_ERR2511955_ds_1.fastq -2 /mnt/prostlocal/lasse/snakemake_pipelines/multiAssembly/output/test_eco_ERR2511955_ds/fastp/test_eco_ERR2511955_ds_2.fastq -A /mnt/prostlocal/lasse/snakemake_pipelines/multiAssembly/output/test_eco_ERR2511955_ds/PE_assemblies/trinity/Trinity.fasta --transcripts /mnt/prostlocal/lasse/transcripto/reference_genomes/eco_str_k_12_substr_mg1655.ASM584v2.cdna.all.fa

# python3 detonate.py -t 60 -o eco_ERR779269 -1 /mnt/prostlocal/lasse/snakemake_pipelines/multiAssembly/output/eco_ERR779269/fastp/eco_ERR779269_1.fastq -2 /mnt/prostlocal/lasse/snakemake_pipelines/multiAssembly/output/eco_ERR779269/fastp/eco_ERR779269_2.fastq -A /mnt/prostlocal/lasse/snakemake_pipelines/multiAssembly/output/eco_ERR779269/PE_assemblies/final/trinity.fasta --transcripts /mnt/prostlocal/lasse/transcripto/reference_genomes/eco_str_k_12_substr_mg1655.ASM584v2.cdna.all.fa

# python3 detonate.py -t 60 -o eco_ERR779269_ds -1 /mnt/prostlocal/lasse/snakemake_pipelines/multiAssembly/output/eco_ERR779269_ds/fastp/eco_ERR779269_ds_1.fastq -2 /mnt/prostlocal/lasse/snakemake_pipelines/multiAssembly/output/eco_ERR779269_ds/fastp/eco_ERR779269_ds_2.fastq -A /mnt/prostlocal/lasse/snakemake_pipelines/multiAssembly/output/eco_ERR779269_ds/PE_assemblies/final/trinity.fasta --transcripts /mnt/prostlocal/lasse/transcripto/reference_genomes/eco_str_k_12_substr_mg1655.ASM584v2.cdna.all.fa

parser = argparse.ArgumentParser(description='Detonate assembly evaluation pipeline.\n python3 detonate.py -t 40 -o ECO-TEST -L 99 -1 /mnt/prostlocal/lasse/snakemake_pipelines/multiAssembly/output/ecoSample/fastp/ecoSample_1.fastq -2 /mnt/prostlocal/lasse/snakemake_pipelines/multiAssembly/output/ecoSample/fastp/ecoSample_2.fastq -A /mnt/prostlocal/lasse/snakemake_pipelines/merger/output/ecoSample.PE/grouper.fasta --transcripts /mnt/prostlocal/lasse/transcripto/reference_genomes/eco_str_k_12_substr_mg1655.ASM584v2.cdna.all.fa')
parser.add_argument('-t', '--threads', dest='threads', metavar='INT', default=6, help='Threads to use. (default: 6)')
parser.add_argument('--transcripts', required=True, dest='refTranscripts', help='Reference transcripts.')
parser.add_argument('-o', '--output', required=True, dest='outDir', help='Output directory.')
#parser.add_argument('-L', required=False, dest='readLength', metavar='MAX_READ_LEN', help='Gets calculated on the way. For single-end data, L represents the average read length. For paired-end data, L represents the average fragment length. It should be a positive integer (real value will be rounded to the nearest integer)')
parser.add_argument('-q', '--quiet', dest='quiet', action='store_true', help='shut up DETONATE!')
parser.add_argument('-1', '--leftReads', dest='leftReads', metavar='RIGHT_READS', help='Left reads.')
parser.add_argument('-2', '--rightReads', dest='rightReads', metavar='LEFT_READS', help='Right reads.')
parser.add_argument('-s', dest='reads', metavar='READS', help='single end reads')
parser.add_argument('-A', '--assemblies', required=True, dest='assemblies', help='Comma separated list of assemblies.')
args = parser.parse_args()


# Create logger
FORMAT = '{asctime} [{levelname}]: {message}'
logging.basicConfig(level=logging.INFO, format=FORMAT, style='{', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)

# Don't allow single and paired end reads input at the same time
if args.reads and args.leftReads and args.rightReads:
    logger.info('Use -s for single end assemblies/reads\nor\nuse -1/-2 for paired end assemblies/reads')
    logger.info('Exiting...')
    sys.exit(1)

# Quiet option for detonate if the user wishes so
silencio = '-q' if args.quiet else ''

assemblies = [a.strip() for a in args.assemblies.split(',')]
logger.info('Creating output folder.')
if not exists(args.outDir):
    os.mkdir(args.outDir)

# estimate transcript length distributions
logger.info('RSEM-EVAL-ESTIMATE-TRANSCRIPT-LENGHT-DISTRIBUTION')
os.system(f"rsem-eval-estimate-transcript-length-distribution {args.refTranscripts} {args.outDir}/transcript_length_parameters.txt")

if args.reads:
    logger.info('SINGLE END ASSEMBLY')

    ##########################
    ## ESTIMATE THE TRUE ASSEMBLY
    # construct estimate of "true" assembly with RSEM
    logger.info('RSEM-PREPARE-REFERENCE')
    os.system(f'rsem-prepare-reference --bowtie {args.refTranscripts} {args.outDir}/rsem_ref')
    logger.info('RSEM-CALCULATE-EXPRESSION')
    os.system('rsem-calculate-expression -p {args.threads} {args.reads} {args.outDir}/rsem_ref {args.outDir}/rsem_expr')

    # command for Samtools 0.18 from the conda detonate package doesnt support -T and --threads
    logger.info('SAMTOOLS SORT')
    os.system(f'samtools sort {args.outDir}/rsem_expr.transcript.bam {args.outDir}/rsem_expr.transcript.sorted')

    # Now we are ready to estimate the “true” assembly:
    logger.info('REF-EVAL-ESTIMATE-TRUE-ASSEMBLY')
    os.system(f'ref-eval-estimate-true-assembly --reference {args.outDir}/rsem_ref --expression {args.outDir}/rsem_expr --assembly {args.outDir}/ta --alignment-policy best')

    ##########################
    ## COMPUTE KMER-COMPRESSION SCORE FOR EACH ASSEMBLY

    # This time, we will run RSEM to estimate the expression levels of each sequence in the estimated “true” assembly, as follows
    logger.info('RSEM-PREPARE-REFERENCE-2')
    os.system(f'rsem-prepare-reference --bowtie {args.outDir}/ta_0.fa {args.outDir}/ta_0_ref')
    logger.info('RSEM-CALCULATE-EXPRESSION-2')
    os.system(f'rsem-calculate-expression -p {args.threads} {args.reads} {args.outDir}/ta_0_ref {args.outDir}/ta_0_expr')

    # Estimate the num-reads parameter, this is the number of reads times the read length equals the total number of nucleotides in the read set.
    logger.info('Calculating number of reads.')
    output = os.popen("awk '{s++}END{print s/4}' " + args.leftReads)
    numReads = str(output.read())
    logger.info(f"Number of reads: {numReads}")

    # Get read length and insert size info via the mapped bam file from bowtie
    logger.info('Calculating read length and insert size.')
    p = subprocess.run(f"samtools stats {args.outDir}/rsem_expr.transcript.sorted.bam | grep -e '^SN\tmaximum length' | cut -f 3", stdout=subprocess.PIPE, shell=True)
    readLength = str(p.stdout.decode('utf-8'))
    logger.info(f"Read length: {readLength}")

    for assembly in assemblies:
        assembly = abspath(assembly)  # absolute path
        tool = basename_without_ext(assembly)  # basename without ext

        logger.info(f'RSEM-EVAL-CALCULATE-SCORE: {tool}')
        os.system(f'rsem-eval-calculate-score -p {args.threads} --transcript-length-parameters {args.outDir}/transcript_length_parameters.txt {args.reads} {assembly} {args.outDir}/rsem_eval_{tool} {readLength}')

        # We now compute the KC score of each assembly as follows:
        logger.info(f'REF-EVAL-KC (KC score): {tool}')

        output = os.popen(f'ref-eval --scores kc --A-seqs {assembly} --B-seqs {args.outDir}/ta_0.fa --B-expr {args.outDir}/ta_0_expr.isoforms.results --kmerlen {args.readLength} --readlen {args.readLength} --num-reads {numReads}')

        kcFile = f'{args.outDir}/kc_{tool}.txt'
        with open(kcFile, "w+") as writer:
            writer.write(output.read())

        # COMPUTE ALIGNMENT-BASED SCORES FOR EACH ASSEMBLY WITH BLAT
        logger.info(f'BLAT (Alignment based scores): {tool}')
        os.system(f'blat -minIdentity=80 {args.outDir}/ta_0.fa {assembly} {args.outDir}/{tool}_to_ta_0.psl')
        os.system(f'blat -minIdentity=80 {assembly} {args.outDir}/ta_0.fa {args.outDir}/ta_0_to_{tool}.psl')

        # We can now compute the contig and nucleotide scores as follows:
        logger.info(f'REF-EVAL (Contig/Nucleotide scores): {tool}')
        output = os.popen(f'ref-eval --scores contig,nucl --weighted no --A-seqs {assembly} --B-seqs {args.outDir}/ta_0.fa --A-to-B {args.outDir}/{tool}_to_ta_0.psl --B-to-A {args.outDir}/ta_0_to_{tool}.psl --min-frac-identity 0.90')

        contigNuclFile = f'{args.outDir}/contig_nucl_{tool}.txt'
        with open(contigNuclFile, "w+") as writer:
            writer.write(output.read())

else:
    logger.info('PAIRED END ASSEMBLY')

    ##########################
    ## ESTIMATE THE TRUE ASSEMBLY
    # construct estimate of "true" assembly with RSEM
    logger.info('RSEM-PREPARE-REFERENCE')
    os.system(f"rsem-prepare-reference {silencio} --bowtie {args.refTranscripts} {args.outDir}/rsem_ref")
    
    logger.info('RSEM-CALCULATE-EXPRESSION')
    os.system(f"rsem-calculate-expression {silencio} -p {args.threads} --paired-end {args.leftReads} {args.rightReads} {args.outDir}/rsem_ref {args.outDir}/rsem_expr")

    # command for Samtools 0.18 from the conda detonate package doesnt support -T and --threads
    logger.info('SAMTOOLS SORT')
    os.system(f"samtools sort {args.outDir}/rsem_expr.transcript.bam -o {args.outDir}/rsem_expr.transcript.sorted.bam -T {args.outDir}/tmp --threads {args.threads}")

    # Now we are ready to estimate the “true” assembly:
    logger.info('REF-EVAL-ESTIMATE-TRUE-ASSEMBLY')

    os.system(f"ref-eval-estimate-true-assembly --reference {args.outDir}/rsem_ref --expression {args.outDir}/rsem_expr --assembly {args.outDir}/ta --alignment-policy best")

    ##########################
    ## COMPUTE KMER-COMPRESSION SCORE FOR EACH ASSEMBLY
    # This time, we will run RSEM to estimate the expression levels of each sequence in the estimated “true” assembly, as follows
    logger.info('RSEM-PREPARE-REFERENCE-2')
    os.system(f"rsem-prepare-reference {silencio} --bowtie {args.outDir}/ta_0.fa {args.outDir}/ta_0_ref")

    logger.info('RSEM-CALCULATE-EXPRESSION-2')
    os.system(f"rsem-calculate-expression {silencio} -p {args.threads} --paired-end {args.leftReads} {args.rightReads} {args.outDir}/ta_0_ref {args.outDir}/ta_0_expr")

    # Estimate the num-reads parameter, this is the number of reads times the read length equals the total number of nucleotides in the read set.
    logger.info('Calculating number of reads.')
    output = os.popen("awk '{s++}END{print s/4}' " + args.leftReads)
    numReads = str(output.read())
    logger.info(f"Number of reads: {numReads}")

    # Get read length and insert size info via the mapped bam file from bowtie
    logger.info('Calculating read length and insert size.')
    p = subprocess.run(f"samtools stats {args.outDir}/rsem_expr.transcript.sorted.bam | grep -e '^SN\tmaximum length' -e '^SN\tinsert size average' | cut -f 3", stdout=subprocess.PIPE, shell=True)
    output = p.stdout.decode('utf-8').split()
    readLength = int(round(float(output[0])))
    insertSize = float(output[1])
    # for paired end reads: fragment size = 2xreadLength + insertSize
    fragmentLength = (2 * readLength) + insertSize
    logger.info(f"Read length: {readLength}, insert size: {insertSize}, fragment length: {fragmentLength}")


####### BLAT36 and DETONATE

    for assembly in assemblies:
        assembly = abspath(assembly)  # absolute path
        tool = basename_without_ext(assembly)  # basename without ext
        
        args.readLength = 400
        # RUN RSEM_EVAL ON EACH ASSEMBLY
        logger.info(f'RSEM-EVAL-CALCULATE-SCORE: {tool}')
        os.system(f"rsem-eval-calculate-score {silencio} --paired-end {args.leftReads} {args.rightReads} {assembly} {args.outDir}/rsem_eval_{tool} {fragmentLength} --transcript-length-parameters {args.outDir}/transcript_length_parameters.txt -p {args.threads}")


        # We now compute the KC score of each assembly as follows:
        logger.info(f'REF-EVAL-KC (KC score): {tool}')
        output = os.popen(f"ref-eval --scores kc --A-seqs {assembly} --B-seqs {args.outDir}/ta_0.fa --B-expr {args.outDir}/ta_0_expr.isoforms.results --kmerlen {readLength} --readlen {readLength} --num-reads {numReads}")

        kcFile = f"{args.outDir}/kc_{tool}.txt"
        with open(kcFile, "w+") as writer:
            writer.write(output.read())

        # COMPUTE ALIGNMENT-BASED SCORES FOR EACH ASSEMBLY WITH BLAT
        logger.info(f'BLAT (Alignment based scores): {tool}')
        #for assembly in assemblies:
        os.system(f"blat -minIdentity=80 {args.outDir}/ta_0.fa {assembly} {args.outDir}/{tool}_to_ta_0.psl")
        os.system(f"blat -minIdentity=80 {assembly} {args.outDir}/ta_0.fa {args.outDir}/ta_0_to_{tool}.psl")

        # We can now compute the contig and nucleotide scores as follows:
        logger.info(f'REF-EVAL (Contig/Nucleotide scores): {tool}')
        output = os.popen(f"ref-eval --scores contig,nucl --weighted no --A-seqs {assembly} --B-seqs {args.outDir}/ta_0.fa --A-to-B {args.outDir}/{tool}_to_ta_0.psl --B-to-A {args.outDir}/ta_0_to_{tool}.psl --min-frac-identity 0.90")
        contigNuclFile = f"{args.outDir}/contig_nucl_{tool}.txt"

        with open(contigNuclFile, "w+") as writer:
            writer.write(output.read())

