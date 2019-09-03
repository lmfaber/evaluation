rule PE_detonate:
	input:
		R1 = lambda wildcards: config['samples'][wildcards.sample]['reads']['left'],
        R2 = lambda wildcards: config['samples'][wildcards.sample]['reads']['right'],
		refTranscripts = lambda wildcards: config['samples'][wildcards.sample]['reference']['transcripts'],
	output:
		outDir = directory('output/pe/{sample}/detonate')
	threads: config['threads']
	params:
		assemblies = lambda wildcards: config['samples'][wildcards.sample]['paths'].replace(' ', ''),
	# conda: '../envs/detonate.yml'
	shell: 'python3 scripts/detonate.py --threads {threads} \
										   --transcripts {input.refTranscripts} \
										   --output {output.outDir} \
										   --quiet \
										   --leftReads {input.R1} \
										   --rightReads {input.R2} \
										   --assemblies {params.assemblies}'


# import json, os
# def numberOfReads(wildcards):
# 	fastpFile = config['samples'][wildcards.sample]['fastp']
# 	with open(fastpFile, 'r') as jsonFile:
# 		data = json.load(jsonFile)
# 		numberOfReads = data['summary']['after_filtering']['total_reads']
# 	return(numberOfReads)

# def readLength(wildcards):
# 	fastpFile = config['samples'][wildcards.sample]['fastp']
# 	with open(fastpFile, 'r') as jsonFile:
# 		data = json.load(jsonFile)
# 		readLength = str(data['summary']['after_filtering']['read1_mean_length'])
# 	return(readLength)

# def fragmentLength(wildcards):
# 	fastpFile = config['samples'][wildcards.sample]['fastp']
# 	with open(fastpFile, 'r') as jsonFile:
# 		data = json.load(jsonFile)
# 		readLength = str(data['summary']['after_filtering']['read1_mean_length'])
# 		insertSize = data['insert_size']['peak']
# 	fragmentLength = str((2 * int(readLength)) + int(insertSize))
# 	return(fragmentLength)


# rule rsem_prepare_reference:
# 	input: 
# 		referenceTranscripts = lambda wildcards: config['samples'][wildcards.sample]['reference']['transcripts'],
# 		R1 = lambda wildcards: config['samples'][wildcards.sample]['reads']['left'],
# 		R2 = lambda wildcards: config['samples'][wildcards.sample]['reads']['right']
# 	output: 'output/pe/{sample}/rsem_expr.transcript.bam'
# 	params: dir = 'output/pe/{sample}'
# 	threads: config['threads']
# 	# conda: '../envs/rsem.yml'
# 	shell: '''
# 	rsem-prepare-reference -q --bowtie {input.referenceTranscripts} {params.dir}/rsem_ref
# 	rsem-calculate-expression -q -p {threads} --paired-end {input.R1} {input.R2} {params.dir}/rsem_ref {params.dir}/rsem_expr
# 	'''

# rule samtools_sort:
# 	input: rules.rsem_prepare_reference.output
# 	output:
# 		sortedBam = 'output/pe/{sample}/rsem_expr.transcript.sorted.bam'
# 	params: tempDir = 'output/pe/{sample}/tmp'
# 	threads: config['threads']
# 	# conda: 'envs/rsem.yml'
# 	shell: 'samtools sort {input} -o {output.sortedBam} -T {params.tempDir} --threads {threads}'

# rule ref_eval_estimate_true_assembly:
# 	input: rules.samtools_sort.output.sortedBam
# 	output: 'output/pe/{sample}/ta_0.fa'
# 	params: dir = 'output/pe/{sample}'
# 	# conda: 'envs/detonate.yml'
# 	shell: 'ref-eval-estimate-true-assembly --reference {params.dir}/rsem_ref --expression {params.dir}/rsem_expr --assembly {params.dir}/ta --alignment-policy best'

# ##########################
# ## COMPUTE KMER-COMPRESSION SCORE FOR EACH ASSEMBLY
# # This time, we will run RSEM to estimate the expression levels of each sequence in the estimated “true” assembly, as follows
# rule rsem_prepare_reference_2:
# 	input: 
# 		TA = rules.ref_eval_estimate_true_assembly.output,
# 		R1 = lambda wildcards: config['samples'][wildcards.sample]['reads']['left'],
# 		R2 = lambda wildcards: config['samples'][wildcards.sample]['reads']['right']
# 	output: 'output/pe/{sample}/ta_0_ref.transcripts.fa'
# 	params: dir = 'output/pe/{sample}'
# 	threads: config['threads']
# 	# conda: 'envs/rsem.yml'
# 	shell: '''
# 	rsem-prepare-reference -q --bowtie {input.TA} {params.dir}/ta_0_ref
# 	rsem-calculate-expression -q -p {threads} --paired-end {input.R1} {input.R2} {params.dir}/ta_0_ref {params.dir}/ta_0_expr
# 	'''

# ############## AB HIER BLAT36 UND DETONATE
# rule transcript_length_distribution:
# 	input: 
# 		referenceTranscripts = lambda wildcards: config['samples'][wildcards.sample]['reference']['transcripts'],
# 		rsem = rules.rsem_prepare_reference_2.output
# 	output: 'output/pe/{sample}/transcript_length_parameters.txt'
# 	# conda: 'envs/detonate.yml'
# 	shell: 'rsem-eval-estimate-transcript-length-distribution {input.referenceTranscripts} {output}'



# rule PE_detonate:
# 	input: 
# 		assembly = lambda wildcards: config['samples'][wildcards.sample]['paths'].split(', '),
# 		R1 = lambda wildcards: config['samples'][wildcards.sample]['reads']['left'],
# 		R2 = lambda wildcards: config['samples'][wildcards.sample]['reads']['right'],
# 		referenceTranscripts = lambda wildcards: config['samples'][wildcards.sample]['reference']['transcripts'],
# 		rsem = rules.rsem_prepare_reference_2.output,
# 		transcript_length_parameters = rules.transcript_length_distribution.output
# 	output: 
# 		kc = 'output/pe/{sample}/kc_{tool}.txt',
# 		nucl = 'output/pe/{sample}/contig_nucl_{tool}.txt',
# 		score = 'output/pe/{sample}/rsem_eval_{tool}.score'
# 	params:
# 		dir = 'output/pe/{sample}',
# 		tool = '{tool}',
# 		fragmentLength = lambda wildcards: fragmentLength(wildcards),
# 		readLength = lambda wildcards: readLength(wildcards),
# 		numReads = lambda wildcards: numberOfReads(wildcards)
# 	threads: config['threads']
# 	# conda: '../envs/detonate.yml'
# 	shell: '''	
# 	# RUN RSEM_EVAL ON EACH ASSEMBLY
# 	rsem-eval-calculate-score -q --paired-end {input.R1} {input.R2} {input.assembly} {params.dir}/rsem_eval_{params.tool} {params.fragmentLength} --transcript-length-parameters {input.transcript_length_parameters} -p {threads}

# 	# We now compute the KC score of each assembly as follows:
# 	ref-eval --scores kc --A-seqs {input.assembly} --B-seqs {params.dir}/ta_0.fa --B-expr {params.dir}/ta_0_expr.isoforms.results --kmerlen {params.readLength} --readlen {params.readLength} --num-reads {params.numReads} | tee {params.dir}/kc_{params.tool}.txt

# 	# COMPUTE ALIGNMENT-BASED SCORES FOR EACH ASSEMBLY WITH BLAT
# 	blat -minIdentity=80 {params.dir}/ta_0.fa {input.assembly} {params.dir}/{params.tool}_to_ta_0.psl
# 	blat -minIdentity=80 {input.assembly} {params.dir}/ta_0.fa {params.dir}/ta_0_to_{params.tool}.psl

# 	# We can now compute the contig and nucleotide scores as follows:
# 	ref-eval --scores contig,nucl --weighted no --A-seqs {input.assembly} --B-seqs {params.dir}/ta_0.fa --A-to-B {params.dir}/{params.tool}_to_ta_0.psl --B-to-A {params.dir}/ta_0_to_{params.tool}.psl --min-frac-identity 0.90 | tee {params.dir}/contig_nucl_{params.tool}.txt
# 	'''
