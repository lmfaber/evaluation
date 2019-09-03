import os.path 
import json 

def gen_names(wildcards, sep = ' ', array=False):
	"""Generates a 'sep'-separated string with all input names
	
	Keyword Arguments:
		sep {str} -- Separator (default: {' '})
	"""
	singlePaths = config['samples'][wildcards.sample]['paths'].split(', ')
	if array:
		return([basename_without_ext(a) for a in singlePaths])
	else:
		return(sep.join(basename_without_ext(a) for a in singlePaths))
	
def basename_without_ext(path):
    """Accepts a path and returns the basename without extension.
    Arguments:
        path {[type]} -- path of a file
    Returns:
        [type] -- basename without extension of the given file
    """
    return os.path.splitext(os.path.basename(path))[0]

def gen_paths(wildcards, sep = ' ', array=False):
	"""
	Generates a 'sep'-separated string with all input paths. Either as array output or as string with sep.
	
	Keyword Arguments:
		sep {str} -- Separator (default: {' '})

	"""
	singlePaths = config['samples'][wildcards.sample]['paths'].split(', ')
	singlePaths = [os.path.abspath(a) for a in singlePaths]
	if array:
		return(singlePaths)
	else:
		return(sep.join(singlePaths))

SE_TOOLS, PE_TOOLS = {}, {}
SE_SAMPLES, PE_SAMPLES = [], []

for i, sample in enumerate(config['samples']):
	# Paired end data
	if 'left' in config['samples'][sample]['reads']:
		PE_SAMPLES.append(sample)
		paths = config['samples'][sample]['paths']
		basenames = [os.path.splitext(os.path.basename(A))[0] for A in paths.split(', ')]
		PE_TOOLS[sample] = basenames
	# Single end data
	elif 'left' not in config['samples'][sample]['reads']:
		SE_SAMPLES.append(sample)
		paths = config['samples'][sample]['paths']
		basenames = [os.path.splitext(os.path.basename(A))[0] for A in paths.split(', ')]
		SE_TOOLS[sample] = basenames
	else:
		print(f'Could not determine if paired or single end data for {sample}')

# print(SE_TOOLS)
# print(PE_TOOLS)
# print(SE_SAMPLES)
# print(PE_SAMPLES)


rule all:
	input:
		# Paired end output
		pe_metrics = expand('output/pe/{sample}/metrics.csv', sample=PE_SAMPLES),
		se_metrics = expand('output/se/{sample}/metrics.csv', sample=SE_SAMPLES),
		# Dummy detonate files
		# detonate = [expand('output/pe/{sample}/detonate/', sample=sample, tool=PE_TOOLS[sample])],

		# kc = [expand('output/pe/{sample}/kc_{tool}.txt', sample=sample, tool=PE_TOOLS[sample]) for sample in PE_SAMPLES],
		# nucl = [expand('output/pe/{sample}/contig_nucl_{tool}.txt', sample=sample, tool=PE_TOOLS[sample]) for sample in PE_SAMPLES],
		# score = [expand('output/pe/{sample}/rsem_eval_{tool}.score', sample=sample, tool=PE_TOOLS[sample]) for sample in PE_SAMPLES],
		# pe_transrate = [expand('output/pe/{sample}/transrate/assemblies.csv', sample=PE_SAMPLES)],
		# sensitivity = [expand('output/pe/{sample}/rnaquast/comparison_output/sensitivity.txt', sample=PE_SAMPLES)],

		# Single end output
		# SE_pre = [expand('output/se/{sample}/assemblies/{tool}.fasta', sample=sample, tool=SE_TOOLS[sample]) for sample in SE_SAMPLES],
		# SE_kc = [expand('output/se/{sample}/kc_{tool}.txt', sample=sample, tool=SE_TOOLS[sample]) for sample in SE_SAMPLES],
		# SE_nucl = [expand('output/se/{sample}/contig_nucl_{tool}.txt', sample=sample, tool=SE_TOOLS[sample]) for sample in SE_SAMPLES],
		# SE_score = [expand('output/se/{sample}/rsem_eval_{tool}.score', sample=sample, tool=SE_TOOLS[sample]) for sample in SE_SAMPLES],

		
include: 'rules/PE_detonate.smk'
include: 'rules/PE_transrate.smk'
include: 'rules/PE_rnaquast.smk'
include: 'rules/PE_busco.smk'
include: 'rules/PE_Ex90N50.smk'
include: 'rules/PE_hisat2.smk'

include: 'rules/SE_detonate.smk'
include: 'rules/SE_rnaquast.smk'

rule PE_extract_metrices:
	input:
		rnaquast_sensitivity = rules.PE_rnaquast.output.sensitivity,
		detonateDir = rules.PE_detonate.output,
		transrate = rules.PE_transrate.output,
		busco = rules.PE_busco.output,
		EX90N50 = rules.PE_ex90n50.output,
		mapping_rates = rules.PE_mapping.output
	output: 
		raw = 'output/pe/{sample}/metrics.csv',
		minRaw = 'output/pe/{sample}/reduced_metrics.csv',
		normalized = 'output/pe/{sample}/norm_metrics.csv',
		heatmap = 'output/pe/{sample}/heatmap.svg'
	params:
		assemblies = lambda wildcards: gen_names(wildcards, array=True)
	conda: 'envs/python3.yml'
	script: 'scripts/extract_metrices.py'

rule SE_extract_metrices:
		input:
			rnaquast_sensitivity = rules.SE_rnaquast.output.sensitivity,
			detonateDir = rules.SE_detonate.output
		output: 
			raw = 'output/se/{sample}/metrics.csv',
			minRaw = 'output/se/{sample}/reduced_metrics.csv',
			normalized = 'output/se/{sample}/norm_metrics.csv',
			heatmap = 'output/se/{sample}/heatmap.svg'
		conda: 'envs/python3.yml'
		script: 'scripts/extract_metrices.py'

