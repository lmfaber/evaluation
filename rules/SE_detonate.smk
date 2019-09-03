rule SE_detonate:
	input:
		R1 = lambda wildcards: config['samples'][wildcards.sample]['reads'],
		refTranscripts = lambda wildcards: config['samples'][wildcards.sample]['reference']['transcripts'],
	output:
		outDir = directory('output/se/{sample}/detonate')
	threads: config['threads']
	params:
		assemblies = lambda wildcards: config['samples'][wildcards.sample]['paths'].replace(' ', ''),
	# conda: '../envs/detonate.yml'
	shell: 'python3 ../scripts/detonate.py --threads {threads} \
										   --transcripts {input.refTranscripts} \
										   --output {output.outDir} \
										   --quiet \
										   --s {input.R1} \
										   --assemblies {params.assemblies}'