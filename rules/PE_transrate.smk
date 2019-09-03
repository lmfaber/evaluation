rule PE_transrate:
    input:
        R1 = lambda wildcards: config['samples'][wildcards.sample]['reads']['left'],
        R2 = lambda wildcards: config['samples'][wildcards.sample]['reads']['right'],
        refTranscripts = lambda wildcards: config['samples'][wildcards.sample]['reference']['transcripts'],
        contigs = lambda wildcards: config['samples'][wildcards.sample]['paths'].split(', ')
    output: 'output/pe/{sample}/transrate/assemblies.csv'
    threads: config['threads']
    params:
        assemblies = lambda wildcards: config['samples'][wildcards.sample]['paths'].replace(' ', ''),
        outDir = 'output/pe/{sample}/transrate'
    conda: '../envs/transrate.yml'
    shell: 'transrate --left {input.R1} \
                      --right {input.R2} \
                      --reference {input.refTranscripts} \
                      --threads {threads} \
                      --output {params.outDir} \
                      --assembly {params.assemblies}'