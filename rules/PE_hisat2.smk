rule PE_mapping:
    input:
        R1 = lambda wildcards: config['samples'][wildcards.sample]['reads']['left'],
        R2 = lambda wildcards: config['samples'][wildcards.sample]['reads']['right'],
        assemblies = lambda wildcards: config['samples'][wildcards.sample]['paths'].split(', ')
    output: 'output/pe/{sample}/hisat2/mapping_rates.csv'
    params: 
        out_dir = 'output/pe/{sample}/hisat2'
    threads: config['threads']
    conda: '../envs/hisat2.yml'
    script: '../scripts/hisat2.py'