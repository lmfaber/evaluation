rule PE_ex90n50:
    input:
        fasta = lambda wildcards: config['samples'][wildcards.sample]['paths'].split(', '),
        R1 = lambda wildcards: config['samples'][wildcards.sample]['reads']['left'],
        R2 = lambda wildcards: config['samples'][wildcards.sample]['reads']['right']
    output: directory('output/pe/{sample}/EX90N50')
    conda: '../envs/ex90n50.yml'
    threads: config['threads']
    script: '../scripts/ex90n50.py'