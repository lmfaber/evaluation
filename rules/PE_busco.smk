rule PE_busco:
    input:
        fasta = lambda wildcards: config['samples'][wildcards.sample]['paths'].split(', '),
        database = lambda wildcards: config['samples'][wildcards.sample]['busco_database']
    output: directory('output/pe/{sample}/busco')
    threads: config['threads']
    params:
        sample = '{sample}'
    conda: '../envs/busco.yml'
    script: '../scripts/busco.py'