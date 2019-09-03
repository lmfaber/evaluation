rule SE_rnaquast:
    input:
        R1 = lambda wildcards: config['samples'][wildcards.sample]['reads'],
        refGenome = lambda wildcards: config['samples'][wildcards.sample]['reference']['genome'],
        refGTF = lambda wildcards: config['samples'][wildcards.sample]['reference']['gtf'],
    output: sensitivity = 'output/se/{sample}/rnaquast/comparison_output/sensitivity.txt'
    threads: config['threads']
    params:
        names = lambda wildcards: gen_names(wildcards),
        paths = lambda wildcards: config['samples'][wildcards.sample]['paths'].replace(',', ''),
        outDir = 'output/se/{sample}/rnaquast',
        prokaryote = lambda wildcards: '--prokaryote' if config['samples'][wildcards.sample]['prokaryote'] else ''
    conda: '../envs/rnaquast.yml'
    shell: 'rnaQUAST.py -s {input.R1} \
						-r {input.refGenome} \
						--gtf {input.refGTF} \
						-o {params.outDir} \
						-t {threads} \
						-c {params.paths} \
						-l {params.names} \
                        {params.prokaryote} \
						--gene_mark\
						--disable_infer_genes\
						--disable_infer_transcripts'