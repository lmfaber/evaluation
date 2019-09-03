rule PE_rnaquast:
    input:
        R1 = lambda wildcards: config['samples'][wildcards.sample]['reads']['left'],
        R2 = lambda wildcards: config['samples'][wildcards.sample]['reads']['right'],
        refGenome = lambda wildcards: config['samples'][wildcards.sample]['reference']['genome'],
        refGTF = lambda wildcards: config['samples'][wildcards.sample]['reference']['gtf'],
    output: sensitivity = 'output/pe/{sample}/rnaquast/comparison_output/sensitivity.txt'
    threads: config['threads']
    params:
        names = lambda wildcards: gen_names(wildcards),
        paths = lambda wildcards: config['samples'][wildcards.sample]['paths'].replace(',', ''),
        outDir = 'output/pe/{sample}/rnaquast',
        prokaryote = lambda wildcards: '--prokaryote' if config['samples'][wildcards.sample]['prokaryote'] else ''
    conda: '../envs/rnaquast.yml'
    shell: 'rnaQUAST.py -1 {input.R1} \
						-2 {input.R2} \
						-r {input.refGenome} \
						--gtf {input.refGTF} \
						-o {params.outDir} \
						-t {threads} \
						-c {params.paths} \
						-l {params.names} \
                        {params.prokaryote} \
                        --debug \
						--gene_mark \
						--disable_infer_genes \
						--disable_infer_transcripts'