infile = config['runs']

runs = []
with open(infile, 'r') as f:
    for line in f:
        runs.append(line.strip('\n'))

rule all:
    input: expand('{run}/quant.sf', run=runs)

rule fetch_fastq:
    output: '{run}.fastq'
    shell: "/usr/bin/fasterq-dump {wildcards.run}"

rule run_salmon:
    input: '{run}.fastq'
    output:
        outdir = directory('{run}'),
        counts = '{run}/quant.sf'
    params:
        salmon_index = '/data/salmon_index/grch38_salmon_index/',
        genemap = '/data/salmon_index/tx_gene_map.tsv',
    shell:
        '''
        /opt/salmon-latest_linux_x86_64/bin/salmon quant \
          -i {params.salmon_index} \
          -l A \
          -r {input} \
          -o {output} \
          -g {params.genemap} \
          -p 24 \
          --validateMappings
        '''
