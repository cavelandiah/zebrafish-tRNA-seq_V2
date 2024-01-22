
rule build_segemehl_index:
    input:
        refs = 'resources/references/{ref_set}_refs.fa'
    output:
        index = 'resources/index/{ref_set}.idx'
    shell:
        "segemehl.x -x {output.index} -d {input.refs}"

rule build_segemehl_BS_index:
    input:
        refs = 'resources/references/{ref_set}_refs.fa'
    output:
        index_ct = 'resources/index/{ref_set}.ctidx',
        index_ga = 'resources/index/{ref_set}.gaidx'
    shell:
        "segemehl.x -x {output.index_ct} -y {output.index_ga} -F 1 -d {input.refs}"


rule mapping_segemehl:
    input:
        index =  'resources/index/{ref_set}.idx',
        index_ct = 'resources/index/{ref_set}.ctidx',
        index_ga = 'resources/index/{ref_set}.gaidx',
        fastqT = 'resources/filteres-reads/a_u_t_q_{sample}.fastq',
        #fastqF = 'resources/filteres-reads/a_u_f_q_{sample}.fastq',
        fastq = 'resources/filteres-reads/a_u_q_{sample}.fastq',
        fasta_files = ['resources/references/{ref_set}_refs.fa'],
    output:
        sam = 'resources/mapping/{ref_set}/polyTtrimming_{polyT}/{sample}.sam',
        unmapped_fastq = 'resources/mapping/{ref_set}/polyTtrimming_{polyT}/{sample}_unmapped.fastq',
        mapping_stats = "qc/mapping/{ref_set}_polyTtrimming_{polyT}_{sample}.tsv"
    #resources:
        #mem_mb = 500
    threads: 2
    run:
        import subprocess
        import pandas as pd

        def get_mapping_stats(data, filepath, sample_name):
            """Get statistics on mapping from segemehl stderr output."""
            name = sample_name
            for i, line in enumerate(data.split('\n')):
                if 'Mapping stats' in line:
                    mapping_info_start = i
            mapping_stats = data.split('\n')[mapping_info_start+1:mapping_info_start+3]
            mapping_stats = [m.split('\t') for m in mapping_stats]
            mapping_stats[1][0] = name
            mapping_stats[0][0] = name
            with open(filepath, 'w') as f:
                for line in mapping_stats:
                    f.write('\t'.join(line)+'\n')
            return mapping_stats

        # select input fastq files with correct poly T treatment
        if wildcards.polyT == "OFF":
            fastq = input.fastq
        elif wildcards.polyT == "ONtttt":
            fastq = input.fastqT
        elif wildcards.polyT == "ONfixed":
            fastq = input.fastqF
        else:
            print('no valide polyT wildcard')

        # select parameter set (see config)
        treatment = sample_dict[wildcards.sample]['treatment']
        if wildcards.ref_set == 'all':
            params = config['pre_mapping'][treatment]
        else:
            params = config['mapping'][treatment]


        # check if bisulfite mapping should be used
        if sample_dict[wildcards.sample]['treatment'] == 'BS':
            call_args = ['segemehl.x',
                         '-i', input.index_ct,
                         '-j',input.index_ga,
                         '-d' ]  + input.fasta_files + [
                         '-q', fastq,
                         '-A', params['accuracy'],
                         '-M', params['maxinterval'],
                         '-D', params['difference'],
                         '-E', params['evalue'],
                         '-t' , str(threads),
                         '-o', output.sam,
                         '-u', output.unmapped_fastq,
                         '-F', params['bisulfite'],
                        ]
        else:
            call_args = ['segemehl.x',
                         '-i', input.index,
                         '-d' ]  + input.fasta_files + [
                         '-q', fastq,
                         '-A', params['accuracy'],
                         '-M', params['maxinterval'],
                         '-D', params['difference'],
                         '-E', params['evalue'],
                         '-t' , str(threads),
                         '-o', output.sam,
                         '-u', output.unmapped_fastq,
                        ]


        results = subprocess.run(call_args, capture_output=True)
        data = results.stderr.decode()
        print(data)
        mapping_stat = get_mapping_stats(data, output.mapping_stats, wildcards.sample)
