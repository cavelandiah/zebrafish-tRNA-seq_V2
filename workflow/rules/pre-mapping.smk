rule build_index:
    input:
        genomic = config['genomic_refs'],
        mitochondrial = config['mitochondrial_refs'],
        #housekeeping = config['housekeeping_refs']
    output:
        index = config['index']
    shell:
        "segemehl.x -x {output.index} -d {input.genomic} {input.mitochondrial}" # {input.housekeeping}"

rule build_BS_index:
    input:
        genomic = config['genomic_refs'],
        mitochondrial = config['mitochondrial_refs'],
        #housekeeping = config['housekeeping_refs']
    output:
        index_ct = config['bs_index']['ct'],
        index_ga = config['bs_index']['ga']
    shell:
        "segemehl.x -x {output.index_ct} -y {output.index_ga} -F 1 -d {input.genomic} {input.mitochondrial}" # {input.housekeeping}"


rule mapping:
    input:
        index = config['index'],
        index_ct = config['bs_index']['ct'],
        index_ga = config['bs_index']['ga'],
        fastqT = 'resources/filteres-reads/a_u_t_q_{sample}.fastq',
        fastqF = 'resources/filteres-reads/a_u_f_q_{sample}.fastq',
        fastq = 'resources/filteres-reads/a_u_q_{sample}.fastq',
        #fasta_files = [config['genomic_refs'], config['mitochondrial_refs'], config['housekeeping_refs']],
        fasta_files = [config['genomic_refs'], config['mitochondrial_refs']],
    output:
        sam = 'resources/mapping/all/polyTtrimming_{polyT}/{sample}.sam',
        unmapped_fastq = 'resources/mapping/all/polyTtrimming_{polyT}/{sample}_unmapped.fastq',
        mapping_stats = "qc/mapping/all_polyTtrimming_{polyT}_{sample}.tsv"
    #resources:
        #mem_mb = 500
    threads: 2
    run:
        import subprocess
        import pandas as pd

        def get_mapping_stats(data, filepath, sample_name):
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

        if wildcards.polyT == "OFF":
            fastq = input.fastq
        elif wildcards.polyT == "ONtttt":
            fastq = input.fastqT
        elif wildcards.polyT == "ONfixed":
            fastq = input.fastqF
        else:
            print('no valide polyT wildcard')

        call_args = ['segemehl.x',
                     '-i', input.index,
                     '-d' ]  + input.fasta_files + [
                     '-q', fastq,
                     '--accuracy', '85',
                     '-M', '100',
                     '-D', '1',
                     '-t' , str(threads),
                     '-o', output.sam,
                     '-u', output.unmapped_fastq,
                    ]
        if sample_dict[wildcards.sample]['treatment'] == 'BS':
            print(wildcards.sample, 'BS')
            call_args = ['segemehl.x',
                         '-i', input.index_ct,
                         '-j',input.index_ga,
                         '-d' ]  + input.fasta_files + [
                         '-q', fastq,
                         '--accuracy', '90',
                         '-M', '100',
                         '-D', '1',
                         #'-E', '7',
                         '-t' , str(threads),
                         '-o', output.sam,
                         '-u', output.unmapped_fastq,
                         '-F', '1',
                        ]
        if sample_dict[wildcards.sample]['treatment'] == 'MOCK':
            call_args = ['segemehl.x',
                                 '-i', input.index,
                                 '-d' ]  + input.fasta_files + [
                                 '-q', fastq,
                                  '--accuracy', '80',
                                 '-M', '100',
                                 '-D', '1',
                                 '-t' , str(threads),
                                 '-o', output.sam,
                                 '-u', output.unmapped_fastq,
                                ]


        results = subprocess.run(call_args, capture_output=True)
        data = results.stderr.decode()
        print(data)
        mapping_stat = get_mapping_stats(data, output.mapping_stats, wildcards.sample)
