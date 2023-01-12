
rule build_index_selected:
    input:
        selected_refs = config['seleceted_refs']
    output:
        index = config['selected_index']
    shell:
        "segemehl.x -x {output.index} -d {input.selected_refs}"

rule build_BS_index_selected:
    input:
        selected_refs = config['seleceted_refs']
    output:
        index_ct = config['selected_bs_index']['ct'],
        index_ga = config['selected_bs_index']['ga']
    shell:
        "segemehl.x -x {output.index_ct} -y {output.index_ga} -F 1 -d {input.selected_refs}"



rule mapping_selected:
    input:
        index = config['selected_index'],
        index_ct = config['selected_bs_index']['ct'],
        index_ga = config['selected_bs_index']['ga'],
        fastqT = 'resources/filteres-reads/a_u_t_q_{sample}.fastq',
        fastq = 'resources/filteres-reads/a_u_q_{sample}.fastq',
        fasta_files = [config['seleceted_refs']],
    output:
        sam = 'resources/mapping/selected/polyTtrimming_{polyT}/{sample}.sam',
        unmapped_fastq = 'resources/mapping/selected/polyTtrimming_{polyT}/{sample}_unmapped.fastq',
        mapping_stats = "qc/mapping/selected_polyTtrimming_{polyT}_{sample}.tsv"
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
                     '-M', '200',
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
                         '-M', '200',
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
                                 '-M', '200',
                                 '-D', '1',
                                 '-t' , str(threads),
                                 '-o', output.sam,
                                 '-u', output.unmapped_fastq,
                                ]


        results = subprocess.run(call_args, capture_output=True)
        data = results.stderr.decode()
        print(data)
        mapping_stat = get_mapping_stats(data, output.mapping_stats, wildcards.sample)


rule mapping_summary:
    input:
        all_ref_tsvs = expand( "qc/mapping/all_polyTtrimming_{polyT}_{sample}.tsv", sample=dm_samples, polyT=['ONfixed']),
        all_ref_tsvs_fixed = expand( "qc/mapping/all_polyTtrimming_{polyT}_{sample}.tsv", sample=dm_samples, polyT=['ONtttt']),
        selected_refs_tsvs = expand( "qc/mapping/selected_polyTtrimming_{polyT}_{sample}.tsv", sample=dm_samples, polyT=['ONtttt']),
    output:
        'results/qc/mapping-summary.pdf'
    run:
        import pandas as pd
        import matplotlib.pyplot as plt

        data = []
        tsvs = input.all_ref_tsvs +input.all_ref_tsvs_fixed+ input.selected_refs_tsvs
        for f in tsvs:
            df = pd.read_csv(f, sep = '\t')
            df.drop(columns=df.columns[0], axis=1, inplace=True)
            filename = f.split('/')[-1].replace('.tsv', '')
            sample = filename.split('_')[-1]
            mapping_id = ' '.join([sample_dict[sample]['treatment'], str(sample_dict[sample]['timepoint']), sample, filename.split('_')[0], filename.split('_')[2]])
            df['mapping_id'] = mapping_id
            df['mapped fraction [all reads]'] = df['mapped']/df['total']
            df['unique fraction [all reads]'] = df['unique']/df['total']
            df['multi fraction [all reads]'] = df['multi']/df['total']

            df['unique fraction [all mapped]'] = df['unique']/df['mapped']
            df['multi fraction [all mapped]'] = df['multi']/df['mapped']

            df['mapping_id'] = mapping_id
            data.append(df)
        df =  pd.concat(data)

        df.sort_values(by = 'mapping_id', inplace = True )

        fig, axs = plt.subplots(figsize=(len(df)*0.6*CM, 12*CM), nrows= 3, ncols=1)
        fig.subplots_adjust(hspace=0.0)
        df.plot(x = 'mapping_id', y= ['total','mapped','unique']  ,ax=axs[0])
        axs[1].axhline(y=0.9, color='grey', linestyle='--')
        axs[1].axhline(y=0.8, color='lightgrey', linestyle='--')

        df.plot(x = 'mapping_id', y = ['mapped fraction [all reads]', 'unique fraction [all reads]'],ax=axs[1])

        df.plot(x = 'mapping_id', y = ['multi fraction [all mapped]'],ax=axs[2])

        axs[1].set_ylim(0,1)
        axs[2].set_ylim(0,1)

        axs[0].set_xticks(list(range(0,len(df))))
        axs[1].set_xticks(list(range(0,len(df))))
        axs[2].set_xticks(list(range(0,len(df))))

        axs[0].grid(axis='x')
        axs[1].grid(axis='x')
        axs[2].grid(axis='x')

        axs[0].set_xticklabels(['' for l in axs[0].get_xticklabels()], rotation = 90)
        axs[1].set_xticklabels(['' for l in axs[1].get_xticklabels()], rotation = 90)
        axs[2].set_xticklabels(df['mapping_id'].to_list(), rotation = 90)
        fig.savefig(str(output), bbox_inches="tight")
