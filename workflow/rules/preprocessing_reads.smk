
rule raw_fastqc: # only used in qc overview, not for preprocessing before mapping
    input:
        config['raw_reads_folder']+"/{sample}.R1.fastq"
    output:
        html="qc/fastqc/raw_{sample}.html",
        json="qc/fastqc/raw_{sample}.json"
    run:
        import os
        import subprocess
        import shutil
        call_args = ['fastp',
                     '-i', input[0],
                     '--trim_tail1=1',
                     #'--trim_front1=5',
                     '--length_required=20',
                     '--trim_poly_x',
                     '-p',
                     '-h', output.html]
        subprocess.run(call_args)
        shutil.copyfile('fastp.json', output.json)

        html="qc/fastqc/a_u_t_q_{sample}.html",
        json="qc/fastqc/a_u_t_q_{sample}.json",
        fastq = "resources/filteres-reads/a_u_t_q_{sample}.fastq",


rule trim_adapter_only:
    input:
        config['raw_reads_folder']+"/{sample}.R1.fastq"
    output:
        fastq = "resources/filteres-reads/a_{sample}.fastq",
        json = "qc/trimming/a_{sample}_cutadapt.json"
    run:
        import os
        import subprocess
        adapter3 = 'CTGTAGGCACCATCAAT;min_overlap=6'
        if sample_dict[wildcards.sample]['length'] < 100:
            adapter3 = 'CTGTAGGCACCATCAAT'
        call_args = ['cutadapt',
                     #'-j '+ str(cores), # number of cores
                     '-q 10', #3' end quality trimm; before adapter trimming
                     '-e 0.15', # maximum error rate in adapter
                     '--minimum-length=29', # use -m 1 if downsteram tools have problems with zero length reads
                     '-a', adapter3, # 3' adapter
                     '--discard-untrimmed',
                     '-o', output.fastq,
                     #'--info-file=' +  os.path.join(fastp_reports, 'cutadapy_info_'+raw_fastq+'.tsv'),
                     '--json=' +output.json,
                     input[0],
                      ]
        subprocess.run(call_args)

rule get_umi:
    input:
        "resources/filteres-reads/a_{sample}.fastq"
    output:
        fastq = "resources/filteres-reads/a_u_{sample,[A-Za-z0-9]+}.fastq",
        log = "qc/trimming/a_u_{sample}_umi-tools.log"
    shell:
        'umi_tools extract --stdin={input} --extract-method=regex --bc-pattern="(?P<umi_1>.{{3}})([AGCTUN].*)(?P<umi_2>.{{6}})$" -L {output.log} --stdout {output.fastq}'


rule trim_cca: # currently only used for statistics, not in preprocessing for mapping
    input:
        "resources/filteres-reads/a_u_{sample}.fastq"
    output:
        fastq = "resources/filteres-reads/a_u_c_{sample,[A-Za-z0-9]+}.fastq",
        json = "qc/trimming/a_u_c_{sample}_cutadapt.json"
    run:
        import os
        import subprocess
        to_trim = 'CCAX'
        if sample_dict[wildcards.sample]['treatment'] == 'BS':
            to_trim = 'TTAX'
        call_args = ['cutadapt',
                 #'-j '+ str(cores), # number of cores
                 '-a', to_trim+';min_overlap=3',
                 '--minimum-length=17', # use -m 1 if downsteram tools have problems with zero length reads
                 '-o', output.fastq,
                     #'--info-file=' +  os.path.join(fastp_reports, 'cutadapy_info_'+raw_fastq+'.tsv'),
                 '--json=' +output.json,
                 input[0],
                 ]
        subprocess.run(call_args)


rule quality_filter:
    input:
        "resources/filteres-reads/a_u_{sample}.fastq",
    output:
        html="qc/fastqc/a_u_q_{sample}.html",
        json="qc/fastqc/a_u_q_{sample}.json",
        fastq = "resources/filteres-reads/a_u_q_{sample}.fastq",
    run:
        import os
        import subprocess
        import shutil
        call_args = ['fastp',
                     '-i', input[0],
                     '--length_required=20',
                     '--trim_poly_x',
                     'poly_x_min_len=7',
                     '-p',
                     '-A',
                     '--low_complexity_filter',
                     '--complexity_threshold=30',
                     '-o', output.fastq,
                     '-h', output.html]
        subprocess.run(call_args)
        shutil.copyfile('fastp.json', output.json)

rule dupl_fastqc: #TODO currently not used
    input:
        "resources/filteres-reads/a_{sample}.fastq"
    output:
        html="qc/fastqc/dupl_{sample}.html",
        json="qc/fastqc/dupl_{sample}_fastq.json"
    run:
        import os
        import subprocess
        import shutil
        call_args = ['fastp',
                     '-i', input[0],
                     '--length_required=20',
                     '-p',
                     #'--quiet',
                     #https://pypi.org/project/umitools/
                     #https://github.com/weng-lab/umitools
                     #https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4933-1
                     #https://www.nature.com/articles/s41598-020-71323-0
                     #https://umi-tools.readthedocs.io/en/latest/reference/count_tab.html
                     #https://pubmed.ncbi.nlm.nih.gov/29368091/
                     '-U', '--umi_loc=reahttps://pypi.org/project/umitools/d1', '--umi_len=3',
                     '-h', output.html]
        subprocess.run(call_args)
        shutil.copyfile('fastp.json', output.json)


rule trim_5prime_polyT:
    input:
        "resources/filteres-reads/a_u_{sample}.fastq",
        "resources/filteres-reads/a_u_{sample}.fastq",
    output:
        fastq = "resources/filteres-reads/a_u_t_{sample}.fastq",
        json = "qc/trimming/a_u_t_{sample}_cutadapt.json"
    run:
        import os
        import subprocess

        to_trim = 'XTTTTTTTTTTTTTTT'
        if sample_dict[wildcards.sample]['treatment'] == "BS":
            to_trim = 'XTTT'

        call_args = ['cutadapt',
                 #'-j '+ str(cores), # number of cores
                 '-g', to_trim+';min_overlap=1',
                 '-e 0.0001', # maximum error rate in adapter
                 '--minimum-length=20', # use -m 1 if downsteram tools have problems with zero length reads
                 '-o', output.fastq,
                     #'--info-file=' +  os.path.join(fastp_reports, 'cutadapy_info_'+raw_fastq+'.tsv'),
                 '--json=' +output.json,
                 input[0],
                 ]
        subprocess.run(call_args)

rule trim_5prime_fixed:
    input:
        "resources/filteres-reads/a_u_{sample}.fastq",
        "resources/filteres-reads/a_u_{sample}.fastq",
    output:
        fastq = "resources/filteres-reads/a_u_f_{sample}.fastq",
        json = "qc/trimming/a_u_f_{sample}_cutadapt.json"
    run:
        import os
        import subprocess

        call_args = ['cutadapt',
                 #'-j '+ str(cores), # number of cores
                 '--cut=4',
                 '--minimum-length=17', # use -m 1 if downsteram tools have problems with zero length reads
                 '-o', output.fastq,
                 '--json=' +output.json,
                 input[0],
                 ]
        subprocess.run(call_args)


rule quality_filter_after_t_trimming:
    input:
        "resources/filteres-reads/a_u_t_{sample}.fastq",
    output:
        html="qc/fastqc/a_u_t_q_{sample}.html",
        json="qc/fastqc/a_u_t_q_{sample}.json",
        fastq = "resources/filteres-reads/a_u_t_q_{sample}.fastq",
    run:
        import os
        import subprocess
        import shutil
        call_args = ['fastp',
                     '-i', input[0],
                     '--length_required=20',
                     '--trim_poly_x',
                     'poly_x_min_len=7',
                     '-p',
                     '-A',
                     '--low_complexity_filter',
                     '--complexity_threshold=30',
                     '-o', output.fastq,
                     '-h', output.html]
        subprocess.run(call_args)
        shutil.copyfile('fastp.json', output.json)

rule quality_filter_after_fixed_5prime_trimming:
    input:
        "resources/filteres-reads/a_u_f_{sample}.fastq",
    output:
        html="qc/fastqc/a_u_f_q_{sample}.html",
        json="qc/fastqc/a_u_f_q_{sample}.json",
        fastq = "resources/filteres-reads/a_u_f_q_{sample}.fastq",
    run:
        import os
        import subprocess
        import shutil
        call_args = ['fastp',
                     '-i', input[0],
                     '--length_required=17',
                     '--trim_poly_x',
                     'poly_x_min_len=7',
                     '-p',
                     '-A',
                     '--low_complexity_filter',
                     '--complexity_threshold=30',
                     '-o', output.fastq,
                     '-h', output.html]
        subprocess.run(call_args)
        shutil.copyfile('fastp.json', output.json)

rule preprocessing_summary: #TODO fastqc results summary
    input:
        raw_reports = expand("qc/fastqc/raw_{sample}.json", sample=SAMPLES),
        adapter_trimming_reports = expand("qc/trimming/a_{sample}_cutadapt.json", sample=SAMPLES),
        CCA_trimming_report = expand("qc/trimming/a_u_c_{sample}_cutadapt.json", sample=SAMPLES),
        polyT_trimming_report = expand("qc/trimming/a_u_t_{sample}_cutadapt.json", sample=SAMPLES),
        qc_without_T_trimming_reports  = expand("qc/fastqc/a_u_q_{sample}.json", sample=SAMPLES),
        qc_after_T_trimming_reports = expand("qc/fastqc/a_u_t_q_{sample}.json", sample=SAMPLES),
    output:
        pdf_summary = 'results/qc/preproccessing_summary.pdf'
    run:
        import os
        import matplotlib.pyplot as plt

        def get_total_reads(trimming_reports_dir, raw_fastq):
            json_file = os.path.join(trimming_reports_dir, 'a_'+raw_fastq+'_cutadapt.json')
            with open(json_file) as f:
                data=json.load(f)
            return data['read_counts']['input']

        def get_adapter_reads(trimming_reports_dir, raw_fastq):
            json_file = os.path.join(trimming_reports_dir, 'a_'+raw_fastq+'_cutadapt.json')
            with open(json_file) as f:
                data=json.load(f)
            return data['read_counts']['read1_with_adapter']

        def get_short_reads_after_adapter_trimming(trimming_reports, raw_fastq):
            json_file = os.path.join(trimming_reports, 'a_'+raw_fastq+'_cutadapt.json')
            with open(json_file) as f:
                data=json.load(f)
            return data['read_counts']['filtered']['too_short']

        def get_notshort_reads_after_adapter_trimming(trimming_reports, raw_fastq):
            json_file = os.path.join(trimming_reports, 'a_'+raw_fastq+'_cutadapt.json')
            with open(json_file) as f:
                data=json.load(f)
            return data['read_counts']['output']

        def get_notshort_reads_after_polyT(trimming_reports, raw_fastq):
            json_file = os.path.join(trimming_reports, 'a_u_t_'+raw_fastq+'_cutadapt.json')
            with open(json_file) as f:
                data=json.load(f)
            return data['read_counts']['output']

        def get_fastp_reads(fastp_reports, raw_fastq):
            json_file = os.path.join(fastp_reports, 'a_u_t_q_'+raw_fastq+'.json')
            with open(json_file) as f:
                data=json.load(f)
            return data['filtering_result']['passed_filter_reads']

        def get_cca_reads(trimming_reports, raw_fastq):
            json_file = os.path.join(trimming_reports, 'a_u_c_'+raw_fastq+'_cutadapt.json')
            with open(json_file) as f:
                data=json.load(f)
            return data['read_counts']['read1_with_adapter']


        fastp_reports = '/'.join(input.qc_after_T_trimming_reports[0].split('/')[0:-1])
        print(fastp_reports)
        trimming_reports = '/'.join(input.polyT_trimming_report[0].split('/')[0:-1])
        print(trimming_reports)

        df = pd.read_csv(config["samples_tsv"], sep="\t",)# index_col = 'fastq')

        fig, axs = plt.subplots(figsize=((len(df)*1+1.5)*CM, 9*CM), nrows= 1, ncols=1)
        fig.tight_layout()
        #fig.subplots_adjust(hspace=0.6)
        #df.reset_index(inplace = True)

        df['raw reads'] = df.apply(lambda row: get_total_reads(trimming_reports, row['fastq']), axis = 1)
        df['adapter trimmed reads'] = df.apply(lambda row: get_adapter_reads(trimming_reports, row['fastq']), axis = 1)
        df['to short reads (after adapter trimming)'] = df.apply(lambda row: get_short_reads_after_adapter_trimming(trimming_reports, row['fastq']), axis = 1)
        df['remaining reads after adapter trimming (>=29nts)'] = df.apply(lambda row: get_notshort_reads_after_adapter_trimming(trimming_reports, row['fastq']), axis = 1)
        df['remaining after polyT trimming'] = df.apply(lambda row: get_notshort_reads_after_polyT(trimming_reports, row['fastq']), axis = 1)
        df['fastp qc passed reads'] = df.apply(lambda row: get_fastp_reads(fastp_reports, row['fastq']), axis = 1)
        df['CCA/TTA containing reads (after adapter and umi trimming)'] = df.apply(lambda row: get_cca_reads(trimming_reports, row['fastq']), axis = 1)

        print(df.head())

        df.plot(x='fastq', y=['raw reads',
                                   'adapter trimmed reads',
                                   'remaining reads after adapter trimming (>=29nts)',
                                   'remaining after polyT trimming',
                                   'fastp qc passed reads',
                                   'CCA/TTA containing reads (after adapter and umi trimming)',
                                  ],
                     ax = axs)

        axs.set_xticks(ticks = list(range(0,len(df))))
        xlabels = df['fastq'].to_list()
        time_names = df['timepoint_name'].to_list()
        treatment = df['treatment'].to_list()
        xlabels = [x+'\n'+time_names[e]+', '+treatment[e] for e,x in enumerate(xlabels)]
        axs.set_xticklabels(xlabels, rotation=90)
        axs.set_ylabel('read count')

        fig.savefig(os.path.join(output.pdf_summary),bbox_inches='tight')
