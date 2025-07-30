
# Author: Maria Waldl • code@waldl.org
# Version: 2024-01-24
# Corrected: 2025-07-30


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
                     '--low_complexity_filter',
                     '--complexity_threshold=30',
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

rule trim_barcode_3_end:
    input:
        fastq = "resources/filteres-reads/a_{sample}.fastq",
    output:
        fastq = "resources/filteres-reads/a_b_{sample}.fastq",
        json = "qc/trimming/a_b_{sample}_cutadapt.json"
    run:
        import os
        import subprocess
        
        barcode = sample_dict[wildcards.sample]['barcode'] 
        barcode3 = barcode.upper()
        call_args = ['cutadapt',
                     #'-j '+ str(cores), # number of cores
                     '-q 10', #3' end quality trimm; before adapter trimming
                     '-e 0.15', # maximum error rate in adapter
                     '--minimum-length=17', # use -m 1 if downsteram tools have problems with zero length reads
                     '-a', barcode3, # 3' adapter
                     '--discard-untrimmed',
                     '-o', output.fastq,
                     #'--info-file=' +  os.path.join(fastp_reports, 'cutadapy_info_'+raw_fastq+'.tsv'),
                     '--json=' +output.json,
                     input.fastq,
                      ]
        subprocess.run(call_args)


rule get_umi:
    input:
        "resources/filteres-reads/a_b_{sample}.fastq",
        #"resources/filteres-reads/a_{sample}.fastq"
    output:
        fastq = "resources/filteres-reads/a_u_{sample,[A-Za-z0-9_]+}.fastq",
        log = "qc/trimming/a_u_{sample}_umi-tools.log"
    shell:
        # Author: Maria Waldl • code@waldl.org
        # Corrected: 2025-07-30
        'umi_tools extract --stdin={input} --extract-method=regex --bc-pattern="^(?P<umi_1>.{3}).+(?P<umi_2>.{6})$" -L {output.log} --stdout {output.fastq}'


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
                 '--minimum-length=20', # use -m 1 if downsteram tools have problems with zero length reads
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
