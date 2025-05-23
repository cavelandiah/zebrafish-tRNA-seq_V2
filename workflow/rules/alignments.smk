# Author: Maria Waldl • code@waldl.org
# Version: 2024-01-24

rule align_to_cm:
    input:
        rfam_alignment = config['rfam_alignment'],
        cm = config['rfam_cm'],
        mitochondrial_fasta = config['mitochondrial_refs'],
        genomic_split_hc_fasta = config['genomic_split_hc_refs'],
        genomic_split_additional_fasta = config['genomic_split_additional_refs'],
    output:
        #canonical_summary = 'resources/references/canonical_tRNA/position_map.txt',
        mitochondrial_alignment = config['mitochondrial_alignment'],
        genomic_hc_alignment = config['genomic_hc_alignment'],
        genomic_additional_alignment = config['genomic_additional_alignment'],
    threads: 4
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        import os
        import subprocess
        import shutil

        call_args = ['cmalign',
                     '--mapali', input.rfam_alignment,
                     '--outformat', 'AFA',
                     #'--outformat', 'Stockholm',
                     #'--outformat', 'Pfam',
                     '-o', output.mitochondrial_alignment,
                     '--cpu', str(threads),
                     '-g',
                    input.cm,
                    input.mitochondrial_fasta,
                    ]

        subprocess.run(call_args)

        call_args = ['cmalign',
                     '--mapali', input.rfam_alignment,
                     '--outformat', 'AFA',
                     #'--outformat', 'Stockholm',
                     #'--outformat', 'Pfam',
                     '-o', output.genomic_hc_alignment,
                     '--cpu', str(threads),
                     '-g',
                    input.cm,
                    input.genomic_split_hc_fasta ,
                    ]

        subprocess.run(call_args)

        call_args = ['cmalign',
                     '--mapali', input.rfam_alignment,
                     '--outformat', 'AFA',
                     #'--outformat', 'Stockholm',
                     #'--outformat', 'Pfam',
                     '-o', output.genomic_additional_alignment,
                     '--cpu', str(threads),
                     '-g',
                    input.cm,
                    input.genomic_split_additional_fasta,
                    ]

        subprocess.run(call_args)


rule get_min_cov_refs_for_alignment:
    input:
        cov_sum = 'resources/coverage/pre-filter_umi/raw/min_coverage_summary_DM.tsv'
    output:
        keep_refs = 'resources/min_coverage_refs/pre-filter_umi/min_cov_refs_alignment.yaml'
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        import pandas as pd
        import yaml
        df = pd.read_csv(input.cov_sum, sep="\t")
        refs = []

        for criterium, cutoff in config['min_raw_abundance_count_for_alignment']:
            selected = df[df[criterium]>=cutoff]['RNAME'].to_list()
            print(criterium, cutoff, len(selected))
            refs += selected
        refs = list(set(refs))
        with open(output.keep_refs, 'w') as file:
            outputs = yaml.dump(refs, file)
        print(len(df))
        print( len(refs))


rule get_min_raw_coverage_fasta: # cm align can not align more than 10 000 seqeunces at a time so prefiltering is needed,TODO: change to reuse functions in coverage.smk
    input:
        min_cov_ids = 'resources/min_coverage_refs/pre-filter_umi/min_cov_refs_alignment.yaml',
        tRNAs_fa = 'resources/references/all_refs.fa'
    output:
        fasta = config['min_raw_abundance_refs']
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        import yaml
        from Bio import SeqIO
        from Bio.SeqRecord import SeqRecord

        with open(input.min_cov_ids) as file:
            refs = yaml.safe_load(file)

        filtered_sequences = []
        for record in SeqIO.parse(input.tRNAs_fa, "fasta"):
            if record.id in refs:
                seq = record.seq[0:-3]
                record = SeqRecord(seq, id = record.id, name = '', description= '')
                filtered_sequences.append(record)

        SeqIO.write(filtered_sequences, output.fasta, "fasta")





rule align_min_coverage_refs_to_cm:
    input:
        rfam_alignment = config['rfam_alignment'],
        cm = config['rfam_cm'],
        min_coverage_fasta = config['min_raw_abundance_refs'],
    output:
        min_coverage_alignment = config['min_raw_abundance_alignment'].replace('.fa', '')+'without_CCA.fa', # Rfam alignment contains no CCA
    threads: 4
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        import os
        import subprocess
        import shutil

        call_args = ['cmalign',
                     '--mapali', input.rfam_alignment,
                     '--outformat', 'AFA',
                     #'--outformat', 'Stockholm',
                     #'--outformat', 'Pfam',
                     '-o', output.min_coverage_alignment ,
                     '--cpu', str(threads),
                     '-g',
                    input.cm,
                    input.min_coverage_fasta,
                    ]

        subprocess.run(call_args)

rule add_CCA_to_alignment:
    input:
        #resources/references/alignment/min_raw_abundance_tRNAswithout_CCA.fa
        minimal_coverage_alignment = config['min_raw_abundance_alignment'].replace('.fa', '')+'without_CCA.fa',
        selected_ids = 'resources/min_coverage_refs/pre-filter_'+config['reads_filter']+'/min_cov_refs.yaml',
    output:
        #resources/references/alignment/min_raw_abundance_tRNAs.fa
        minimal_coverage_alignment = config['min_raw_abundance_alignment'],
        selected_alignment = config['selected_alignment'],
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        from Bio import SeqIO
        from Bio.SeqRecord import SeqRecord
        import yaml

        with open(input.selected_ids) as file:
            selected_refs = yaml.safe_load(file)

        extended_sequences = []
        with open(output.selected_alignment, 'w') as selected_out:
            for record in SeqIO.parse(input.minimal_coverage_alignment, "fasta"):
                seq = record.seq+'CCA'
                record = SeqRecord(seq, id = record.id, name = '', description= '')
                extended_sequences.append(record)
                # get additional file with only selected refs and without line breaks in sequence
                if record.id in selected_refs:
                    selected_out.write('>' + record.id + '\n')
                    selected_out.write(str(seq) + '\n')
                if '/' in record.id:
                    selected_out.write('>' + record.id + '\n')
                    selected_out.write(str(seq) + '\n')

        SeqIO.write(extended_sequences, output.minimal_coverage_alignment, "fasta")


rule get_sorted_selected_alignment:
    input:
        selected_alignment = config['selected_alignment'],
        selected_ids = 'resources/min_coverage_refs/pre-filter_'+config['reads_filter']+'/min_cov_refs.yaml',
    output:
        selected_alignment = config['selected_alignment']+'_sorted.fa',
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        from Bio import SeqIO
        from Bio.SeqRecord import SeqRecord
        import yaml

        with open(input.selected_ids) as file:
            selected_refs = yaml.safe_load(file)

        with open(output.selected_alignment, 'w') as selected_out:
            records = list(SeqIO.parse(input.selected_alignment, "fasta"))
            records = sorted(records, key = lambda record: record.id)
            for record in records:
                if record.id in selected_refs:
                    selected_out.write('>' + record.id + '\n')
                    selected_out.write(str(record.seq) + '\n')

rule align_manual_refs_to_cm:
    input:
        rfam_alignment = config['rfam_alignment'],
        cm = config['rfam_cm'],
        fasta = 'resources/references/non-redundant/manual-noCCA.fa',
    output:
        manual_alignment = 'resources/references/alignment/manual_tRNAs_noCCA.fa'
    threads: 4
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        import os
        import subprocess
        import shutil

        call_args = ['cmalign',
                     '--mapali', input.rfam_alignment,
                     '--outformat', 'AFA',
                     #'--outformat', 'Stockholm',
                     #'--outformat', 'Pfam',
                     '-o', output.manual_alignment ,
                     '--cpu', str(threads),
                     '-g',
                    input.cm,
                    input.fasta,
                    ]

        subprocess.run(call_args)

rule add_CCA_to_manual_alignment:
    input:
        manual_alignment = 'resources/references/alignment/manual_tRNAs_noCCA.fa'
    output:
        manual_alignment = 'resources/references/alignment/manual_tRNAs.fa'
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        from Bio import SeqIO
        from Bio.SeqRecord import SeqRecord

        with open(output.manual_alignment, 'w') as output_fasta:
            for record in SeqIO.parse(input.manual_alignment, "fasta"):
                seq = record.seq+'CCA'
                output_fasta.write('>' + record.id + '\n')
                output_fasta.write(str(seq) + '\n')
