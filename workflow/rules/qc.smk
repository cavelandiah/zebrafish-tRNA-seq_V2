
# Author: Maria Waldl • code@waldl.org
# Version: 2024-01-24

import json
import pandas as pd
import matplotlib.pyplot as plt


def print_mapping_summary(input_file, mapping_ids,title, summary_tsv, summary_plot):
    # Author: Maria Waldl • code@waldl.org
    # Version: 2024-01-24
    data = []
    tsvs = input_file
    for i, f in enumerate(tsvs):
        df = pd.read_csv(f, sep="\t")

        df.rename(columns={df.columns[0]:'sample'},inplace=True)

        mapping_id = mapping_ids[i]

        df["mapping_id"] = mapping_id
        df["mapped fraction [all reads]"] = df["mapped"] / df["total"]
        df["unique fraction [all reads]"] = df["unique"] / df["total"]
        df["multi fraction [all reads]"] = df["multi"] / df["total"]
        df["unique fraction [all mapped]"] = df["unique"] / df["mapped"]
        df["multi fraction [all mapped]"] = df["multi"] / df["mapped"]

        data.append(df)
    df = pd.concat(data)
    df.sort_values(by="mapping_id", inplace=True)

    fig, axs = plt.subplots(figsize=(len(df) * 0.6 * CM, 12 * CM), nrows=3, ncols=1)
    fig.subplots_adjust(hspace=0.0)

    df.plot(x="mapping_id", y=["total", "mapped", "unique"], ax=axs[0])
    axs[1].axhline(y=0.9, color="grey", linestyle="--")
    axs[1].axhline(y=0.8, color="lightgrey", linestyle="--")

    df.plot(
        x="mapping_id",
        y=["mapped fraction [all reads]", "unique fraction [all reads]"],
        ax=axs[1],
    )

    df.plot(x="mapping_id", y=["multi fraction [all mapped]"], ax=axs[2])

    axs[1].set_ylim(0, 1)
    axs[2].set_ylim(0, 1)

    axs[0].set_xticks(list(range(0, len(df))))
    axs[1].set_xticks(list(range(0, len(df))))
    axs[2].set_xticks(list(range(0, len(df))))

    axs[0].grid(axis="x")
    axs[1].grid(axis="x")
    axs[2].grid(axis="x")

    axs[0].set_xticklabels(["" for l in axs[0].get_xticklabels()], rotation=90)
    axs[1].set_xticklabels(["" for l in axs[1].get_xticklabels()], rotation=90)
    axs[2].set_xticklabels(df["mapping_id"].to_list(), rotation=90)
    plt.suptitle(title)

    fig.savefig(summary_plot, bbox_inches="tight")
    df.set_index('mapping_id', inplace=True)
    df.to_csv(summary_tsv, sep="\t")


rule mapping_summary_compareDM_ref_sets:
    input:
        tsvs = expand( "qc/mapping/{ref_set}_polyTtrimming_{polyT}_{sample}.tsv", ref_set = ['manual'], polyT = config['poly_T_processing'], sample=dm_samples),
    output:
        pdf = 'qc/summary/mapping-compare-ref_sets.pdf',
        tsv = 'qc/summary/mapping-compare-ref_sets.tsv'
    run:
        title = f"comparison of refernce sets on demethy;ated samples (5'T trimming method: {config['poly_T_processing']})"
        mapping_ids =  [tsv.split('/')[-1].replace('.tsv', '') for tsv in input.tsvs]
        mapping_ids = [{'ref_set':mid.split('_')[0], 'Ttrimming':mid.split('_')[2], 'sample': mid.split('_', 3)[-1]} for mid in mapping_ids]
        mapping_ids = [ ' '.join([sample_dict[mid['sample']]['treatment'],str(sample_dict[mid['sample']]['timepoint']), mid['sample'],  mid['ref_set']]) for mid in mapping_ids]
        print_mapping_summary(input.tsvs, mapping_ids,title, output.tsv, output.pdf)

rule mapping_summary_all_samples:
    input:
        tsvs = expand( "qc/mapping/{ref_set}_polyTtrimming_{polyT}_{sample}.tsv", ref_set='{ref_set}', polyT=config['poly_T_processing'], sample=SAMPLES),
    output:
        pdf = 'qc/mapping-summary/mapping-{ref_set}.pdf',
        tsv = 'qc/mapping-summary/mapping-{ref_set}.tsv',
    run:
        title = f"summary {wildcards.ref_set} (5'T trimming method: {config['poly_T_processing']})"
        mapping_ids =  [tsv.split('/')[-1].replace('.tsv', '') for tsv in input.tsvs]
        mapping_ids = [{'ref_set':mid.split('_')[0], 'Ttrimming':mid.split('_')[2], 'sample': mid.split('_', 3)[-1]} for mid in mapping_ids]
        mapping_ids = [ ' '.join([sample_dict[mid['sample']]['treatment'],str(sample_dict[mid['sample']]['timepoint']), mid['sample']]) for mid in mapping_ids]
        print_mapping_summary(input.tsvs, mapping_ids,title, output.tsv, output.pdf)


rule mapping_summary_compare_Ttrimming:
    input:
        tsvs = expand( "qc/mapping/{ref_set}_polyTtrimming_{polyT}_{sample}.tsv", ref_set = '{ref_set}' ,sample=SAMPLES, polyT=['ONtttt']),
    output:
        pdf = 'qc/summary/mapping-Ttrimming-{ref_set}.pdf',
        tsv = 'qc/summary/mapping-Ttrimming-{ref_set}.tsv'
    run:
        title = f"comparison of 5'T trimming method: on references set '{wildcards.ref_set}'"
        mapping_ids =  [tsv.split('/')[-1].replace('.tsv', '') for tsv in input.tsvs]
        mapping_ids = [{'ref_set':mid.split('_')[0], 'Ttrimming':mid.split('_')[2], 'sample': mid.split('_')[3]} for mid in mapping_ids]
        mapping_ids = [ ' '.join([sample_dict[mid['sample']]['treatment'],str(sample_dict[mid['sample']]['timepoint']), mid['sample'], mid['Ttrimming']]) for mid in mapping_ids]
        print_mapping_summary(input.tsvs, mapping_ids,title, output.tsv, output.pdf)


rule preprocessing_summary: #TODO fastqc results summary
    input:
        raw_reports = expand("qc/fastqc/raw_{sample}.json", sample=SAMPLES),
        adapter_trimming_reports = expand("qc/trimming/a_{sample}_cutadapt.json", sample=SAMPLES),
        CCA_trimming_report = expand("qc/trimming/a_u_c_{sample}_cutadapt.json", sample=SAMPLES),
        polyT_trimming_report = expand("qc/trimming/a_u_t_{sample}_cutadapt.json", sample=SAMPLES),
        qc_without_T_trimming_reports  = expand("qc/fastqc/a_u_q_{sample}.json", sample=SAMPLES),
        qc_after_T_trimming_reports = expand("qc/fastqc/a_u_t_q_{sample}.json", sample=SAMPLES),
    output:
        pdf_summary = 'qc/preprocessing_summary/preproccessing.pdf',
        tsv_summaty = 'qc/preprocessing_summary/preprocessing.tsv'
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

        def get_fastp_reads(fastp_reports_path, raw_fastq, data_field, prefix):
            ''' Possibel data_field :
            'passed_filter_reads' 'low_quality_reads' "too_many_N_reads" "low_complexity_reads" "too_short_reads"
            '''
            json_file = os.path.join(fastp_reports_path, prefix+raw_fastq+'.json')
            with open(json_file) as f:
                data=json.load(f)

            if data_field not in data['filtering_result'].keys():
                print(prefix, data_field, raw_fastq)
                return 0
            return data['filtering_result'][data_field]

        def get_cca_reads(trimming_reports, raw_fastq):
            json_file = os.path.join(trimming_reports, 'a_u_c_'+raw_fastq+'_cutadapt.json')
            with open(json_file) as f:
                data=json.load(f)
            return data['read_counts']['read1_with_adapter']

        fastp_reports = '/'.join(input.qc_after_T_trimming_reports[0].split('/')[0:-1])
        trimming_reports = '/'.join(input.polyT_trimming_report[0].split('/')[0:-1])

        df = pd.read_csv(config["samples_tsv"], sep="\t",)# index_col = 'fastq')



        df['raw reads'] = df.apply(lambda row: get_total_reads(trimming_reports, row['fastq']), axis = 1)
        df['adapter trimmed reads'] = df.apply(lambda row: get_adapter_reads(trimming_reports, row['fastq']), axis = 1)
        df['to short reads (after adapter trimming)'] = df.apply(lambda row: get_short_reads_after_adapter_trimming(trimming_reports, row['fastq']), axis = 1)
        df['passed size selection (>=29nts)'] = df.apply(lambda row: get_notshort_reads_after_adapter_trimming(trimming_reports, row['fastq']), axis = 1)
        df['remaining after polyT trimming + size selection (>=20nts)'] = df.apply(lambda row: get_notshort_reads_after_polyT(trimming_reports, row['fastq']), axis = 1)
        df['fastp qc passing reads'] = df.apply(lambda row: get_fastp_reads(fastp_reports, row['fastq'], 'passed_filter_reads', 'a_u_t_q_'), axis = 1)
        df['fastp qc low quality (<Q15) reads'] = df.apply(lambda row: get_fastp_reads(fastp_reports, row['fastq'], 'low_quality_reads', 'a_u_t_q_'), axis = 1)
        df['fastp qc too many N reads'] = df.apply(lambda row: get_fastp_reads(fastp_reports, row['fastq'], 'too_many_N_reads', 'a_u_t_q_'), axis = 1)
        df['fastp qc too low complexity reads'] = df.apply(lambda row: get_fastp_reads(fastp_reports, row['fastq'], 'low_complexity_reads', 'a_u_t_q_'), axis = 1)
        df['fastp qc too short reads (<20nts)'] = df.apply(lambda row: get_fastp_reads(fastp_reports, row['fastq'], 'too_short_reads', 'a_u_t_q_'), axis = 1)
        df['raw fastp qc passing reads'] = df.apply(lambda row: get_fastp_reads(fastp_reports, row['fastq'], 'passed_filter_reads', 'raw_'), axis = 1)
        df['raw fastp qc low quality (<Q15) reads'] = df.apply(lambda row: get_fastp_reads(fastp_reports, row['fastq'], 'low_quality_reads', 'raw_'), axis = 1)
        #df['raw fastp qc too low complexity reads'] = df.apply(lambda row: get_fastp_reads(fastp_reports, row['fastq'], 'low_complexity_reads', 'raw_'), axis = 1)
        df['raw fastp qc too many N reads'] = df.apply(lambda row: get_fastp_reads(fastp_reports, row['fastq'], 'too_many_N_reads', 'raw_'), axis = 1)
        df['CCA/TTA containing reads (after adapter and umi trimming)'] = df.apply(lambda row: get_cca_reads(trimming_reports, row['fastq']), axis = 1)


        fig, axs = plt.subplots(figsize=((len(df)*0.7+1.5)*CM, 8*CM), nrows= 1, ncols=1)
        fig.tight_layout()
        #fig.subplots_adjust(hspace=0.6)
        #df.reset_index(inplace = True)
        df.plot(x='fastq', y=['raw reads',
                                   'adapter trimmed reads',
                                   'passed size selection (>=29nts)',
                                   'remaining after polyT trimming + size selection (>=20nts)',
                                   'fastp qc passing reads',
                                   #'CCA/TTA containing reads (after adapter and umi trimming)',
                                  ],
                     ax = axs)

        axs.set_xticks(ticks = list(range(0,len(df))))
        xlabels = df['fastq'].to_list()
        time_names = df['timepoint_name'].to_list()
        treatment = df['treatment'].to_list()
        xlabels = [x+' '+time_names[e]+', '+treatment[e] for e,x in enumerate(xlabels)]
        axs.set_xticklabels(xlabels, rotation=90)
        axs.set_ylabel('read count')

        df_with_stats = pd.concat([df, df.describe()])
        df_with_stats.to_csv(output.tsv_summaty, sep="\t", index=False)


        fig.savefig(os.path.join(output.pdf_summary),bbox_inches='tight')
