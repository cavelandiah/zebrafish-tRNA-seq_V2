
# Author: Maria Waldl • code@waldl.org
# Version: 2024-01-24

rule get_raw_mapped_manual_and_remove_GAmapped:
    input:
        sam = 'resources/mapping/manual/polyTtrimming_ONtttt/{sample}.sam'
    output:
        sam = 'resources/filtered-mappings/pre-filter_none/manual/mapped/{sample}.sam'
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24

        import pandas as pd
        import numpy as np

        if sample_dict[wildcards.sample]['treatment'] == 'BS':
            cols = [
                ["QNAME", 0, "category"],
                #["FLAG", 1, np.int16],
                ["RNAME", 2, "category"],
                ["POS", 3, np.int16],
                #["MAPQ", 4, np.int16],
                ["CIGAR", 5, "category"],
                #["MRNM/RNEXT", 6, "category"],
                #["MPOS/PNEXT", 7, np.int16],
                #["ISIZE/TLEN", 8, np.int16],
                ["SEQ", 9, "category"],
                #["QUAL", 10, "category"],
                #["TAGs0", 11, "category"],
                #["TAGs1", 12, "category"],
                #["TAGs2", 13, "category"],
                #["TAGs3", 14, "category"],
                #["TAGs4", 15, "category"],
                #["TAGs5", 16, "category"],
                ["TAGs6", 17, "category"],
                #["TAGs7", 18, "category"],
                #["TAGs8", 19, "category"]
               ]
        else:
            cols = [
                ["QNAME", 0, "category"],
                #["FLAG", 1, np.int16],
                ["RNAME", 2, "category"],
                ["POS", 3, np.int16],
                #["MAPQ", 4, np.int16],
                ["CIGAR", 5, "category"],
                #["MRNM/RNEXT", 6, "category"],
                #["MPOS/PNEXT", 7, np.int16],
                #["ISIZE/TLEN", 8, np.int16],
                ["SEQ", 9, "category"],
                #["QUAL", 10, "category"],
                #["TAGs0", 11, "category"],
                #["TAGs1", 12, "category"],
                #["TAGs2", 13, "category"],
                #["TAGs3", 14, "category"],
                #["TAGs4", 15, "category"],
                #["TAGs5", 16, "category"],
                #["TAGs6", 17, "category"],
                #["TAGs7", 18, "category"],
                #["TAGs8", 19, "category"]
               ]

        col_nr = [ col[1] for col in cols ]
        col_names = [ col[0] for col in cols ]
        data_typ_dict = {}
        for col in cols:
            data_typ_dict[col[0]] = col[2]

        df = pd.read_csv(
        input.sam,
        sep="\t",
        comment="@",
        usecols=col_nr,
        names=col_names,
        dtype = data_typ_dict,
        header=None,
        )
        df = df[df["RNAME"] != "*"]

        if sample_dict[wildcards.sample]['treatment'] == 'BS':
            # drop any mappings to GA index
            df['CT'] = df.apply(lambda row: row["TAGs6"].endswith('CT'), axis =1 )
            df = df[df["CT"]]
            df.drop(columns = ['CT', 'TAGs6'], inplace = True)

        df.to_csv(output.sam, sep="\t", index=False)

rule filter_umi_on_manual:
    input:
        sam = 'resources/filtered-mappings/pre-filter_none/manual/mapped/{sample}.sam'
    output:
        sam = 'resources/filtered-mappings/pre-filter_umi/manual/mapped/{sample}.sam'
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        import pandas as pd
        df = pd.read_csv(
        input.sam,
        sep="\t",
        )
        header = True
        mode = 'w'
        i=0
        for read, g_df in df.groupby('SEQ'):
            g_df['umi'] = g_df.apply(lambda row: row['QNAME'].split('_')[-1], axis =1)
            g_df.drop_duplicates(subset=['umi', 'RNAME'], inplace = True, keep= 'first')

            g_df.to_csv(output.sam,
                         sep = '\t',
                         index=False,
                         header=header,
                         mode=mode)
            header = False
            mode = 'a'

rule get_mapped_random_manual:
    input:
        sam = 'resources/filtered-mappings/pre-filter_{read_filter}/manual/mapped/{sample}.sam',
    output:
        sam = 'resources/filtered-mappings/pre-filter_{read_filter}/manual/random/{sample}.sam'
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        import pandas as pd
        df = pd.read_csv(
        input.sam,
        sep="\t",
        )
        df = df.sample(frac=1).drop_duplicates("QNAME", keep="last")
        df.to_csv(output.sam,
                     sep = '\t',
                     index=False)


rule get_mapped_unique_manual:
    input:
        sam = 'resources/filtered-mappings/pre-filter_{read_filter}/manual/mapped/{sample}.sam',
    output:
        sam = 'resources/filtered-mappings/pre-filter_{read_filter}/manual/unique/{sample}.sam'
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        import pandas as pd
        df = pd.read_csv(
        input.sam,
        sep="\t",
        )
        df.drop_duplicates("QNAME", keep=False, inplace = True)
        df.to_csv(output.sam,
                     sep = '\t',
                     index=False)



rule get_coverage_per_manual_ref:
    input:
        unique = 'resources/filtered-mappings/pre-filter_{read_filter}/manual/unique/{sample}.sam',
        random = 'resources/filtered-mappings/pre-filter_{read_filter}/manual/random/{sample}.sam',
        all = 'resources/filtered-mappings/pre-filter_{read_filter}/manual/mapped/{sample}.sam'
    output:
        tsv = 'resources/coverage/pre-filter_{read_filter}/manual/per-ref/{sample}.tsv'
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        df_uni = pd.read_csv(
        input.unique,
        sep="\t",
        )
        df_uni = df_uni.groupby(['RNAME']).count()[['SEQ']]
        df_uni.rename(columns={"SEQ": "unique count"}, inplace = True)

        df_random =  pd.read_csv(
        input.random,
        sep="\t",
        )
        df_random = df_random.groupby(['RNAME']).count()[['SEQ']]
        df_random.rename(columns={"SEQ": "random count"}, inplace = True)

        df_all =  pd.read_csv(
        input.all,
        sep="\t",
        )
        df_all = df_all.groupby(['RNAME']).count()[['SEQ']]
        df_all.rename(columns={"SEQ": "all count"}, inplace = True)



        df = pd.concat([df_uni, df_random, df_all], axis=1)

        total_unique_read_count = df['unique count'].sum()
        total_random_read_count = df['random count'].sum()
        total_all_read_count = df['all count'].sum()
        df['unique fraction (norm all)'] = df.apply(lambda row: row['unique count']/total_random_read_count, axis =1)
        df['unique fraction'] = df.apply(lambda row: row['unique count']/total_unique_read_count, axis =1)
        df['random fraction'] = df.apply(lambda row: row['random count']/total_random_read_count, axis =1)
        df['all fraction'] = df.apply(lambda row: row['all count']/total_all_read_count, axis =1)
        df['unique/random'] = df.apply(lambda row: row['unique count']/row['random count'], axis =1)
        df['random/all'] = df.apply(lambda row: row['random count']/row['all count'], axis =1)
        df.reset_index(inplace = True)
        df['anticodon'] = df.apply(lambda row: anticodonfamily_from_rname(row['RNAME']), axis =1)
        df.sort_values('anticodon', inplace = True)
        df.set_index('RNAME', inplace = True)
        df.fillna(0, inplace = True)
        df.to_csv(output.tsv,
                         sep = '\t',
                         index=True)



rule get_coverage_summary_manual:
    input:
        expand('resources/coverage/pre-filter_{read_filter}/manual/per-ref/{sample}.tsv', sample=dm_samples, read_filter='{read_filter}')
    output:
        tsv = 'resources/coverage/pre-filter_{read_filter}/manual/coverage_summary_DM.tsv'
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        summary_df =  pd.DataFrame()
        for file in input:
            sample = file.split('/')[-1].replace('.tsv','')
            df = pd.read_csv(
            file,
            sep="\t",
            )
            df.set_index('RNAME', inplace = True)
            df.columns = [sample + ' ' +  column if column != 'anticodon' else column for column in df.columns]
            if len(summary_df) == 0:
                summary_df = df
            else:
                df.drop(['anticodon'], inplace = True, axis = 1)
                summary_df = pd.concat([summary_df, df], axis=1)
        summary_df.to_csv(output.tsv,
                         sep = '\t',
                         index=True)

rule get_coverage_summary_manual_all:
    input:
        expand('resources/coverage/pre-filter_{read_filter}/manual/per-ref/{sample}.tsv', sample=SAMPLES, read_filter='{read_filter}')
    output:
        tsv = 'resources/coverage/pre-filter_{read_filter}/manual/coverage_summary_all.tsv'
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        summary_df =  pd.DataFrame()
        for file in input:
            sample = file.split('/')[-1].replace('.tsv','')
            df = pd.read_csv(
            file,
            sep="\t",
            )
            df.set_index('RNAME', inplace = True)
            df.columns = [sample + ' ' +  column if column != 'anticodon' else column for column in df.columns]
            if len(summary_df) == 0:
                summary_df = df
            else:
                df.drop(['anticodon'], inplace = True, axis = 1)
                summary_df = pd.concat([summary_df, df], axis=1)
        summary_df.to_csv(output.tsv,
                         sep = '\t',
                         index=True)
