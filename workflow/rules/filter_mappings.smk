
# Author: Maria Waldl • code@waldl.org
# Version: 2024-01-24

rule get_raw_mapped_and_remove_BS_GAtransitions:
    input:
        sam = 'resources/mapping/all/polyTtrimming_ONtttt/{sample}.sam'
    output:
        sam = 'resources/filtered-mappings/pre-filter_none/raw/mapped/{sample}.sam'
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



rule filter_cigar:
    input:
        sam = 'resources/filtered-mappings/pre-filter_none/raw/mapped/{sample}.sam'
    output:
        sam = 'resources/filtered-mappings/pre-filter_cigar/raw/mapped/{sample}.sam'
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        import pandas as pd
        df = pd.read_csv(
        input.sam,
        sep="\t",
        )
        header = True
        for read, g_df in df.groupby('QNAME'):
            cigars = list(set(g_df['CIGAR'].to_list()))
            if len(cigars)>1:
                continue

            mode = 'a'
            if header:
                mode = 'w'
            g_df.to_csv(output.sam,
                         sep = '\t',
                         index=False,
                         header=header,
                         mode='a')
            header = False

rule filter_umi:
    input:
        sam = 'resources/filtered-mappings/pre-filter_none/raw/mapped/{sample}.sam'
    output:
        sam = 'resources/filtered-mappings/pre-filter_umi/raw/mapped/{sample}.sam'
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        import pandas as pd
        df = pd.read_csv(
        input.sam,
        sep="\t",
        )
        header = True
        i=0
        for read, g_df in df.groupby('SEQ'):
            g_df['umi'] = g_df.apply(lambda row: row['QNAME'].split('_')[-1], axis =1)
            g_df.drop_duplicates(subset=['umi', 'RNAME'], inplace = True, keep= 'first')
            mode = 'a'
            if header:
                mode = 'w'
            g_df.to_csv(output.sam,
                         sep = '\t',
                         index=False,
                         header=header,
                         mode='a')
            header = False

rule filter_umicigar:
    input:
        sam =  'resources/filtered-mappings/pre-filter_umi/raw/mapped/{sample}.sam'
    output:
        sam = 'resources/filtered-mappings/pre-filter_umicigar/raw/mapped/{sample}.sam'
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        import pandas as pd
        df = pd.read_csv(
        input.sam,
        sep="\t",
        )
        header = True
        for read, g_df in df.groupby('QNAME'):
            cigars = list(set(g_df['CIGAR'].to_list()))
            if len(cigars)>1:
                continue

            mode = 'a'
            if header:
                mode = 'w'
            g_df.to_csv(output.sam,
                         sep = '\t',
                         index=False,
                         header=header,
                         mode='a')
            header = False


rule get_uniquely_mapped:
    input:
        sam = 'resources/filtered-mappings/pre-filter_{reads_filter}/raw/mapped/{sample}.sam'
    output:
        sam = 'resources/filtered-mappings/pre-filter_{reads_filter}/raw/unique/{sample}.sam'
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

rule get_random_mapped:
    input:
        sam = 'resources/filtered-mappings/pre-filter_{reads_filter}/raw/mapped/{sample}.sam'
    output:
        sam = 'resources/filtered-mappings/pre-filter_{reads_filter}/raw/random/{sample}.sam'
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

rule get_coverage_per_ref:
    input:
        unique = 'resources/filtered-mappings/pre-filter_{reads_filter}/raw/unique/{sample}.sam',
        random = 'resources/filtered-mappings/pre-filter_{reads_filter}/raw/random/{sample}.sam',
        all = 'resources/filtered-mappings/pre-filter_{reads_filter}/raw/mapped/{sample}.sam'
    output:
        tsv = 'resources/coverage/pre-filter_{reads_filter}/raw/per-ref/{sample}.tsv'
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



rule max_coverage_mapped:
    input:
        sam = 'resources/filtered-mappings/pre-filter_umi/raw/mapped/{sample}.sam',
        raw_abundance = 'resources/coverage/pre-filter_umi/raw/min_coverage_summary_DM.tsv'
    output:
        sam = 'resources/filtered-mappings/pre-filter_umi/raw/max-coverage_cutoff-{cutoff}_high-confidence-{hc}/{sample}.sam',
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        import pandas as pd

        #raw_abundance_dict[row['RNAME']]['raw']
        raw_abundance_dict = pd.read_csv( input.raw_abundance, sep="\t" ).set_index('RNAME').to_dict(orient='index')


        df = pd.read_csv(
        input.sam,
        sep="\t",
        #nrows=100000
        )
        header = True
        mode = 'w'
        for query, gdf in df.groupby("QNAME"):
            gdf = gdf.copy()
            if len(gdf) ==1:
                gdf.to_csv(output.sam, sep = '\t', mode=mode, index=False, header=header)
                header = False
                mode = 'a'
                continue
            gdf['hc'] = gdf.apply(lambda row: True if '[H]' in row['RNAME']else False, axis =1)
            if gdf['hc'].sum() >= 1:
                gdf= gdf[gdf['hc']==True]
            gdf['raw_abundance'] = gdf.apply(lambda row: raw_abundance_dict[row['RNAME']]['max all fraction'], axis =1)
            #gdf['unique_abundance'] = gdf.apply(lambda row: raw_abundance_dict[row['RNAME']]['raw'], axis =1)
            gdf.sort_values(by= 'raw_abundance', ascending=False,inplace = True)
            gdf = gdf.iloc[:1,:]
            gdf.drop(columns = ['hc','raw_abundance' ], inplace = True)
            gdf.to_csv(output.sam, sep = '\t', mode=mode, index=False, header=header)
            header = False
            mode = 'a'


rule get_coverage_per_selected_ref_with_max:
    input:
        unique = 'resources/filtered-mappings/pre-filter_{reads_filter}/raw/unique/{sample}.sam',
        random = 'resources/filtered-mappings/pre-filter_{reads_filter}/raw/random/{sample}.sam',
        all = 'resources/filtered-mappings/pre-filter_{reads_filter}/raw/mapped/{sample}.sam',
        max_cov_hc = 'resources/filtered-mappings/pre-filter_{reads_filter}/raw/max-coverage_cutoff-1_high-confidence-on/{sample}.sam',
    output:
        tsv = 'resources/coverage/pre-filter_{reads_filter}/raw/per-ref-with-max/{sample}.tsv'
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

        df_max =  pd.read_csv(
        input.max_cov_hc,
        sep="\t",
        )
        df_max = df_max.groupby(['RNAME']).count()[['SEQ']]
        df_max.rename(columns={"SEQ": "max count"}, inplace = True)

        df = pd.concat([df_uni, df_random, df_all, df_max], axis=1)

        total_unique_read_count = df['unique count'].sum()
        total_random_read_count = df['random count'].sum()
        total_all_read_count = df['all count'].sum()
        total_max_count = df['max count'].sum()
        df['unique fraction'] = df.apply(lambda row: row['unique count']/total_unique_read_count, axis =1)
        df['random fraction'] = df.apply(lambda row: row['random count']/total_random_read_count, axis =1)
        df['all fraction'] = df.apply(lambda row: row['all count']/total_all_read_count, axis =1)
        df['max fraction'] = df.apply(lambda row: row['max count']/total_max_count, axis =1)
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

rule get_coverage_summary:
    input:
        expand('resources/coverage/pre-filter_{reads_filter}/raw/per-ref/{sample}.tsv', reads_filter='{reads_filter}',sample=dm_samples)
    output:
        tsv = 'resources/coverage/pre-filter_{reads_filter}/raw/coverage_summary_DM.tsv'
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

rule get_coverage_summary_with_max:
    input:
        expand('resources/coverage/pre-filter_{reads_filter}/raw/per-ref-with-max/{sample}.tsv', reads_filter='{reads_filter}', sample=dm_samples)
    output:
        tsv = 'resources/coverage/pre-filter_{reads_filter}/raw/coverage_summary_DM-with-max.tsv'
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

rule get_min_coverage_summary:
    input:
        'resources/coverage/pre-filter_{reads_filter}/raw/coverage_summary_DM.tsv'
    output:
        'resources/coverage/pre-filter_{reads_filter}/raw/min_coverage_summary_DM.tsv'
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        df = pd.read_csv(input[0], sep="\t")
        df.set_index('RNAME', inplace = True)
        criteria = ['unique count', 'random count', 'all count', 'unique fraction', 'random fraction', 'all fraction', 'unique/random', 'random/all']
        for score in criteria:
            columns = [col for col in list(df.columns) if score in col]
            df['min '+score] = df.apply(lambda row: min([row[col] for col in columns ]), axis=1)
            df['max '+score] = df.apply(lambda row: max([row[col] for col in columns ]), axis=1)
            df.drop(columns, axis =1, inplace = True)
        df.to_csv(output[0], sep = '\t', index=True)

rule get_min_coverage_summary_with_max:
    input:
        'resources/coverage/pre-filter_{reads_filter}/raw/coverage_summary_DM-with-max.tsv'
    output:
        'resources/coverage/pre-filter_{reads_filter}/raw/min_coverage_summary_DM-with-max.tsv'
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        df = pd.read_csv(input[0], sep="\t")
        df.set_index('RNAME', inplace = True)
        criteria = ['unique count', 'random count', 'all count' ,'max count', 'unique fraction', 'random fraction', 'all fraction','max fraction', 'unique/random', 'random/all']
        for score in criteria:
            columns = [col for col in list(df.columns) if score in col]
            df['min '+score] = df.apply(lambda row: min([row[col] for col in columns ]), axis=1)
            df['max '+score] = df.apply(lambda row: max([row[col] for col in columns ]), axis=1)
            df.drop(columns, axis =1, inplace = True)
        df.to_csv(output[0], sep = '\t', index=True)



rule get_min_cov_refs:
    input:
        cov_sum = 'resources/coverage/pre-filter_{reads_filter}/raw/min_coverage_summary_DM-with-max.tsv'
    output:
        keep_refs = 'resources/min_coverage_refs/pre-filter_{reads_filter}/min_cov_refs.yaml'
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        import pandas as pd
        import yaml
        df = pd.read_csv(input.cov_sum, sep="\t")
        refs = []

        for criterium, cutoff in config['min_coverage_per_ref']:
            selected = df[df[criterium]>=cutoff]['RNAME'].to_list()
            print(criterium, cutoff, len(selected))
            refs += selected
        refs = list(set(refs))
        with open(output.keep_refs, 'w') as file:
            outputs = yaml.dump(refs, file)
        print(len(df))
        print( len(refs))

rule get_raw_mapped_selected:
    input:
        #sam = 'resources/mapping/{sample}.sam',
        sam = 'resources/mapping/selected/polyTtrimming_ONtttt/{sample}.sam'
    output:
        sam = 'resources/filtered-mappings/pre-filter_none/selected/mapped/{sample}.sam'
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        import pandas as pd
        import numpy as np
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
        df.to_csv(output.sam, sep="\t", index=False)

rule filter_umi_on_selected:
    input:
        sam = 'resources/filtered-mappings/pre-filter_none/selected/mapped/{sample}.sam'
    output:
        sam = 'resources/filtered-mappings/pre-filter_umi/selected/mapped/{sample}.sam'
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

rule get_mapped_min_cov_ref_random:
    input:
        sam = 'resources/filtered-mappings/pre-filter_{reads_filter}/selected/mapped/{sample}.sam',
    output:
        sam = 'resources/filtered-mappings/pre-filter_{reads_filter}/selected/random/{sample}.sam'
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


rule get_mapped_min_cov_ref_unique:
    input:
        sam = 'resources/filtered-mappings/pre-filter_{reads_filter}/selected/mapped/{sample}.sam',
    output:
        sam = 'resources/filtered-mappings/pre-filter_{reads_filter}/selected/unique/{sample}.sam'
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



rule get_coverage_per_selected_ref:
    input:
        unique = 'resources/filtered-mappings/pre-filter_{reads_filter}/selected/unique/{sample}.sam',
        random = 'resources/filtered-mappings/pre-filter_{reads_filter}/selected/random/{sample}.sam',
        all = 'resources/filtered-mappings/pre-filter_{reads_filter}/selected/mapped/{sample}.sam'
    output:
        tsv = 'resources/coverage/pre-filter_{reads_filter}/selected/per-ref/{sample}.tsv'
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



rule get_coverage_summary_selected:
    input:
        expand('resources/coverage/pre-filter_{reads_filter}/selected/per-ref/{sample}.tsv', reads_filter='{reads_filter}',sample=dm_samples)
    output:
        tsv = 'resources/coverage/pre-filter_{reads_filter}/selected/coverage_summary_DM.tsv'
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

rule get_coverage_summary_selected_mock:
    input:
        expand('resources/coverage/pre-filter_{reads_filter}/selected/per-ref/{sample}.tsv', reads_filter='{reads_filter}',sample=SAMPLES)
    output:
        tsv = 'resources/coverage/pre-filter_{reads_filter}/selected/coverage_summary_all.tsv'
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


rule get_min_coverage_summary_selected:
    input:
        'resources/coverage/pre-filter_{reads_filter}/{ref_set}/coverage_summary_DM.tsv'
    output:
        'resources/coverage/pre-filter_{reads_filter}/{ref_set}/min_coverage_summary_DM.tsv'
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        df = pd.read_csv(input[0], sep="\t")
        df.set_index('RNAME', inplace = True)
        criteria = ['unique count', 'random count', 'all count', 'unique fraction', 'random fraction', 'all fraction', 'unique/random', 'random/all']
        for score in criteria:
            columns = [col for col in list(df.columns) if score in col]
            df['min '+score] = df.apply(lambda row: min([row[col] for col in columns ]), axis=1)
            df['max '+score] = df.apply(lambda row: max([row[col] for col in columns ]), axis=1)
            df.drop(columns, axis =1, inplace = True)
        df.to_csv(output[0], sep = '\t', index=True)
