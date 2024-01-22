rule get_abundance_summary:
    input:
        tsv = 'resources/coverage/pre-filter_{reads_filter}/{ref_set}/coverage_summary_{treatment}.tsv',
        cluster_map = 'resources/cluster/pre-filter_{reads_filter}/{ref_set}/clusters-ed-{e_cutoff}-mm-{m_cutoff}_{c_treatment}.yaml',
        cluster_names = 'resources/cluster/pre-filter_{reads_filter}/{ref_set}/anticodonbased_clusternames-ed-{e_cutoff}-mm-{m_cutoff}_{c_treatment}.yaml'
    output:
        tsv = 'resources/coverage/pre-filter_{reads_filter}/{ref_set}/clusters-ed-{e_cutoff}-mm-{m_cutoff}_{c_treatment}_coverage_summary_{treatment}.tsv',
        html = 'resources/coverage/pre-filter_{reads_filter}/{ref_set}/clusters-ed-{e_cutoff}-mm-{m_cutoff}_{c_treatment}_coverage_summary_{treatment}.html'
    run:
        import pandas as pd
        import yaml
        import numpy as np

        with open(input.cluster_map) as file:
            cluster_dict = yaml.safe_load(file)

        with open(input.cluster_names) as file:
            cluster_name_dict = yaml.safe_load(file)

        df = pd.read_csv(
        input.tsv,
        sep="\t",
        index_col = 'RNAME'
        )

        # only keep random fraction as abundance score
        # rename 'random fraction' to just 'fraction'
        columns_to_keep =  [c for c in df.columns if 'random fraction' in c]
        for c in df.columns:
            if c not in columns_to_keep:
                df.drop(columns = [c], inplace = True)
            else:
                df.rename(columns = {c:c.replace('random ', '')}, inplace=True)

        # get list of samples and their corresponding time point
        samples = [c.split(' ')[0] for c in list(df.columns)]
        samples_timepoints = [ [sample, sample_dict[sample]['timepoint']] for sample in samples]

        # get mean read fraction per timepoint
        timepoints = list(set([tp for [s,tp] in samples_timepoints]))
        for timepoint in timepoints:
            per_timepoint_col = f't{timepoint:02d} fraction'
            sample_cols = [ f'{s} fraction' for [s,tp] in samples_timepoints if tp ==timepoint]
            df[per_timepoint_col] = df.loc[:, sample_cols].mean(axis = 1)

        # get overall mean, min and max fraction as mean of all timepoints
        timpoint_fraction_cols = [ f't{timepoint:02d} fraction' for timepoint in timepoints]
        print(timpoint_fraction_cols)
        df['mean fraction'] = df.loc[:, timpoint_fraction_cols].mean(axis = 1)
        df['min fraction'] = df.loc[:, timpoint_fraction_cols].min(axis = 1)
        df['max fraction'] = df.loc[:, timpoint_fraction_cols].max(axis = 1)


        df.reset_index(inplace = True)

        # annotate cluster
        df['cluster id'] = df.apply(lambda row: cluster_dict[row['RNAME']], axis = 1 )
        df['anticodon ratio'] = df.apply(lambda row: cluster_name_dict[row['cluster id']].split('[')[1].split(']')[0], axis = 1 )
        df['cluster name'] = df.apply(lambda row: cluster_name_dict[row['cluster id']].split('[')[0], axis = 1 )

        # annotate anticodon
        df['anticodon'] = df.apply(lambda row: anticodonfamily_from_rname(row['RNAME']), axis = 1)
        df['RNAME'] = df.apply(lambda row: row['RNAME'].split('(')[0], axis = 1)


        df_cluster = df.groupby(['cluster id', 'cluster name']).sum()
        df_cluster.sort_values(by= ['cluster name', 'cluster id'], inplace = True)


        # get mean read fraction per timepoint
        for timepoint in timepoints:
            per_timepoint_col = f't{timepoint:02d} fraction'
            sample_cols = [ f'{s} fraction' for [s,tp] in samples_timepoints if tp ==timepoint]
            df_cluster[per_timepoint_col] = df_cluster.loc[:, sample_cols].mean(axis = 1)

        # get overall mean, min and max fraction as mean of all timepoints
        timpoint_fraction_cols = [ f't{timepoint:02d} fraction' for timepoint in timepoints]
        print(timpoint_fraction_cols)
        df_cluster['mean fraction'] = df_cluster.loc[:, timpoint_fraction_cols].mean(axis = 1)
        df_cluster['min fraction'] = df_cluster.loc[:, timpoint_fraction_cols].min(axis = 1)
        df_cluster['max fraction'] = df_cluster.loc[:, timpoint_fraction_cols].max(axis = 1)

        clustered_cols = [ c + ' clustered' for c in df_cluster.columns]

        df_merged = pd.merge(df_cluster, df, on= ['cluster name', 'cluster id'], how = 'outer', suffixes = (' clustered', ''))
        df_merged.set_index(['cluster name', 'cluster id', 'anticodon ratio','anticodon', 'RNAME'],inplace = True)
        df_merged.to_csv(output.tsv, sep = '\t')

        df_merged.reset_index(inplace = True)
        df_merged.set_index(['cluster name', 'cluster id', 'anticodon ratio']+clustered_cols+ ['anticodon', 'RNAME'],inplace = True)
        df_merged.to_html(output.html)


rule get_abundance_summary_all:
    input:
        tsv = 'resources/coverage/pre-filter_{reads_filter}/{ref_set}/coverage_summary_all.tsv',
        cluster_map = 'resources/cluster/pre-filter_{reads_filter}/{ref_set}/clusters-ed-{e_cutoff}-mm-{m_cutoff}_{c_treatment}.yaml',
        cluster_names = 'resources/cluster/pre-filter_{reads_filter}/{ref_set}/anticodonbased_clusternames-ed-{e_cutoff}-mm-{m_cutoff}_{c_treatment}.yaml'
    output:
        tsv = 'resources/coverage/pre-filter_{reads_filter}/{ref_set}/clusters-ed-{e_cutoff}-mm-{m_cutoff}_{c_treatment}_coverage_summary_all.tsv'
    run:
        import pandas as pd
        import yaml
        import numpy as np

        with open(input.cluster_map) as file:
            cluster_dict = yaml.safe_load(file)

        with open(input.cluster_names) as file:
            cluster_name_dict = yaml.safe_load(file)

        df = pd.read_csv(
        input.tsv,
        sep="\t",
        index_col = 'RNAME'
        )

        # only keep random fraction as abundance score
        # rename 'random fraction' to just 'fraction'
        columns_to_keep =  [c for c in df.columns if 'random fraction' in c]
        for c in df.columns:
            if c not in columns_to_keep:
                df.drop(columns = [c], inplace = True)
            else:
                df.rename(columns = {c:c.replace('random ', '')}, inplace=True)

        # get list of samples and their corresponding time point
        samples = [c.split(' ')[0] for c in list(df.columns)]
        samples_timepoints = [ [sample, sample_dict[sample]['timepoint']] for sample in samples]

        # get mean read fraction per timepoint
        timepoints = list(set([tp for [s,tp] in samples_timepoints]))
        for timepoint in timepoints:
            per_timepoint_col = f't{timepoint:02d} fraction'
            sample_cols = [ f'{s} fraction' for [s,tp] in samples_timepoints if tp ==timepoint]
            df[per_timepoint_col] = df.loc[:, sample_cols].mean(axis = 1)

        # get overall mean, min and max fraction as mean of all timepoints
        timpoint_fraction_cols = [ f't{timepoint:02d} fraction' for timepoint in timepoints]
        print(timpoint_fraction_cols)
        df['mean fraction'] = df.loc[:, timpoint_fraction_cols].mean(axis = 1)
        df['min fraction'] = df.loc[:, timpoint_fraction_cols].min(axis = 1)
        df['max fraction'] = df.loc[:, timpoint_fraction_cols].max(axis = 1)


        df.reset_index(inplace = True)

        # annotate cluster
        df['cluster'] = df.apply(lambda row: cluster_dict[row['RNAME']], axis = 1 )
        df['cluster name'] = df.apply(lambda row: cluster_name_dict[row['cluster']], axis = 1 )

        # annotate anticodon
        df['anticodon'] = df.apply(lambda row: anticodonfamily_from_rname(row['RNAME']), axis = 1)

        df.set_index(['cluster name', 'anticodon', 'RNAME'],inplace = True)
        print(df.head(40))
        df.to_csv(output.tsv, sep = '\t')
