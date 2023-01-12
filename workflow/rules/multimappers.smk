
rule get_multimappers_per_sample:
    input:
        sam = 'resources/filtered-mappings/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/mapped/{sample}.sam'
    output:
        pair_tsv = 'resources/multimappers/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/per-ref/pair_{sample}.tsv'
    run:
        import pandas as pd
        df = pd.read_csv(input.sam, sep="\t")
        read_count = len(list(set(df['QNAME'].to_list())))
        normalization = 1000000/read_count
        df = df.groupby('RNAME').agg({'QNAME':list})
        refs = list(df.index)
        data = []
        for ref1 in refs:
            ref1_reads = df.at[ref1,'QNAME']
            for ref2 in refs:
                ref2_reads = df.at[ref2,'QNAME']
                shared_reads = list(set(ref1_reads) & set(ref2_reads))
                data.append({'ref1':ref1, 'ref2':ref2, 'reads1':len(ref1_reads), 'reads2':len(ref2_reads), 'shared':len(shared_reads)})
        multimapper_df = pd.DataFrame.from_records(data)
        for c in ['reads1', 'reads2', 'shared']:
            multimapper_df['norm_'+c] = multimapper_df[c]*normalization

        multimapper_df.to_csv(output.pair_tsv,sep = '\t',index=False)

rule get_multimappers_per_sample_cluster:
    input:
        sam = 'resources/filtered-mappings/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/mapped/{sample}.sam',
        clusters =  'resources/cluster/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-editdist-{dist}.yaml',
        cluster_names = 'resources/cluster/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusternames-editdist-{dist}.yaml', #TODO
    output:
        cluster_tsv = 'resources/multimappers/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/per-cluster-editdist-{dist}/pair_{sample}.tsv'
    run:
        import pandas as pd
        import yaml

        with open(input.clusters) as file:
            cluster_dict = yaml.safe_load(file)
        with open(input.cluster_names) as file:
            cluster_name_dict = yaml.safe_load(file)

        df = pd.read_csv(input.sam, sep="\t")
        read_count = len(list(set(df['QNAME'].to_list())))
        normalization = 1000000/read_count
        df['CLUSTER'] = df.apply(lambda row: cluster_dict[row['RNAME']] if row['RNAME'] in cluster_dict.keys() else 'Nan' , axis = 1)
        df = df.groupby('CLUSTER').agg({'QNAME':list})
        refs = list(df.index)
        data = []
        for ref1 in refs:
            ref1_reads = list(set(df.at[ref1,'QNAME']))
            for ref2 in refs:
                ref2_reads = list(set(df.at[ref2,'QNAME']))
                shared_reads = list(set(ref1_reads) & set(ref2_reads))
                data.append({'cluster1':ref1, 'cluster2':ref2, 'reads1':len(ref1_reads), 'reads2':len(ref2_reads), 'shared':len(shared_reads)})
        multimapper_df = pd.DataFrame.from_records(data)
        for c in ['reads1', 'reads2', 'shared']:
            multimapper_df['norm_'+c] = multimapper_df[c]*normalization
        multimapper_df['cluster1_multimapper_fraction'] =  multimapper_df.apply(lambda row: row['shared']/row['reads1'] if row['reads1']>0 else 0 , axis =1)
        multimapper_df['cluster1_name'] = multimapper_df.apply(lambda row: cluster_name_dict[row['cluster1']] if row['cluster1'] in cluster_name_dict.keys() else 'Nan',axis=1)
        multimapper_df['cluster2_name'] = multimapper_df.apply(lambda row: cluster_name_dict[row['cluster2']] if row['cluster2'] in cluster_name_dict.keys() else 'Nan',axis=1)

        multimapper_df.to_csv(output.cluster_tsv ,sep = '\t',index=False)

rule get_multimappers_summary_per_ref:
    input:
        expand('resources/multimappers/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/per-ref/pair_{sample}.tsv', sample=dm_samples)
    output:
        summary_tsv = 'resources/multimappers/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/multimapper_all.tsv',
        pivot_tsv = 'resources/multimappers/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/multimapper_pivot.tsv',
        cluster_map = 'resources/cluster/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/sch.fcluster.clusters.yaml',
        pdf = 'results/multimappers/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/multimapper.pdf'
    run:
        import pandas as pd
        import yaml
        import numpy as np
        import seaborn as sns
        import statistics
        import scipy.cluster.hierarchy as sch
        import matplotlib.pyplot as plt
        def get_cluster_sorted_id(df):
            X = df.to_numpy()
            X = np.nan_to_num(X, copy=True, nan=0.0, posinf=None, neginf=None)
            d = sch.distance.pdist(X)  # vector of ('55' choose 2) pairwise distances
            L = sch.linkage(d, method="complete")
            ind = sch.fcluster(L, 0.5 * d.max(), "distance")
            # get new colums order
            columns = [df.columns.tolist()[i] for i in list((np.argsort(ind)))]
            # get cluster dictionary
            cluster_dict = {}
            for i,ref in enumerate(df.columns.tolist()):
                cluster_dict[ref] = int(ind[i])
            return columns, cluster_dict

        summary_df =  pd.DataFrame()
        for file in input:
            sample = file.split('/')[-1].replace('.tsv','')
            df = pd.read_csv(
            file,
            sep="\t",
            )
            df.set_index(['ref1','ref2' ], inplace = True)
            df['ref1_multimaper_fraction'] = df.apply(lambda row: row['shared']/row['reads1'] if row['reads1']>0 else 0 , axis =1)
            df.columns = [sample + ' ' +  column  for column in df.columns]
            if len(summary_df) == 0:
                summary_df = df
            else:
                summary_df = pd.concat([summary_df, df], axis=1)

        norm_shared_columns = [c for c in summary_df.columns if 'norm_shared' in c]
        summary_df['max_norm_shared'] = summary_df.apply(lambda row: max([row[c] for c in norm_shared_columns]), axis =1 )
        summary_df['mean_norm_shared'] = summary_df.apply(lambda row: statistics.mean([row[c] for c in norm_shared_columns]), axis =1 )

        norm_ref1_columns = [c for c in summary_df.columns if 'norm_reads1' in c]
        summary_df['max_norm_ref1'] = summary_df.apply(lambda row: max([row[c] for c in norm_ref1_columns]), axis =1 )
        summary_df['mean_norm_ref1'] = summary_df.apply(lambda row: statistics.mean([row[c] for c in norm_ref1_columns]), axis =1 )

        norm_ref2_columns = [c for c in summary_df.columns if 'norm_reads2' in c]
        summary_df['max_norm_ref2'] = summary_df.apply(lambda row: max([row[c] for c in norm_ref2_columns]), axis =1 )
        summary_df['mean_norm_ref2'] = summary_df.apply(lambda row: statistics.mean([row[c] for c in norm_ref2_columns]), axis =1 )

        ref1_multimapper_fraction_columns = [c for c in summary_df.columns if 'ref1_multimaper_fraction' in c]
        summary_df['max_ref1_multimaper_fraction'] = summary_df.apply(lambda row: max([row[c] for c in ref1_multimapper_fraction_columns]), axis =1 )
        summary_df['mean_ref1_multimaper_fraction'] = summary_df.apply(lambda row: statistics.mean([row[c] for c in ref1_multimapper_fraction_columns]), axis =1 )

        summary_df['max_log_norm_shared'] = summary_df.apply(lambda row: np.log10(row['max_norm_shared']) if row['max_norm_shared']>1 else 0, axis=1 )
        summary_df['inverse max_log_norm_shared'] = summary_df.apply(lambda row: 1/row['max_log_norm_shared'] if row['max_log_norm_shared']>0 else 1, axis=1 )

        summary_df['min_norm_shared'] = summary_df.apply(lambda row: min([row[c] for c in norm_shared_columns]), axis =1 )
        summary_df['min_log_norm_shared'] = summary_df.apply(lambda row: np.log10(row['min_norm_shared']) if row['min_norm_shared']>1 else 0, axis=1 )
        summary_df.to_csv(output.summary_tsv,
                         sep = '\t',
                         index=True)

        # get pivot table and cluster mapping
        pivot_df = summary_df.reset_index().pivot_table(values = 'max_log_norm_shared', columns = 'ref1', index = 'ref2', fill_value=1)
        cluster_sorted_ref_ids, ref_cluster_map = get_cluster_sorted_id(pivot_df)
        pivot_df = pivot_df.reindex(cluster_sorted_ref_ids, axis=1)
        pivot_df = pivot_df.reindex(cluster_sorted_ref_ids, axis=0)
        pivot_df.to_csv(output.pivot_tsv,
                         sep = '\t',
                         index=True)

        with open(output.cluster_map, 'w') as file:
            outputs = yaml.dump(ref_cluster_map, file)


        # plot multimapper heatmap
        fig, ax = plt.subplots(figsize=(70, 70))
        sns.heatmap(pivot_df, ax = ax, square = True,linewidths=0.00005,)
        fig.savefig(output.pdf, bbox_inches="tight")

rule plot_multimappers_between_to_level_clusters:
    input:
        expand('resources/multimappers/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/per-cluster-ed-{e_cutoff}-mm-{m_cutoff}_DM/pair_{sample}.tsv', sample=SAMPLES, e_cutoff='{e_cutoff}', m_cutoff='{m_cutoff}')
    output:
        pdf = 'results/multimappers/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/cluster-ed-{e_cutoff}-mm-{m_cutoff}_DM/{treatment}/multimappers_between_clusters.pdf',
        tsv = 'resources/multimappers/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/cluster-ed-{e_cutoff}-mm-{m_cutoff}_DM/{treatment}/multimappers_between_clusters.tsv'
    wildcard_constraints:
        e_cutoff="[0-5]",
        m_cutoff="\d+",
    run:
        import yaml
        import os
        import numpy as np
        import seaborn as sns
        import statistics
        import matplotlib.pyplot as plt


        summary_df =  pd.DataFrame()
        for file in input:
            sample = file.split('/')[-1].replace('.tsv','').replace('pair_', '')
            if sample_dict[sample]['treatment']!=wildcards.treatment:
                continue
            sample_for_pivot = sample
            df = pd.read_csv(
            file,
            sep="\t",
            )
            df.set_index(['cluster1','cluster2' ], inplace = True)
            df.columns = [ sample + ' ' +  column for column in df.columns ]
            if len(summary_df) == 0:

                summary_df = df
            else:
                summary_df = pd.concat([summary_df, df], axis=1)

        sample = sample_for_pivot
        #summary_df.reset_index(inplace = True)
        norm_shared_columns = [c for c in summary_df.columns if 'norm_shared' in c]
        summary_df['max_norm_shared'] = summary_df.apply(lambda row: max([row[c] for c in norm_shared_columns]), axis =1 )
        summary_df['mean_norm_shared'] = summary_df.apply(lambda row: statistics.mean([row[c] for c in norm_shared_columns]), axis =1 )
        cluster1_multimaper_fraction_columns = [c for c in summary_df.columns if 'cluster1_multimapper_fraction' in c]

        summary_df['max_ref1_multimaper_fraction'] = summary_df.apply(lambda row: max([row[c] for c in cluster1_multimaper_fraction_columns]), axis =1 )
        summary_df['mean_ref1_multimaper_fraction'] = summary_df.apply(lambda row: statistics.mean([row[c] for c in cluster1_multimaper_fraction_columns]), axis =1 )
        summary_df['min_ref1_multimaper_fraction'] = summary_df.apply(lambda row: min([row[c] for c in cluster1_multimaper_fraction_columns]), axis =1 )


        fig, axs = plt.subplots(figsize=(240*CM, 60*CM),nrows=1, ncols=4, layout='tight')
        pivot_df = summary_df.reset_index().pivot_table(values = 'max_norm_shared', columns =  sample + ' cluster1_name', index = sample + ' cluster2_name', fill_value=0)
        sns.heatmap(pivot_df, ax = axs[0], square = True,linewidths=0.1, cmap='Blues', linecolor='black')
        axs[0].set_title('max RPM shared')
        pivot_df = summary_df.reset_index().pivot_table(values = 'mean_norm_shared', columns =  sample + ' cluster1_name', index = sample + ' cluster2_name', fill_value=0)
        sns.heatmap(pivot_df, ax = axs[1], square = True,linewidths=0.1, cmap='Blues', linecolor='black')
        axs[1].set_title('mean RPM shared')
        pivot_df = summary_df.reset_index().pivot_table(values = 'max_ref1_multimaper_fraction', columns =  sample + ' cluster1_name', index = sample + ' cluster2_name', fill_value=0)
        sns.heatmap(pivot_df, ax = axs[2], square = True,linewidths=0.1, cmap='Blues', linecolor='black')
        axs[2].set_title('max ref1 fraction shared')
        pivot_df = summary_df.reset_index().pivot_table(values = 'mean_ref1_multimaper_fraction', columns =  sample + ' cluster1_name', index = sample + ' cluster2_name', fill_value=0)
        sns.heatmap(pivot_df, ax = axs[3], square = True,linewidths=0.1, cmap='Blues', linecolor='black')
        axs[3].set_title('mean ref1 fraction shared')
        fig.savefig(output.pdf)
        summary_df.to_csv(output.tsv,
                         sep = '\t',
                         index=True)

rule plot_multimappers_between_clusters:
    input:
        expand('resources/multimappers/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/per-cluster-editdist-{dist}/pair_{sample}.tsv', sample=SAMPLES, dist='{dist}')
    output:
        pdf = 'results/multimappers/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/cluster-editdist-{dist}/{treatment}/multimappers_between_clusters.pdf',
        tsv = 'resources/multimappers/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/cluster-editdist-{dist}/{treatment}/multimappers_between_clusters.tsv'
    wildcard_constraints:
        dist="[0-5]",
    run:
        import yaml
        import os
        import numpy as np
        import seaborn as sns
        import statistics
        import matplotlib.pyplot as plt


        summary_df =  pd.DataFrame()
        for file in input:
            sample = file.split('/')[-1].replace('.tsv','').replace('pair_', '')
            if sample_dict[sample]['treatment']!=wildcards.treatment:
                continue
            sample_for_pivot = sample
            df = pd.read_csv(
            file,
            sep="\t",
            )
            df.set_index(['cluster1','cluster2' ], inplace = True)
            df.columns = [ sample + ' ' +  column for column in df.columns ]
            if len(summary_df) == 0:

                summary_df = df
            else:
                summary_df = pd.concat([summary_df, df], axis=1)

        sample = sample_for_pivot
        #summary_df.reset_index(inplace = True)
        norm_shared_columns = [c for c in summary_df.columns if 'norm_shared' in c]
        summary_df['max_norm_shared'] = summary_df.apply(lambda row: max([row[c] for c in norm_shared_columns]), axis =1 )
        summary_df['mean_norm_shared'] = summary_df.apply(lambda row: statistics.mean([row[c] for c in norm_shared_columns]), axis =1 )
        cluster1_multimaper_fraction_columns = [c for c in summary_df.columns if 'cluster1_multimapper_fraction' in c]

        summary_df['max_ref1_multimaper_fraction'] = summary_df.apply(lambda row: max([row[c] for c in cluster1_multimaper_fraction_columns]), axis =1 )
        summary_df['mean_ref1_multimaper_fraction'] = summary_df.apply(lambda row: statistics.mean([row[c] for c in cluster1_multimaper_fraction_columns]), axis =1 )
        summary_df['min_ref1_multimaper_fraction'] = summary_df.apply(lambda row: min([row[c] for c in cluster1_multimaper_fraction_columns]), axis =1 )


        fig, axs = plt.subplots(figsize=(240*CM, 60*CM),nrows=1, ncols=4, layout='tight')
        pivot_df = summary_df.reset_index().pivot_table(values = 'max_norm_shared', columns =  sample + ' cluster1_name', index = sample + ' cluster2_name', fill_value=0)
        sns.heatmap(pivot_df, ax = axs[0], square = True,linewidths=0.1, cmap='Blues', linecolor='black')
        axs[0].set_title('max RPM shared')
        pivot_df = summary_df.reset_index().pivot_table(values = 'mean_norm_shared', columns =  sample + ' cluster1_name', index = sample + ' cluster2_name', fill_value=0)
        sns.heatmap(pivot_df, ax = axs[1], square = True,linewidths=0.1, cmap='Blues', linecolor='black')
        axs[1].set_title('mean RPM shared')
        pivot_df = summary_df.reset_index().pivot_table(values = 'max_ref1_multimaper_fraction', columns =  sample + ' cluster1_name', index = sample + ' cluster2_name', fill_value=0)
        sns.heatmap(pivot_df, ax = axs[2], square = True,linewidths=0.1, cmap='Blues', linecolor='black')
        axs[2].set_title('max ref1 fraction shared')
        pivot_df = summary_df.reset_index().pivot_table(values = 'mean_ref1_multimaper_fraction', columns =  sample + ' cluster1_name', index = sample + ' cluster2_name', fill_value=0)
        sns.heatmap(pivot_df, ax = axs[3], square = True,linewidths=0.1, cmap='Blues', linecolor='black')
        axs[3].set_title('mean ref1 fraction shared')
        fig.savefig(output.pdf)
        summary_df.to_csv(output.tsv,
                         sep = '\t',
                         index=True)
