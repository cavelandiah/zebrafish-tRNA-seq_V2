# Author: Maria Waldl • code@waldl.org
# Version: 2024-01-24


import seaborn as sns
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
import numpy as np
import pandas as pd
import statistics
import yaml


def plot_multimappers_between_clusters(df, figure_path, title):
    # Author: Maria Waldl • code@waldl.org
    # Version: 2024-01-24

    fig, axs = plt.subplots(
        figsize=(240 * CM, 60 * CM), nrows=1, ncols=4, layout="tight"
    )
    fig.suptitle(title)

    pivot_df = df.pivot_table(
        values="max_norm_shared",
        index="cluster1 name",
        columns="cluster2 name",
        fill_value=0,
    )
    sns.heatmap(
        pivot_df,
        ax=axs[0],
        square=True,
        linewidths=0.1,
        cmap="Blues",
        linecolor="black",
    )
    axs[0].set_title("number of muti-mapping reads in RPM\n(maximum over all samples)")

    pivot_df = df.pivot_table(
        values="mean_norm_shared",
        index="cluster1 name",
        columns="cluster2 name",
        fill_value=0,
    )
    sns.heatmap(
        pivot_df,
        ax=axs[1],
        square=True,
        linewidths=0.1,
        cmap="Blues",
        linecolor="black",
    )
    axs[1].set_title("number of muti-mapping reads in RPM\n(mean over all samples)")

    df["max_ref1_multimaper_percentage"] = df["max_ref1_multimaper_fraction"] * 100
    pivot_df = df.pivot_table(
        values="max_ref1_multimaper_fraction",
        index="cluster1 name",
        columns="cluster2 name",
        fill_value=0,
    )
    sns.heatmap(
        pivot_df,
        ax=axs[2],
        square=True,
        linewidths=0.1,
        cmap="Blues",
        linecolor="black",
    )
    axs[2].set_title(
        "percentage of reads from cluster 1 that also mapp to cluster 2\n(RPM on cluster 1 = 100%; maximum over all samples)"
    )

    df["mean_ref1_multimaper_percentage"] = df["mean_ref1_multimaper_fraction"] * 100
    pivot_df = df.pivot_table(
        values="mean_ref1_multimaper_fraction",
        index="cluster1 name",
        columns="cluster2 name",
        fill_value=0,
    )
    sns.heatmap(
        pivot_df,
        ax=axs[3],
        square=True,
        linewidths=0.1,
        cmap="Blues",
        linecolor="black",
    )
    axs[3].set_title(
        "percentage of reads from cluster 1 that also mapp to cluster 2\n(RPM on cluster 1 = 100%; mean over all samples)"
    )

    for ax in axs.flat:
        ax.set_xlabel("cluster 2")
        ax.set_ylabel("cluster 1")

    fig.savefig(figure_path)


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
    for i, ref in enumerate(df.columns.tolist()):
        cluster_dict[ref] = int(ind[i])
    return columns, cluster_dict


def plot_multimappers_between_refs(df, figure_path, pivot_tsv, cluster_map):

    # get pivot table and cluster mapping
    pivot_df = df.reset_index().pivot_table(
        values="max_log_norm_shared", columns="ref1", index="ref2", fill_value=1
    )
    cluster_sorted_ref_ids, ref_cluster_map = get_cluster_sorted_id(pivot_df)
    pivot_df = pivot_df.reindex(cluster_sorted_ref_ids, axis=1)
    pivot_df = pivot_df.reindex(cluster_sorted_ref_ids, axis=0)
    pivot_df.to_csv(pivot_tsv, sep="\t", index=True)

    # save cluster info
    with open(cluster_map, "w") as file:
        outputs = yaml.dump(ref_cluster_map, file)

    # plot multimapper heatmap
    fig, ax = plt.subplots(figsize=(70, 70))
    sns.heatmap(
        pivot_df,
        ax=ax,
        square=True,
        linewidths=0.00005,
    )
    ax.set_title(
        "Shared RPM between references on a log10 scale (maximum over all DM samples)"
    )
    fig.savefig(figure_path, bbox_inches="tight")


def summarize_multimappers_between_clusters(tsv_files, cluster_name_file, summary_file):
    # Author: Maria Waldl • code@waldl.org
    # Version: 2024-01-24

    with open(cluster_name_file) as file:
        cluster_name_dict = yaml.safe_load(file)

    summary_df = pd.DataFrame()
    for file in tsv_files:
        sample = file.split("/")[-1].replace(".tsv", "").replace("pair_", "")
        df = pd.read_csv(
            file,
            sep="\t",
        )
        df.set_index(["cluster1", "cluster2"], inplace=True)
        df.drop(columns=["cluster1 name", "cluster2 name"], inplace=True)
        df.columns = [sample + " " + column for column in df.columns]
        if len(summary_df) == 0:
            summary_df = df
        else:
            summary_df = pd.concat([summary_df, df], axis=1)

    summary_df.reset_index(inplace=True)
    summary_df["cluster1 name"] = summary_df.apply(
        lambda row: cluster_name_dict[row["cluster1"]]
        if row["cluster1"] in cluster_name_dict.keys()
        else "Nan",
        axis=1,
    )
    summary_df["cluster2 name"] = summary_df.apply(
        lambda row: cluster_name_dict[row["cluster2"]]
        if row["cluster2"] in cluster_name_dict.keys()
        else "Nan",
        axis=1,
    )
    summary_df.set_index(["cluster1", "cluster2"], inplace=True)

    norm_shared_columns = [c for c in summary_df.columns if "norm_shared" in c]
    summary_df["max_norm_shared"] = summary_df.apply(
        lambda row: max([row[c] for c in norm_shared_columns]), axis=1
    )
    summary_df["mean_norm_shared"] = summary_df.apply(
        lambda row: statistics.mean([row[c] for c in norm_shared_columns]), axis=1
    )

    cluster1_multimaper_fraction_columns = [
        c for c in summary_df.columns if "cluster1 multimapper fraction" in c
    ]
    summary_df["max_ref1_multimaper_fraction"] = summary_df.apply(
        lambda row: max([row[c] for c in cluster1_multimaper_fraction_columns]), axis=1
    )
    summary_df["mean_ref1_multimaper_fraction"] = summary_df.apply(
        lambda row: statistics.mean(
            [row[c] for c in cluster1_multimaper_fraction_columns]
        ),
        axis=1,
    )
    summary_df["min_ref1_multimaper_fraction"] = summary_df.apply(
        lambda row: min([row[c] for c in cluster1_multimaper_fraction_columns]), axis=1
    )
    summary_df.to_csv(summary_file, sep="\t", index=True)


def get_multimappers_between_clusters(sam_file, cluster_file, cluster_name_file, out_file):
    # Author: Maria Waldl • code@waldl.org
    # Version: 2024-01-24
    with open(cluster_file) as file:
        cluster_dict = yaml.safe_load(file)
    with open(cluster_name_file) as file:
        cluster_name_dict = yaml.safe_load(file)

    df = pd.read_csv(sam_file, sep="\t")
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
    multimapper_df['cluster1 multimapper fraction'] =  multimapper_df.apply(lambda row: row['shared']/row['reads1'] if row['reads1']>0 else 0 , axis =1)
    multimapper_df['cluster1 name'] = multimapper_df.apply(lambda row: cluster_name_dict[row['cluster1']] if row['cluster1'] in cluster_name_dict.keys() else 'Nan',axis=1)
    multimapper_df['cluster2 name'] = multimapper_df.apply(lambda row: cluster_name_dict[row['cluster2']] if row['cluster2'] in cluster_name_dict.keys() else 'Nan',axis=1)

    multimapper_df.to_csv(out_file ,sep = '\t',index=False)

def get_multimapper_stats(sam_file,cluster_file, cluster_name_file):
    with open(cluster_file) as file:
        cluster_dict = yaml.safe_load(file)
    with open(cluster_name_file) as file:
        cluster_name_dict = yaml.safe_load(file)
    df = pd.read_csv(sam_file, sep="\t")
    read_count = len(list(set(df['QNAME'].to_list())))

    df['CLUSTER'] = df.apply(lambda row: cluster_dict[row['RNAME']] if row['RNAME'] in cluster_dict.keys() else 'Nan' , axis = 1)
    df.drop_duplicates(subset=["CLUSTER", "QNAME"], keep='first', inplace=True)
    df=df.groupby('QNAME').count()
    counts = df['CLUSTER'].value_counts()
    print(counts)
    uniquemappers = counts[1]
    multimappers = read_count-uniquemappers
    doublemappers = counts[2]
    tripplemappers = counts[3]
    highermappers = read_count - counts[1] - counts[2] - counts[3]

    return read_count, uniquemappers, multimappers, doublemappers, tripplemappers, highermappers


rule get_multimapper_stats_r:
    input:
        sams=expand('resources/filtered-mappings/pre-filter_{reads_filter}/{ref_set}/mapped/{sample}.sam', reads_filter='{reads_filter}', ref_set = '{ref_set}', sample=SAMPLES, dist='{dist}'),
        clusters = 'resources/cluster/pre-filter_{reads_filter}/{ref_set}/clusters-ed-{e_cutoff}-mm-{m_cutoff}_{c_treatment}.yaml',
        cluster_names = 'resources/cluster/pre-filter_{reads_filter}/{ref_set}/anticodonbased_clusternames-ed-{e_cutoff}-mm-{m_cutoff}_{c_treatment}.yaml',
    output:
        mutimapper_summary = 'resources/multimappers/pre-filter_{reads_filter}/{ref_set}/multimapper_stats-{e_cutoff}-mm-{m_cutoff}_{c_treatment}_multimappers_between_clusters.tsv'
    run:
        import pandas as pd
        data = []
        for sam in input.sams:
            sample = sam.split('/')[-1].replace('.sam','')
            print(sample)

            treatment = sample_dict[sample]['treatment']
            time = sample_dict[sample]['timepoint']
            mapped_reads, unique_reads, multimappers, double, tripple, higher =  get_multimapper_stats(sam, input.clusters, input.cluster_names)
            data.append({'sample': sample,
                         'treatment': treatment,
                         'time point': time,
                         'mapped reads': mapped_reads,
                         'unique reads': unique_reads,
                         'multimappers': multimappers,
                         'double mappers': double,
                         'tripple mappers': tripple,
                         '>3 clusters mapper': higher,
                         'unique fraction': unique_reads/mapped_reads,
                         'multimapper fraction': multimappers/mapped_reads,
                         })
        df = pd.DataFrame.from_dict(data)
        df.set_index('sample', inplace =True)
        df.to_csv(output.mutimapper_summary,sep='\t')
        print(df.head())


rule get_multimappers_per_sample:
    input:
        sam = 'resources/filtered-mappings/pre-filter_{reads_filter}/{ref_set}/mapped/{sample}.sam'
    output:
        pair_tsv = 'resources/multimappers/pre-filter_{reads_filter}/{ref_set}/per-sample/pair_{sample}.tsv'
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
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

rule get_multimappers_summary_per_ref:
    input:
        expand('resources/multimappers/pre-filter_{reads_filter}/{ref_set}/per-sample/pair_{sample}.tsv', reads_filter='{reads_filter}', ref_set = '{ref_set}', sample=dm_samples)
    output:
        summary_tsv = 'resources/multimappers/pre-filter_{reads_filter}/{ref_set}/multimapper_between_refs.tsv',
        pivot_tsv = 'resources/multimappers/pre-filter_{reads_filter}/{ref_set}/multimapper_between_refs_pivot.tsv',
        cluster_map = 'resources/cluster/pre-filter_{reads_filter}/{ref_set}/clusters_for_sorting_between_refs_multimapper_plot.yaml',
        pdf = 'results/multimappers/pre-filter_{reads_filter}/{ref_set}/multimapper_between_refs.pdf'
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        import pandas as pd
        import yaml
        import numpy as np
        import seaborn as sns
        import statistics

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

        plot_multimappers_between_refs(summary_df, output.pdf, output.pivot_tsv, output.cluster_map)





rule get_multimappers_between_editdist_cluster_per_sample:
    input:
        sam = 'resources/filtered-mappings/pre-filter_{reads_filter}/{ref_set}/mapped/{sample}.sam',
        clusters =  'resources/cluster/pre-filter_{reads_filter}/{ref_set}/clusters-editdist-{dist}.yaml',
        cluster_names = 'resources/cluster/pre-filter_{reads_filter}/{ref_set}/clusternames-editdist-{dist}.yaml',
    output:
        cluster_tsv = 'resources/multimappers/pre-filter_{reads_filter}/{ref_set}/per-sample-editdist-{dist}/pair_{sample}.tsv'
    run:
        get_multimappers_between_clusters(input.sam, input.clusters, input.cluster_names, output.cluster_tsv)


rule summarize_multimappers_between_editdist_clusters:
    input:
        tsvs=expand('resources/multimappers/pre-filter_{reads_filter}/{ref_set}/per-sample-editdist-{dist}/pair_{sample}.tsv', reads_filter='{reads_filter}', ref_set = '{ref_set}', sample=SAMPLES, dist='{dist}'),
        cluster_names = 'resources/cluster/pre-filter_{reads_filter}/{ref_set}/clusternames-editdist-{dist}.yaml',
    output:
        tsv = 'resources/multimappers/pre-filter_{reads_filter}/{ref_set}/cluster-editdist-{dist}_multimappers_between_clusters_{treatment}.tsv'
    run:
        tsv_files = []
        for file in input.tsvs:
            sample = file.split('/')[-1].replace('.tsv','').replace('pair_', '')
            if sample_dict[sample]['treatment']==wildcards.treatment:
                tsv_files.append(file)
        summarize_multimappers_between_clusters(tsv_files, input.cluster_names, output.tsv)


rule plot_multimappers_between_editdist_clusters:
    input:
        tsv = 'resources/multimappers/pre-filter_{reads_filter}/{ref_set}/cluster-editdist-{dist}_multimappers_between_clusters_{treatment}.tsv'
    output:
        pdf = 'results/multimappers/pre-filter_{reads_filter}/{ref_set}/cluster-editdist-{dist}_multimappers_between_clusters_{treatment}.pdf',
    run:
        import pandas as pd
        df = pd.read_csv(
                    input.tsv,
                    sep="\t",
                    )
        title = f'multimappers between edit distance based clusters (max edit dist.: {wildcards.dist}; treatment: {wildcards.treatment})'
        plot_multimappers_between_clusters(df, output.pdf, title)


rule get_multimappers_between_multimapper_cluster_per_sample:
    input:
        sam = 'resources/filtered-mappings/pre-filter_{reads_filter}/{ref_set}/mapped/{sample}.sam',
        clusters = 'resources/cluster/pre-filter_{reads_filter}/{ref_set}/clusters-ed-{e_cutoff}-mm-{m_cutoff}_{c_treatment}.yaml',
        cluster_names = 'resources/cluster/pre-filter_{reads_filter}/{ref_set}/anticodonbased_clusternames-ed-{e_cutoff}-mm-{m_cutoff}_{c_treatment}.yaml',
    output:
        cluster_tsv = 'resources/multimappers/pre-filter_{reads_filter}/{ref_set}/per-sample-ed-{e_cutoff}-mm-{m_cutoff}_{c_treatment}/pair_{sample}.tsv',
    run:
        get_multimappers_between_clusters(input.sam, input.clusters, input.cluster_names, output.cluster_tsv)




rule summarize_multimappers_between_multimapper_clusters:
    input:
        tsvs = expand('resources/multimappers/pre-filter_{reads_filter}/{ref_set}/per-sample-ed-{e_cutoff}-mm-{m_cutoff}_{c_treatment}/pair_{sample}.tsv', reads_filter='{reads_filter}', ref_set = '{ref_set}', sample=SAMPLES, e_cutoff='{e_cutoff}', m_cutoff='{m_cutoff}', c_treatment = '{c_treatment}'),
        cluster_names = 'resources/cluster/pre-filter_{reads_filter}/{ref_set}/anticodonbased_clusternames-ed-{e_cutoff}-mm-{m_cutoff}_{c_treatment}.yaml',
    output:
        tsv = 'resources/multimappers/pre-filter_{reads_filter}/{ref_set}/cluster-editdist-ed-{e_cutoff}-mm-{m_cutoff}_{c_treatment}_multimappers_between_clusters_{treatment}.tsv'
    run:
        tsv_files = []
        for file in input.tsvs:
            sample = file.split('/')[-1].replace('.tsv','').replace('pair_', '')
            if sample_dict[sample]['treatment']==wildcards.treatment:
                tsv_files.append(file)
        summarize_multimappers_between_clusters(tsv_files, input.cluster_names, output.tsv)


rule plot_multimappers_between_multimapper_clusters:
    input:
        tsv = 'resources/multimappers/pre-filter_{reads_filter}/{ref_set}/cluster-editdist-ed-{e_cutoff}-mm-{m_cutoff}_{c_treatment}_multimappers_between_clusters_{treatment}.tsv'
    output:
        pdf = 'results/multimappers/pre-filter_{reads_filter}/{ref_set}/cluster-editdist-ed-{e_cutoff}-mm-{m_cutoff}_{c_treatment}_multimappers_between_clusters_{treatment}.pdf'
    run:
        import pandas as pd
        df = pd.read_csv(
                    input.tsv,
                    sep="\t",
                    )
        title = f'multimappers between edit distance plus multimapper based clusters (max edit dist.: {wildcards.e_cutoff}; maximum multimapper percentage to merge: {wildcards.m_cutoff}; samples for clustering {wildcards.c_treatment}; treatment of shown samples: {wildcards.treatment})'
        plot_multimappers_between_clusters(df, output.pdf, title)
