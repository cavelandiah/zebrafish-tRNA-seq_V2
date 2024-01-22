


rule duplication_rate:
    input:
        unfiltered = 'resources/coverage/pre-filter_none/{ref_set}/coverage_summary_all.tsv',
        filtered = 'resources/coverage/pre-filter_umi/{ref_set}/coverage_summary_all.tsv',
    output:
        tsv = 'qc/umi/{ref_set}_duplication_rate.tsv',
        pdf = 'qc/umi/{ref_set}_duplication_rate.pdf',
    run:
        import pandas as pd
        import seaborn as sn
        import matplotlib.pyplot as plt

        def rename(col):

            sample = col.split(' ')[0]
            time = sample_dict[sample]['timepoint']
            treatment = sample_dict[sample]['treatment']
            name = f'{treatment} {time} {sample}'
            return name

        def get_coverage_count_df(path):
            df = pd.read_csv(
            path,
            sep="\t",
            index_col = 'RNAME'
            )
            df.drop(columns = [c for c in df.columns if 'random count' not in c], inplace = True)
            column_dict = {}
            for c in df.columns:
                column_dict[c]=rename(c)
            df.rename(columns=column_dict, inplace = True)
            df.loc['sum'] = df.sum()
            return df

        df_none = get_coverage_count_df(input.unfiltered)
        df_umi = get_coverage_count_df(input.filtered)
        print(df_none.head())
        print(df_umi.head())
        df = df_umi.div(df_none)
        df = df.reindex(sorted(df.columns), axis=1)
        df['mean'] = df.mean(axis=1)
        df['min'] = df.min(axis=1)
        df['max'] = df.max(axis=1)

        for treatment in ['DM', 'MOCK', 'BS']:
            treatment_cols = [c for c in df.columns if c.startswith(treatment)]
            df[f'{treatment} mean'] = df[treatment_cols].mean(axis=1)
            df[f'{treatment} min'] = df[treatment_cols].min(axis=1)
            df[f'{treatment} max'] = df[treatment_cols].max(axis=1)


        fig, axs = plt.subplots(figsize=(45*CM, 120*CM), nrows= 1, ncols=1)

        sn.heatmap(df, ax=axs, vmin =0, vmax=1)
        fig.savefig(output.pdf)

        df.to_csv(output.tsv,
                 sep = '\t',
                 index=True)

        print(df.head())
