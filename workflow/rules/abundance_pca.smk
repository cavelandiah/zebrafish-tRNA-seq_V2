# Author: Maria Waldl • code@waldl.org
# Version: 2024-01-24

import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sn
import timeimport pandas as pd


def plot_explained_variance(data, path):

    pca = PCA().fit(data)
    var_exp = pca.explained_variance_ratio_
    cum_var_exp = np.cumsum(var_exp)

    fig, axs = plt.subplots(figsize=(16*CM,5*CM), nrows=1, ncols=2, constrained_layout=True)

    x_values = list(range(1,len(var_exp)+1))

    # plot explained variance for all PCs
    axs[0].bar(x_values,var_exp)
    axs[0].plot(x_values,cum_var_exp, c = 'tab:gray', label = 'cumulative')

    # plot explained variance for the first 'pc_number_to_show' PCss
    x_name = [f'PC{i}' for i in x_values]
    pc_number_to_show = 10
    if len(x_name) < pc_number_to_show:
        pc_number_to_show = len(x_name)
    axs[1].bar(x_name[0:pc_number_to_show ],var_exp[0:pc_number_to_show ])
    axs[1].plot(x_name[0:pc_number_to_show ],cum_var_exp[0:pc_number_to_show ], c = 'tab:gray', label = 'cumulative')
    axs[1].set_xticks(list(range(0,pc_number_to_show )))
    axs[1].set_xticklabels(x_name[0:pc_number_to_show ],  rotation= 90)

    # some tick and label formatting
    for ax in axs.flatten():
        ax.set_ylim(0,1)
        ax.set_yticks([0,0.2, 0.4, 0.6, 0.8, 1])
        ax.set_yticklabels([0,0.2, 0.4, 0.6, 0.8, 1], )
        ax.set_ylabel('explained variance')
        ax.set_xlabel('principle components')
        ax.legend(fontsize=8)

    # print information on how many PCs are needed to cover 80, 95
    # and 99 pecent of the total variance
    cum_var_targets= [0.8, 0.95, 0.99]
    target_pcs = []
    for target_var in cum_var_targets:
        for i,var in enumerate(cum_var_exp):
            if var >= target_var:
                target_pcs.append([target_var, i+1, var])
                break
    textstr = [f'>={pc_info[0]}: PC {pc_info[1]} ({pc_info[2]:.4f})' for pc_info in target_pcs]
    textstr = '\n'.join(['cumulative expl. variance:']+textstr)
    axs[0].text(0.45, 0.05, textstr, fontsize=7, transform=axs[0].transAxes,)

    fig.savefig(path, bbox_inches="tight")

def plot_feature_contribution_to_components(data, features, path):
    pca = PCA().fit(data)
    comp_df = pd.DataFrame(pca.components_, columns = features)
    fig, axs = plt.subplots(figsize=( 1+0.5*(len(features))*CM,0.5*len(comp_df)*CM), nrows=1, ncols=1)

    sn.heatmap(comp_df, ax=axs, center =0, square=False)

    plt.tick_params(axis='both', which='major',labelbottom = False, bottom=False, top = True, labeltop=True)
    x_ticks = [n+0.5 for n in range(0,len(features))]
    axs.set_xticks(x_ticks)
    axs.set_xticklabels([l for l in features], rotation=90)

    y_ticks = [n+0.5 for n in range(0,len(comp_df))]
    axs.set_yticks(y_ticks)
    axs.set_yticklabels([f'PC{i}' for i in range(1,len(comp_df)+1)])
    fig.savefig(path, bbox_inches="tight")

def plot_main_pc_contributions(data, features, path):
    pca = PCA().fit(data)

    fig, axs = plt.subplots(figsize=(13*CM,10*CM*len(pca.components_)), nrows=len(pca.components_), ncols=1, )
    plt.subplots_adjust(hspace=1.5)
    for i,comp in enumerate(pca.components_):
        plot_df = pd.DataFrame(comp, columns = ['comp'] ,index = features, )
        plot_df = plot_df[(plot_df['comp'] <= -0.09) | (plot_df['comp']>=0.09)]
        if len(plot_df) <=0:
            continue
        plot_df.sort_values('comp', ascending=False ,inplace=True)

        plot_df.plot.bar(ax= axs[i])

    fig.savefig(path)


def plot_pc_2D(data,index_names, path, text = True ,offset = 0):
    pc_data = PCA(n_components=2).fit_transform(data)
    pca = PCA().fit(data)
    var_exp = pca.explained_variance_ratio_

    plot_df = pd.DataFrame(data = pc_data, columns = ['principal component 1', 'principal component 2'], index = index_names.to_list())
    plot_df.reset_index(inplace = True)

    plot_df['timepoint'] = plot_df.apply(lambda row: int(row['index'].split(' ')[1]), axis =1)
    plot_df['treatment'] = plot_df.apply(lambda row: row['index'].split(' ')[0], axis =1)


    fig, axs = plt.subplots(figsize=(13*CM,13*CM), nrows=1, ncols=1)
    var_exp = pca.explained_variance_ratio_
    color_dict = {1: 'tab:blue',
                    2: 'tab:orange',
                    3: 'tab:green',
                    4: 'tab:red',
                    5: 'tab:purple',
                    6: 'tab:brown',
                    7:'tab:pink',}



    sn.scatterplot(data = plot_df,
                    x = 'principal component 1',
                    y = 'principal component 2',
                    style = 'treatment',
                    hue = 'timepoint',
                    palette = color_dict,
                    ax = axs,
                    #alpha = 0.8,
                    )
    var_comp1 = var_exp[0]
    var_comp2 = var_exp[1]
    axs.set_xlabel(f'principal component 1 ({var_comp1*100:.1f}%)')
    axs.set_ylabel(f'principal component 2 ({var_comp2*100:.1f}%)')

    if text:
        x_coords = plot_df['principal component 1']
        y_coords = plot_df['principal component 2']
        names = index_names
        for i,name in enumerate(names):
            axs.text(x_coords[i]+offset, y_coords[i]-offset, name, fontsize = 5, c = 'tab:gray')

    fig.savefig(path, bbox_inches="tight")





rule get_abundance_pca:
    input:
        tsv = 'resources/coverage/pre-filter_{reads_filter}/{ref_set}/coverage_summary_all.tsv',
    output:
        pdf = 'results/abundance/pre-filter_{reads_filter}/{ref_set}/pca-qc/abundance_pca_{treatment}.pdf',
        variance_pdf = 'results/abundance/pre-filter_{reads_filter}/{ref_set}/pca-qc/explained_variance_{treatment}.pdf',
        pdf_raw = 'results/abundance/pre-filter_{reads_filter}/{ref_set}/pca-qc/abundance_pca_raw_{treatment}.pdf',
        variance_pdf_raw = 'results/abundance/pre-filter_{reads_filter}/{ref_set}/pca-qc/explained_variance_raw_{treatment}.pdf',
    wildcard_constraints:
        treatment="[BSalDMOCK]+"
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        import pandas as pd
        from sklearn.decomposition import PCA
        from sklearn.preprocessing import StandardScaler
        import numpy as np
        import matplotlib.pyplot as plt
        import seaborn as sn

        def rename(name):
            sample = name.split(" ")[0]
            treatment = sample_dict[sample]["treatment"]
            timepoint = sample_dict[sample]["timepoint"]
            new_name = f"{treatment} {timepoint} {sample}"
            return new_name


        df = pd.read_csv(input.tsv, sep="\t", index_col="RNAME")

        drop_columns = [c for c in df.columns if config['abundance_score'] != ' '.join(c.split(' ')[1:])]
        df.drop(
            columns=drop_columns,
            inplace=True,
        )


        col_rename_dict = {}
        for c in df.columns:
            col_rename_dict[c]=rename(c)
        df.rename( col_rename_dict, axis="columns", inplace = True)
        df.replace(np.nan, 0, inplace = True)



        if wildcards.treatment != 'all':
            col_names = list(df.columns)
            drop_rows = [row for row in df if wildcards.treatment not in row]
            df = df.drop(columns=drop_rows)



        df = df.transpose()

        # raw data
        plot_explained_variance(df, str(output.variance_pdf_raw))
        plot_pc_2D(df,df.index,output.pdf_raw)
        plot_main_pc_contributions(df, df.columns, output.variance_pdf_raw+'_feature_contributions_bar.pdf')
        plot_feature_contribution_to_components(df, df.columns,  output.variance_pdf_raw+'_feature_contributions_heatmap.pdf')


        # normalized data
        data_st =  StandardScaler().fit_transform(df)
        plot_explained_variance(data_st, output.variance_pdf)
        plot_pc_2D(data_st,df.index,output.pdf)
        plot_main_pc_contributions(data_st, df.columns, output.variance_pdf+'_feature_contributions_bar.pdf')
        plot_feature_contribution_to_components(data_st, df.columns,  output.variance_pdf+'_feature_contributions_heatmap.pdf')



rule get_abundance_pca_per_timepoint_cluster:
    input:
        tsv = 'resources/coverage/pre-filter_{reads_filter}/{ref_set}/clusters-ed-{e_cutoff}-mm-{m_cutoff}_{c_treatment}_coverage_summary_all.tsv',
    output:
        pdf = 'results/abundance/pre-filter_{reads_filter}/{ref_set}/pca-cluster-time/ed-{e_cutoff}-mm-{m_cutoff}_{c_treatment}/abundance_pca_{treatment}.pdf',
        variance_pdf = 'results/abundance/pre-filter_{reads_filter}/{ref_set}/pca-cluster-time/ed-{e_cutoff}-mm-{m_cutoff}_{c_treatment}/explained_variance_{treatment}.pdf',
        pdf_raw = 'results/abundance/pre-filter_{reads_filter}/{ref_set}/pca-cluster-time/ed-{e_cutoff}-mm-{m_cutoff}_{c_treatment}/abundance_pca_raw_{treatment}.pdf',
        variance_pdf_raw = 'results/abundance/pre-filter_{reads_filter}/{ref_set}/pca-cluster-time/ed-{e_cutoff}-mm-{m_cutoff}_{c_treatment}/explained_variance_raw_{treatment}.pdf',
    wildcard_constraints:
        treatment="[BSalDMOCK]+"
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        import pandas as pd
        from sklearn.decomposition import PCA
        from sklearn.preprocessing import StandardScaler
        import numpy as np
        import matplotlib.pyplot as plt
        import seaborn as sn
        import time

        def rename(name):
            sample = name.split(" ")[0]
            treatment = sample_dict[sample]["treatment"]
            timepoint = sample_dict[sample]["timepoint"]
            new_name = f"{treatment} {timepoint} {sample}"
            return new_name

        df = pd.read_csv(input.tsv, sep="\t")


        drop_columns = [c for c in df.columns if not c.split(' ')[0] in SAMPLES]
        drop_columns.remove('cluster name')
        df.drop(
            columns=drop_columns,
            inplace=True,
        )

        df = df.groupby("cluster name").sum()

        col_rename_dict = {}
        for c in df.columns:
            col_rename_dict[c]=rename(c)

        df.rename( col_rename_dict, axis="columns", inplace = True)
        df.replace(np.nan, 0, inplace = True)



        if wildcards.treatment != 'all':
            col_names = list(df.columns)
            drop_rows = [row for row in df if wildcards.treatment not in row]

            df = df.drop(columns=drop_rows)

        df = df.transpose()
        df.reset_index(inplace =True)
        df['treatment timepoint'] = df.apply(lambda row: ' '.join(row['index'].split(' ')[0:2]), axis = 1)

        df.drop(columns=['index'], inplace=True)
        df =df.groupby('treatment timepoint').mean()



        # raw data
        plot_explained_variance(df, output.variance_pdf_raw)
        plot_pc_2D(df,df.index,output.pdf_raw)
        plot_main_pc_contributions(df, df.columns, output.variance_pdf_raw+'_feature_contributions_bar.pdf')
        plot_feature_contribution_to_components(df, df.columns,  output.variance_pdf_raw+'_feature_contributions_heatmap.pdf')


        # normalized data
        data_st =  StandardScaler().fit_transform(df)
        plot_explained_variance(data_st, output.variance_pdf)
        plot_pc_2D(data_st,df.index,output.pdf)
        plot_main_pc_contributions(data_st, df.columns, output.variance_pdf+'_feature_contributions_bar.pdf')
        plot_feature_contribution_to_components(data_st, df.columns,  output.variance_pdf+'_feature_contributions_heatmap.pdf')


rule get_abundance_pca_per_sample_cluster:
    input:
        tsv = 'resources/coverage/pre-filter_{reads_filter}/{ref_set}/clusters-ed-{e_cutoff}-mm-{m_cutoff}_{c_treatment}_coverage_summary_all.tsv',
    output:
        pdf = 'results/abundance/pre-filter_{reads_filter}/{ref_set}/pca-cluster-sample/ed-{e_cutoff}-mm-{m_cutoff}_{c_treatment}/abundance_pca_{treatment}.pdf',
        variance_pdf = 'results/abundance/pre-filter_{reads_filter}/{ref_set}/pca-cluster-sample/ed-{e_cutoff}-mm-{m_cutoff}_{c_treatment}/explained_variance_{treatment}.pdf',
        pdf_raw = 'results/abundance/pre-filter_{reads_filter}/{ref_set}/pca-cluster-sample/ed-{e_cutoff}-mm-{m_cutoff}_{c_treatment}/abundance_pca_raw_{treatment}.pdf',
        variance_pdf_raw = 'results/abundance/pre-filter_{reads_filter}/{ref_set}/pca-cluster-sample/ed-{e_cutoff}-mm-{m_cutoff}_{c_treatment}/explained_variance_raw_{treatment}.pdf',
    wildcard_constraints:
        treatment="[BSalDMOCK]+"
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        import pandas as pd
        from sklearn.decomposition import PCA
        from sklearn.preprocessing import StandardScaler
        import numpy as np
        import matplotlib.pyplot as plt
        import seaborn as sn
        import time

        def rename(name):
            sample = name.split(" ")[0]
            treatment = sample_dict[sample]["treatment"]
            timepoint = sample_dict[sample]["timepoint"]
            new_name = f"{treatment} {timepoint} {sample}"
            return new_name



        df = pd.read_csv(input.tsv, sep="\t")


        drop_columns = [c for c in df.columns if not c.split(' ')[0] in SAMPLES]
        drop_columns.remove('cluster name')
        df.drop(
            columns=drop_columns,
            inplace=True,
        )

        df = df.groupby("cluster name").sum()

        col_rename_dict = {}
        for c in df.columns:
            col_rename_dict[c]=rename(c)
        df.rename( col_rename_dict, axis="columns", inplace = True)
        df.replace(np.nan, 0, inplace = True)



        if wildcards.treatment != 'all':
            col_names = list(df.columns)
            drop_rows = [row for row in df if wildcards.treatment not in row]

            df = df.drop(columns=drop_rows)

        df = df.transpose()

        # raw data
        plot_explained_variance(df, str(output.variance_pdf_raw))
        plot_pc_2D(df,df.index,output.pdf_raw)
        plot_main_pc_contributions(df, df.columns, output.variance_pdf_raw+'_feature_contributions_bar.pdf')
        plot_feature_contribution_to_components(df, df.columns,  output.variance_pdf_raw+'_feature_contributions_heatmap.pdf')


        # normalized data
        data_st =  StandardScaler().fit_transform(df)
        plot_explained_variance(data_st, output.variance_pdf)
        plot_pc_2D(data_st,df.index,output.pdf)
        plot_main_pc_contributions(data_st, df.columns, output.variance_pdf+'_feature_contributions_bar.pdf')
        plot_feature_contribution_to_components(data_st, df.columns,  output.variance_pdf+'_feature_contributions_heatmap.pdf')


rule get_all_abundance_pca_plots:
    input:
        'results/abundance/pre-filter_{reads_filter}/{ref_set}/pca-cluster-sample/ed-3-mm-50_DM/abundance_pca_BS.pdf',
        'results/abundance/pre-filter_{reads_filter}/{ref_set}/pca-cluster-sample/ed-3-mm-50_DM/abundance_pca_DM.pdf',
        'results/abundance/pre-filter_{reads_filter}/{ref_set}/pca-cluster-sample/ed-3-mm-50_DM/abundance_pca_MOCK.pdf',
        'results/abundance/pre-filter_{reads_filter}/{ref_set}/pca-cluster-sample/ed-3-mm-50_DM/abundance_pca_all.pdf',
        'results/abundance/pre-filter_{reads_filter}/{ref_set}/pca-cluster-time/ed-3-mm-50_DM/abundance_pca_BS.pdf',
        'results/abundance/pre-filter_{reads_filter}/{ref_set}/pca-cluster-time/ed-3-mm-50_DM/abundance_pca_DM.pdf',
        'results/abundance/pre-filter_{reads_filter}/{ref_set}/pca-cluster-time/ed-3-mm-50_DM/abundance_pca_MOCK.pdf',
        'results/abundance/pre-filter_{reads_filter}/{ref_set}/pca-cluster-time/ed-3-mm-50_DM/abundance_pca_all.pdf',
        'results/abundance/pre-filter_{reads_filter}/{ref_set}/pca-qc/abundance_pca_MOCK.pdf',
        'results/abundance/pre-filter_{reads_filter}/{ref_set}/pca-qc/abundance_pca_DM.pdf',
        'results/abundance/pre-filter_{reads_filter}/{ref_set}/pca-qc/abundance_pca_BS.pdf',
        'results/abundance/pre-filter_{reads_filter}/{ref_set}/pca-qc/abundance_pca_all.pdf',
        'results/abundance/pre-filter_{reads_filter}/{ref_set}/pca-qc/abundance_pca_MOCK.pdf',
        'results/abundance/pre-filter_{reads_filter}/{ref_set}/pca-qc/abundance_pca_DM.pdf',
        'results/abundance/pre-filter_{reads_filter}/{ref_set}/pca-qc/abundance_pca_BS.pdf',
        'results/abundance/pre-filter_{reads_filter}/{ref_set}/pca-qc/abundance_pca_all.pdf',
    output:
        'results/abundance/pcs-done.txt'
    run:
        with open(str(output),'w') as file:
            for f in input:
                file.write(f+'\n')
