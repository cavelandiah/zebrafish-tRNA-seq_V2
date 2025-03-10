# Author: Maria Waldl • code@waldl.org
# Version: 2024-01-24

# TODO:
# Calibrate for automatic labels, is hard-coded to 8 datapoints
# axs.set_xlim(0,8)

rule per_cluster_per_timepoint_position_reference:
    input:
        mismatch_tsvs = expand('resources/coverage_counts/pre-filter_{reads_filter}/{ref_set}/clusters-ed-{e_cutoff}-mm-{m_cutoff}_DM/{sample}.tsv', reads_filter = '{reads_filter}',ref_set='{ref_set}', sample = SAMPLES, e_cutoff='{e_cutoff}', m_cutoff='{m_cutoff}' ),
        clusters = 'resources/cluster/pre-filter_{reads_filter}/{ref_set}/clusters-ed-{e_cutoff}-mm-{m_cutoff}_DM.yaml',
        cluster_names = 'resources/cluster/pre-filter_{reads_filter}/{ref_set}/anticodonbased_clusternames-ed-{e_cutoff}-mm-{m_cutoff}_DM.yaml',
    output:
        merged_pdf = 'results/modifications/pre-filter_{reads_filter}/{ref_set}/clusters-ed-{e_cutoff}-mm-{m_cutoff}_DM/per-cluster-and-timepoint_position-reference-heatmap/{treatment}/{mismatch_type}/all_clusters.pdf'
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        import pandas as pd
        import yaml
        import os
        import matplotlib.pyplot as plt
        import seaborn as sn
        import subprocess
        import numpy as np

        # select column to plot
        plot_column = ''
        if wildcards.mismatch_type == 'mismatch':
            plot_column = 'mismatch fraction'
        elif wildcards.mismatch_type== 'm5C':
            plot_column = 'm5C fraction'
        elif wildcards.mismatch_type== 'stop':
            plot_column = 'read start fraction'
        elif wildcards.mismatch_type== 'stopfraction':
            plot_column = 'stop fraction'
        elif wildcards.mismatch_type== 'mismatchStop':
            plot_column = 'stop+mismatch fraction'
        elif wildcards.mismatch_type == 'positionCoverage':
            plot_column = 'RPM at position'
        else:
            print("can't plot this column")

        # read data files
        data = []
        for tsv_file in input.mismatch_tsvs:
            sample = tsv_file.split('/')[-1].replace('.tsv', '')
            #sample = tsv_file.split('/')[-1].replace('.tsv', '').split('_')[0]
            if sample_dict[sample]['treatment'] != wildcards.treatment:
                continue
            df = pd.read_csv(tsv_file, sep = '\t')
            df['sample'] = sample
            df['timepoint'] = sample_dict[sample]['timepoint']
            df['treatment'] = sample_dict[sample]['treatment']
            data.append(df)
        df =  pd.concat(data)

        # get cluster names
        with open(input.clusters) as file:
            cluster_dict = yaml.safe_load(file)
        with open(input.cluster_names) as file:
            cluster_name_dict = yaml.safe_load(file)

        # generate one plot per cluster
        plot_dir = '/'.join(output.merged_pdf.split('/')[0:-1])
        pdfs = []
        for cluster, cdf in df.groupby('cluster'):
            cluster_name = cluster_name_dict[int(cluster)]

            pdf = os.path.join(plot_dir, cluster_name+'.pdf')
            pdfs.append(pdf)

            # It's only one timepoint
            fig, axs = plt.subplots(figsize=(45,7), nrows=1, ncols=2, gridspec_kw={'width_ratios': [11, 1]})
            #fig, axs = plt.subplots(figsize=(32,45), nrows=7, ncols=2, gridspec_kw={'width_ratios': [11, 1]})
            fig.subplots_adjust(hspace=1)
            plt.rc('axes', titlesize=50)
            plt.rc('axes', labelsize=20)
            plt.suptitle(str(cluster)+': '+cluster_name+'\n'+wildcards.treatment+' '+plot_column, size = 40)
            for timepoint, tdf in cdf.groupby('timepoint'):
                tdf['position'] = tdf.apply(lambda row: ' '.join([str(row['align_pos']).zfill(3), str(row['canonical_pos']).zfill(3)]), axis = 1)
                pivot_df = pd.pivot_table(tdf, values=plot_column, index='RNAME', columns='position', aggfunc='mean', fill_value=None, margins=False, dropna=False, margins_name='All', observed=False, sort=True)
                vmax = 1
                if plot_column == 'RPM at position':
                    vmax = None
                sn.heatmap(pivot_df, ax=axs[0], vmin = 0, vmax = vmax, square = False, linewidths= False, cmap = 'Blues')
                count_df = tdf.groupby(['RNAME']).mean(numeric_only=True)
                for c in count_df.columns:
                    if c != 'RPM':
                        count_df.drop(c,axis=1, inplace=True)
                sn.heatmap(count_df, ax=axs[1], vmin = 0, linewidths= False, cmap = 'Blues')
                axs[0].set_title(str(timepoint)+': ' + sample_dict[tdf['sample'].to_list()[0]]['timepoint_name'], loc='left')
            fig.savefig(pdf, bbox_inches="tight")
            fig.clf()
            plt.close(fig)
            plt.close("all")
        pdfs.sort()
        call_args =['pdfunite'] + pdfs + [output.merged_pdf]
        results = subprocess.run(call_args, capture_output=True)


rule plot_per_timepoint_position_clusters_heatmap:
    input:
        mismatch_tsvs = expand('resources/coverage_counts/pre-filter_{reads_filter}/{ref_set}/clusters-ed-{e_cutoff}-mm-{m_cutoff}_DM/{sample}_per_cluster.tsv', reads_filter = '{reads_filter}',ref_set='{ref_set}',sample = SAMPLES, e_cutoff='{e_cutoff}', m_cutoff='{m_cutoff}' ),
        clusters = 'resources/cluster/pre-filter_{reads_filter}/{ref_set}/clusters-ed-{e_cutoff}-mm-{m_cutoff}_DM.yaml',
        cluster_names = 'resources/cluster/pre-filter_{reads_filter}/{ref_set}/anticodonbased_clusternames-ed-{e_cutoff}-mm-{m_cutoff}_DM.yaml'
    output:
        pdf = 'results/modifications/pre-filter_{reads_filter}/{ref_set}/clusters-ed-{e_cutoff}-mm-{m_cutoff}_DM/per-timepoint_position-cluster_heatmap/{treatment}/{mismatch_type}/summary.pdf'
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        import pandas as pd
        import yaml
        import matplotlib.pyplot as plt
        import seaborn as sn
        import numpy as np
        import subprocess

        plot_column = ''
        if wildcards.mismatch_type == 'mismatch':
            plot_column = 'mismatch fraction'
        elif wildcards.mismatch_type== 'm5C':
            plot_column = 'm5C fraction'
        elif wildcards.mismatch_type== 'stop':
            plot_column = 'read start fraction'
        elif wildcards.mismatch_type== 'stopfraction':
            plot_column = 'stop fraction'
        elif wildcards.mismatch_type== 'mismatchStop':
            plot_column = 'stop+mismatch fraction'
        elif wildcards.mismatch_type == 'positionCoverage':
            plot_column = 'RPM at position'
        elif wildcards.mismatch_type == 'refCfraction' :
            plot_column = 'ref C / all nts'
        else:
            print("can't pllt this column")


        data = []
        for tsv_file in input.mismatch_tsvs:
            sample = tsv_file.split('/')[-1].replace('_per_cluster.tsv', '')
            #sample = tsv_file.split('/')[-1].replace('.tsv', '').split('_')[0]
            if sample_dict[sample]['treatment'] != wildcards.treatment:
                continue
            df = pd.read_csv(tsv_file, sep = '\t')
            df['sample'] = sample
            df['timepoint'] = sample_dict[sample]['timepoint']
            df['treatment'] = sample_dict[sample]['treatment']
            #df['representation'] = df['RPM at position']/df['RPM']
            #df = df[df['representation']>= 0.15]
            df['m5C fraction'] = df.apply(lambda row: row['m5C fraction'] if row['ref C / all nts']>=config['min_C_fraction'] else np.nan, axis =1)
            #for cluster in config['clusters_to_remove']:
            #    df = df[df['cluster']!= cluster]
            data.append(df)
        df =  pd.concat(data)

        #df = df[df['cluster']!= "Nan"]

        with open(input.clusters) as file:
            cluster_dict = yaml.safe_load(file)
        with open(input.cluster_names) as file:
            cluster_name_dict = yaml.safe_load(file)

        df['cluster_name'] = df.apply(lambda row: cluster_name_dict[int(row['cluster'])] if int(row['cluster']) in  cluster_name_dict.keys() else 'todo'  , axis =1)

        plot_dir = '/'.join(output.pdf.split('/')[0:-1])
        pdfs=[]
        plt.rc('axes', titlesize=7)
        plt.rc('axes', labelsize=7)
        for time_point, tdf in df.groupby('timepoint'):

            pdf = os.path.join(plot_dir, str(time_point)+'.pdf')
            pdfs.append(pdf)
            fig, axs = plt.subplots(figsize=(32*CM,20*CM), nrows=1, ncols=1)

            tdf['position'] = tdf.apply(lambda row: ' '.join([str(row['align_pos']).zfill(3), str(row['canonical_pos']).zfill(3)]), axis = 1)
            pivot_df = pd.pivot_table(tdf, values=plot_column, index='cluster_name', columns='position', aggfunc='mean', fill_value=None, margins=False, dropna=True, margins_name='All', observed=False, sort=True)
            vmax = 1
            if plot_column =='RPM at position':
                vmax = None
            g=sn.heatmap(pivot_df, ax=axs, vmin = 0, vmax = vmax, square = True, linewidths= False, cmap = 'Blues')
            axs.set_title('timepoint: ' + str(time_point) +'; ' +wildcards.treatment+'; '+plot_column, size = 7)

            yticks = [t +0.5 for t in range(0,len(pivot_df))]
            g.set_yticks(yticks)
            g.set_yticklabels([c.split('[')[0] for c in pivot_df.index.tolist()], fontsize = 7)

            xticks = [t +0.5 for t in range(0,len(pivot_df.columns))]
            g.set_xticks(xticks)
            g.set_xticklabels([c.split(' ')[1].lstrip('0') for c in pivot_df.columns], fontsize = 7)

            fig.savefig(pdf, bbox_inches="tight")
            fig.clf()
            plt.close(fig)
            plt.close("all")

        pdfs.sort()
        call_args =['pdfunite'] + pdfs + [output.pdf]
        results = subprocess.run(call_args, capture_output=True)


rule per_cluster_sample_position_mismatch_heatmap:
    input:
        mismatch_tsvs = expand('resources/coverage_counts/pre-filter_{reads_filter}/{ref_set}/clusters-ed-{e_cutoff}-mm-{m_cutoff}_DM/{sample}_per_cluster.tsv',reads_filter = '{reads_filter}',ref_set='{ref_set}',sample = SAMPLES, e_cutoff='{e_cutoff}', m_cutoff='{m_cutoff}' ),
        clusters = 'resources/cluster/pre-filter_{reads_filter}/{ref_set}/clusters-ed-{e_cutoff}-mm-{m_cutoff}_DM.yaml',
        cluster_names = 'resources/cluster/pre-filter_{reads_filter}/{ref_set}/anticodonbased_clusternames-ed-{e_cutoff}-mm-{m_cutoff}_DM.yaml'
    output:
        pdf = 'results/modifications/pre-filter_{reads_filter}/{ref_set}/clusters-ed-{e_cutoff}-mm-{m_cutoff}_DM/per-cluster_sample-position_heatmap/{treatment}/{mismatch_type}/summary.pdf'
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        import pandas as pd
        import yaml
        import matplotlib.pyplot as plt
        import seaborn as sn
        import os
        import numpy as np
        import subprocess

        plot_column = ''
        if wildcards.mismatch_type == 'mismatch':
            plot_column = 'mismatch fraction'
        elif wildcards.mismatch_type== 'm5C':
            plot_column = 'm5C fraction'
        elif wildcards.mismatch_type== 'stop':
            plot_column = 'read start fraction'
        elif wildcards.mismatch_type== 'stopfraction':
            plot_column = 'stop fraction'
        elif wildcards.mismatch_type== 'mismatchStop':
            plot_column = 'stop+mismatch fraction'
        elif wildcards.mismatch_type == 'positionCoverage':
            plot_column = 'RPM at position'
        elif wildcards.mismatch_type == 'refCfraction':
            plot_column = 'ref C / all nts'
        else:
            print("can't plot this column")

        print(plot_column)

        data = []
        for tsv_file in input.mismatch_tsvs:
            sample = tsv_file.split('/')[-1].replace('_per_cluster.tsv', '')
            #sample = tsv_file.split('/')[-1].replace('.tsv', '').split('_')[0]
            if sample_dict[sample]['treatment'] != wildcards.treatment:
                continue
            df = pd.read_csv(tsv_file, sep = '\t')
            df['sample'] = sample
            df['timepoint'] = sample_dict[sample]['timepoint']
            df['time point name']  = sample_dict[sample]['timepoint_name']
            df['treatment'] = sample_dict[sample]['treatment']
            #df['representation'] = df['RPM at position']/df['RPM']
            #df = df[df['representation']>= 0.15]
            df['m5C fraction'] = df.apply(lambda row: row['m5C fraction'] if row['ref C / all nts']>=0.15 else np.nan, axis =1)
            data.append(df)
        df =  pd.concat(data)

        #df = df[df['cluster']!= "Nan"]

        with open(input.clusters) as file:
            cluster_dict = yaml.safe_load(file)
        with open(input.cluster_names) as file:
            cluster_name_dict = yaml.safe_load(file)

        df['cluster_name'] = df.apply(lambda row: cluster_name_dict[int(row['cluster'])] if int(row['cluster']) in  cluster_name_dict.keys() else 'todo'  , axis =1)

        plot_dir = '/'.join(output.pdf.split('/')[0:-1])
        pdfs=[]
        plt.rc('axes', titlesize=7)
        plt.rc('axes', labelsize=7)
        for cluster, cdf in df.groupby('cluster_name'):

            pdf = os.path.join(plot_dir, cluster+'.pdf')
            pdfs.append(pdf)
            fig, axs = plt.subplots(figsize=(7*CM,20*CM), nrows=1, ncols=1)
            #fig.subplots_adjust(hspace=0.2)


            cdf['position'] = cdf.apply(lambda row: ' '.join([str(row['align_pos']).zfill(3), str(row['canonical_pos']).zfill(3)]), axis = 1)
            cdf['time'] = cdf.apply(lambda row: ' '.join([str(row['timepoint']), row['time point name'], row['sample']]), axis=1)

            pivot_df = pd.pivot_table(cdf, values=plot_column, index='position', columns='time', aggfunc='mean', fill_value=None, margins=False, dropna=False, margins_name='All', observed=False, sort=True)

            vmax = 1
            if plot_column =='RPM at position':
                vmax = None
            plt.rc('axes', titlesize=7)
            plt.rc('axes', labelsize=7)
            g = sn.heatmap(pivot_df, ax=axs, vmin = 0, vmax = vmax, square = False, linewidths= 0.000001, cmap = 'Blues')
            #g.set_facecolor('xkcd:salmon')

            yticks = [t +0.5 for t in range(0,len(pivot_df))]
            g.set_yticks(yticks)
            g.set_yticklabels(pivot_df.index.tolist(), fontsize = 7)

            xticks = [t +0.5 for t in range(0,len(pivot_df.columns))]
            g.set_xticks(xticks)
            g.set_xticklabels(pivot_df.columns, fontsize = 7)

            axs.set_title(cluster+'\n'+wildcards.treatment+' '+plot_column, size = 7)

            fig.savefig(pdf, bbox_inches="tight")
            fig.clf()
            plt.close(fig)
            plt.close("all")

        pdfs.sort()
        call_args =['pdfunite'] + pdfs + [output.pdf]
        results = subprocess.run(call_args, capture_output=True)

rule per_cluster_timepoint_position_mismatch_heatmap:
    input:
        mismatch_tsvs = expand('resources/coverage_counts/pre-filter_{reads_filter}/{ref_set}/clusters-ed-{e_cutoff}-mm-{m_cutoff}_DM/{sample}_per_cluster.tsv', reads_filter = '{reads_filter}',ref_set='{ref_set}',sample = SAMPLES, e_cutoff='{e_cutoff}', m_cutoff='{m_cutoff}' ),
        clusters = 'resources/cluster/pre-filter_{reads_filter}/{ref_set}/clusters-ed-{e_cutoff}-mm-{m_cutoff}_DM.yaml',
        cluster_names = 'resources/cluster/pre-filter_{reads_filter}/{ref_set}/anticodonbased_clusternames-ed-{e_cutoff}-mm-{m_cutoff}_DM.yaml'
    output:
        pdf = 'results/modifications/pre-filter_{reads_filter}/{ref_set}/clusters-ed-{e_cutoff}-mm-{m_cutoff}_DM/per-cluster_timepoint-position_heatmap/{treatment}/{mismatch_type}/summary.pdf'
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        import pandas as pd
        import yaml
        import matplotlib.pyplot as plt
        import seaborn as sn
        import os
        import numpy as np
        import subprocess

        plot_column = ''
        if wildcards.mismatch_type == 'mismatch':
            plot_column = 'mismatch fraction'
        elif wildcards.mismatch_type== 'm5C':
            plot_column = 'm5C fraction'
        elif wildcards.mismatch_type== 'stop':
            plot_column = 'read start fraction'
        elif wildcards.mismatch_type== 'stopfraction':
            plot_column = 'stop fraction'
        elif wildcards.mismatch_type== 'mismatchStop':
            plot_column = 'stop+mismatch fraction'
        elif wildcards.mismatch_type == 'positionCoverage':
            plot_column = 'RPM at position'
        elif wildcards.mismatch_type == 'refCfraction':
            plot_column = 'ref C / all nts'
        else:
            print("can't plot this column")


        data = []
        for tsv_file in input.mismatch_tsvs:
            sample = tsv_file.split('/')[-1].replace('_per_cluster.tsv', '')
            #sample = tsv_file.split('/')[-1].replace('.tsv', '').split('_')[0]
            if sample_dict[sample]['treatment'] != wildcards.treatment:
                continue
            df = pd.read_csv(tsv_file, sep = '\t')
            df['sample'] = sample
            df['timepoint'] = sample_dict[sample]['timepoint']
            df['time point name']  = sample_dict[sample]['timepoint_name']
            df['treatment'] = sample_dict[sample]['treatment']
            #df['representation'] = df['RPM at position']/df['RPM']
            #df = df[df['representation']>= 0.15]
            df['m5C fraction'] = df.apply(lambda row: row['m5C fraction'] if row['ref C / all nts']>=config['min_C_fraction'] else np.nan, axis =1)
            #df['m5C fraction'] = df.apply(lambda row: row['m5C fraction'] if row['ref C / all nts']*row['representation']>=0.15 else np.nan, axis =1)
            data.append(df)
        df =  pd.concat(data)

        #df = df[df['cluster']!= "Nan"]

        with open(input.clusters) as file:
            cluster_dict = yaml.safe_load(file)
        with open(input.cluster_names) as file:
            cluster_name_dict = yaml.safe_load(file)

        df['cluster_name'] = df.apply(lambda row: cluster_name_dict[int(row['cluster'])] if int(row['cluster']) in  cluster_name_dict.keys() else 'todo'  , axis =1)

        plot_dir = '/'.join(output.pdf.split('/')[0:-1])

        pdfs=[]
        plt.rc('axes', titlesize=7)
        plt.rc('axes', labelsize=7)
        for cluster, cdf in df.groupby('cluster_name'):

            pdf = os.path.join(plot_dir, cluster+'.pdf')
            pdfs.append(pdf)
            fig, axs = plt.subplots(figsize=(3*CM,20*CM), nrows=1, ncols=1)
            fig.subplots_adjust(hspace=0.2)


            cdf['position'] = cdf.apply(lambda row: ' '.join([str(row['align_pos']).zfill(3), str(row['canonical_pos']).zfill(3)]), axis = 1)
            cdf['time'] = cdf.apply(lambda row: ' '.join([str(row['timepoint']), row['time point name']]), axis=1)

            pivot_df = pd.pivot_table(cdf, values=plot_column, index='position', columns='time', aggfunc='mean', fill_value=None, margins=False, dropna=False, margins_name='All', observed=False, sort=True)

            vmax = 1
            if plot_column =='RPM at position':
                vmax = None
            plt.rc('axes', titlesize=7)
            plt.rc('axes', labelsize=7)
            g = sn.heatmap(pivot_df, ax=axs, vmin = 0, vmax = vmax, square = False, linewidths= 0.00001, cmap = 'Blues')
            #g.set_facecolor('xkcd:salmon')
            #g.set_facecolor('darkgrey')


            yticks = [t +0.5 for t in range(0,len(pivot_df))]
            g.set_yticks(yticks)
            g.set_yticklabels([c.split(' ')[1].lstrip('0') for c in pivot_df.index.tolist()], fontsize = 7)

            xticks = [t +0.5 for t in range(0,len(pivot_df.columns))]
            g.set_xticks(xticks)
            g.set_xticklabels(pivot_df.columns, fontsize = 7)

            axs.set_title(cluster+'\n'+wildcards.treatment+' '+plot_column, size = 7)

            fig.savefig(pdf, bbox_inches="tight")
            fig.clf()
            plt.close(fig)
            plt.close("all")

        pdfs.sort()
        call_args =['pdfunite'] + pdfs + [output.pdf]
        results = subprocess.run(call_args, capture_output=True)


rule per_pos_cluster_mismatch_line_plot:
    input:
        mismatch_tsvs = expand('resources/coverage_counts/pre-filter_{reads_filter}/{ref_set}/clusters-ed-{e_cutoff}-mm-{m_cutoff}_DM/{sample}_per_cluster.tsv',reads_filter = '{reads_filter}',ref_set='{ref_set}', sample = SAMPLES,  e_cutoff='{e_cutoff}', m_cutoff='{m_cutoff}' ),
        clusters = 'resources/cluster/pre-filter_{reads_filter}/{ref_set}/clusters-ed-{e_cutoff}-mm-{m_cutoff}_DM.yaml',
        cluster_names = 'resources/cluster/pre-filter_{reads_filter}/{ref_set}/anticodonbased_clusternames-ed-{e_cutoff}-mm-{m_cutoff}_DM.yaml'
    output:
        pdf = 'results/modifications/pre-filter_{reads_filter}/{ref_set}/clusters-ed-{e_cutoff}-mm-{m_cutoff}_DM/per-cluster_per-position-mismatch-fraction/{treatment}/min_cov_{min_cov}/{mismatch_type}/summary.pdf'
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        import pandas as pd
        import yaml
        import matplotlib.pyplot as plt
        import matplotlib
        import seaborn as sn
        import os
        import subprocess

        plot_column = ''
        if wildcards.mismatch_type == 'mismatch':
            plot_column = 'mismatch fraction'
        elif wildcards.mismatch_type== 'm5C':
            plot_column = 'm5C fraction'
        elif wildcards.mismatch_type== 'stop':
            plot_column = 'read start fraction'
        elif wildcards.mismatch_type== 'stopfraction':
            plot_column = 'stop fraction'
        elif wildcards.mismatch_type== 'mismatchStop':
            plot_column = 'stop+mismatch fraction'
        elif wildcards.mismatch_type == 'positionCoverage':
            plot_column = 'RPM at position'
        elif wildcards.mismatch_type == 'refCfraction':
            plot_column = 'ref C / all nts'
        else:
            print("can't pllt this column")


        data = []
        for tsv_file in input.mismatch_tsvs:
            sample = tsv_file.split('/')[-1].replace('_per_cluster.tsv', '')
            #sample = tsv_file.split('/')[-1].replace('.tsv', '').split('_')[0]
            if sample_dict[sample]['treatment'] != wildcards.treatment:
                continue
            df = pd.read_csv(tsv_file, sep = '\t')
            df['sample'] = sample
            df['timepoint'] = sample_dict[sample]['timepoint']
            df['time point name']  = sample_dict[sample]['timepoint_name']
            df['treatment'] = sample_dict[sample]['treatment']
            data.append(df)
        df =  pd.concat(data)

        df = df[df['cluster']!= "Nan"]

        with open(input.clusters) as file:
            cluster_dict = yaml.safe_load(file)
        with open(input.cluster_names) as file:
            cluster_name_dict = yaml.safe_load(file)


        df['cluster_name'] = df.apply(lambda row: cluster_name_dict[int(row['cluster'])] if int(row['cluster']) in  cluster_name_dict.keys() else 'todo'  , axis =1)

        plot_dir = '/'.join(output.pdf.split('/')[0:-1])

        pdfs=[]

        colors = [list(plt.cm.tab10(index)) for index in range(0,10)]
        colors = colors *8

        plot_symbols = config['replicate_markes']
        df['replicate_symbol'] = df.apply(lambda row: plot_symbols[int(sample_dict[row['sample']]['replicate'])-1], axis = 1)

        plt.rc('font', size=7) #controls default text size
        plt.rc('axes', titlesize=7) #fontsize of the title
        plt.rc('axes', labelsize=7) #fontsize of the x and y labels
        plt.rc('xtick', labelsize=7) #fontsize of the x tick labels
        plt.rc('ytick', labelsize=7) #fontsize of the y tick labels
        plt.rc('legend', fontsize=7) #fontsize of the legend

        for cluster, cdf in df.groupby('cluster_name'):

            pdf = os.path.join(plot_dir, cluster+'.pdf')
            pdfs.append(pdf)
            fig, axs = plt.subplots(figsize=(8.5*CM,3.3*CM), nrows=1, ncols=1, layout = 'tight')

            plt.rc('font', size=7) #controls default text size
            plt.rc('axes', titlesize=7) #fontsize of the title
            plt.rc('axes', labelsize=7) #fontsize of the x and y labels
            plt.rc('xtick', labelsize=7) #fontsize of the x tick labels
            plt.rc('ytick', labelsize=7) #fontsize of the y tick labels
            plt.rc('legend', fontsize=7) #fontsize of the legend

            cdf['position'] = cdf.apply(lambda row: str(row['canonical_pos']).zfill(2) + ' ('+str(row['align_pos']).zfill(3)+ ')', axis = 1)
            cdf['time'] = cdf.apply(lambda row: ' '.join([str(row['timepoint']), row['time point name']]), axis=1)

            i = 0
            for pos, p_df in  cdf.groupby('position'):
                cutoff = config['modification_lineplots_cutoff']['other'] #mismatch and others
                if wildcards.mismatch_type == 'stop':
                    cutoff = config['modification_lineplots_cutoff']['stop']
                elif wildcards.mismatch_type == 'm5C':
                    cutoff = config['modification_lineplots_cutoff']['m5C']
                if p_df[plot_column].fillna(value = 0).max() <cutoff and wildcards.mismatch_type != 'positionCoverage':
                    continue
                if p_df['RPM at position'].max() < int(wildcards.min_cov):
                    continue

                p_df['rpm factor'] = p_df['RPM']/p_df['all reads']
                p_df['ref_C_rpm'] = p_df.apply(lambda row: row['ref_C']*row['rpm factor'], axis =1)

                mean_p_df = p_df.groupby('timepoint').mean(numeric_only=True)
                if mean_p_df[plot_column].fillna(value = 0).max() <cutoff and wildcards.mismatch_type != 'positionCoverage':
                    continue
                if wildcards.mismatch_type == 'm5C' and mean_p_df['ref_C_rpm'].fillna(value = 0).max() < config['min_C_fraction']:
                    continue

                mean_p_df.reset_index(inplace = True)
                mean_p_df['time'] = mean_p_df.apply(lambda row: sample_name_dict[row['timepoint']]['timepoint_name'].replace(',', ',\n'), axis =1)
                mean_p_df.plot.line(x = 'timepoint', y = plot_column, ax = axs, color = colors[i], label = pos, linewidth = 0.75)

                for marker, rep_df in p_df.groupby('replicate_symbol'):
                    rep_df.plot.scatter(x = 'timepoint', y = plot_column, ax = axs, color = colors[i], label = pos, alpha = 0.6,  s=7,  marker =marker ,linewidths=0.6 ) #edgecolors='black',  facecolors= None,
                i +=1

            #axs.set_title(cluster+ ' '+wildcards.treatment+' '+plot_column)
            axs.set_title(cluster+ ' '+wildcards.treatment)
            axs.set_xlabel('')

            if plot_column == 'RPM at position':
                axs[1].set_ylim(bottom = 0)
            else:
                axs.set_ylim(0,1)
                axs.set_yticks([0,0.25,0.5,0.75,1.0])
                axs.set_yticklabels( ["{:.2f}".format(l) for l in  axs.get_yticks()])

            axs.set_xlim(0,8)
            axs.set_xticks(list(range(1,8)))
            axs.set_xticklabels([ sample_name_dict[l]['timepoint_name'].replace(', ', ',\n')  for l in list(range(1,8))], rotation=90)

            handles, labels = axs.get_legend_handles_labels()
            new_labels= []
            new_handles = []
            for i,h in enumerate(handles):
                if isinstance(h, matplotlib.lines.Line2D):
                    new_labels.append(labels[i])
                    new_handles.append(h)
            axs.legend(new_handles, new_labels, bbox_to_anchor=(1.,1), loc="upper left")

            fig.savefig(pdf)
            fig.clf()
            plt.close(fig)
            plt.close("all")

        pdfs.sort()
        call_args =['pdfunite'] + pdfs + [output.pdf]
        results = subprocess.run(call_args, capture_output=True)

rule per_pos_cluster_mismatch_line_plot_small_legend:
    input:
        mismatch_tsvs = expand('resources/coverage_counts/pre-filter_{reads_filter}/{ref_set}/clusters-ed-{e_cutoff}-mm-{m_cutoff}_DM/{sample}_per_cluster.tsv', reads_filter = '{reads_filter}',ref_set='{ref_set}',sample = SAMPLES,  e_cutoff='{e_cutoff}', m_cutoff='{m_cutoff}' ),
        clusters = 'resources/cluster/pre-filter_{reads_filter}/{ref_set}/clusters-ed-{e_cutoff}-mm-{m_cutoff}_DM.yaml',
        cluster_names = 'resources/cluster/pre-filter_{reads_filter}/{ref_set}/anticodonbased_clusternames-ed-{e_cutoff}-mm-{m_cutoff}_DM.yaml'
    output:
        pdf = 'results/modifications/pre-filter_{reads_filter}/{ref_set}/clusters-ed-{e_cutoff}-mm-{m_cutoff}_DM/per-cluster_per-position-mismatch-fraction-smaller-legend/{treatment}/min_cov_{min_cov}/{mismatch_type}/summary.pdf'
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        import pandas as pd
        import yaml
        import matplotlib.pyplot as plt
        import matplotlib
        import seaborn as sn
        import os
        import subprocess

        plot_column = ''
        if wildcards.mismatch_type == 'mismatch':
            plot_column = 'mismatch fraction'
        elif wildcards.mismatch_type== 'm5C':
            plot_column = 'm5C fraction'
        elif wildcards.mismatch_type== 'stop':
            plot_column = 'read start fraction'
        elif wildcards.mismatch_type== 'stopfraction':
            plot_column = 'stop fraction'
        elif wildcards.mismatch_type== 'mismatchStop':
            plot_column = 'stop+mismatch fraction'
        elif wildcards.mismatch_type == 'positionCoverage':
            plot_column = 'RPM at position'
        elif wildcards.mismatch_type == 'refCfraction':
            plot_column = 'ref C / all nts'
        else:
            print("can't pllt this column")


        data = []
        for tsv_file in input.mismatch_tsvs:
            sample = tsv_file.split('/')[-1].replace('_per_cluster.tsv', '')
            #sample = tsv_file.split('/')[-1].replace('.tsv', '').split('_')[0]
            if sample_dict[sample]['treatment'] != wildcards.treatment:
                continue
            df = pd.read_csv(tsv_file, sep = '\t')
            df['sample'] = sample
            df['timepoint'] = sample_dict[sample]['timepoint']
            df['time point name']  = sample_dict[sample]['timepoint_name']
            df['treatment'] = sample_dict[sample]['treatment']
            data.append(df)
        df =  pd.concat(data)

        df = df[df['cluster']!= "Nan"]

        with open(input.clusters) as file:
            cluster_dict = yaml.safe_load(file)
        with open(input.cluster_names) as file:
            cluster_name_dict = yaml.safe_load(file)


        df['cluster_name'] = df.apply(lambda row: cluster_name_dict[int(row['cluster'])] if int(row['cluster']) in  cluster_name_dict.keys() else 'todo'  , axis =1)

        plot_dir = '/'.join(output.pdf.split('/')[0:-1])

        pdfs=[]

        colors = [list(plt.cm.tab10(index)) for index in range(0,10)]
        colors = colors *8

        plot_symbols = config['replicate_markes']
        df['replicate_symbol'] = df.apply(lambda row: plot_symbols[int(sample_dict[row['sample']]['replicate'])-1], axis = 1)

        plt.rc('font', size=7) #controls default text size
        plt.rc('axes', titlesize=7) #fontsize of the title
        plt.rc('axes', labelsize=7) #fontsize of the x and y labels
        plt.rc('xtick', labelsize=7) #fontsize of the x tick labels
        plt.rc('ytick', labelsize=7) #fontsize of the y tick labels
        plt.rc('legend', fontsize=5) #fontsize of the legend

        for cluster, cdf in df.groupby('cluster_name'):

            pdf = os.path.join(plot_dir, cluster+'.pdf')
            pdfs.append(pdf)
            fig, axs = plt.subplots(figsize=(8.5*CM,3.3*CM), nrows=1, ncols=1, layout = 'tight')

            plt.rc('font', size=7) #controls default text size
            plt.rc('axes', titlesize=7) #fontsize of the title
            plt.rc('axes', labelsize=7) #fontsize of the x and y labels
            plt.rc('xtick', labelsize=7) #fontsize of the x tick labels
            plt.rc('ytick', labelsize=7) #fontsize of the y tick labels
            plt.rc('legend', fontsize=5) #fontsize of the legend

            cdf['position'] = cdf.apply(lambda row: str(row['canonical_pos']).zfill(2) + ' ('+str(row['align_pos']).zfill(3)+ ')', axis = 1)
            cdf['time'] = cdf.apply(lambda row: ' '.join([str(row['timepoint']), row['time point name']]), axis=1)

            i = 0
            for pos, p_df in  cdf.groupby('position'):
                cutoff = config['modification_lineplots_cutoff']['other'] #mismatch and others
                if wildcards.mismatch_type == 'stop':
                    cutoff = config['modification_lineplots_cutoff']['stop']
                elif wildcards.mismatch_type == 'm5C':
                    cutoff = config['modification_lineplots_cutoff']['m5C']
                if p_df[plot_column].fillna(value = 0).max() <cutoff and wildcards.mismatch_type != 'positionCoverage':
                    continue
                if p_df['RPM at position'].max() < int(wildcards.min_cov):
                    continue

                p_df['rpm factor'] = p_df['RPM']/p_df['all reads']
                p_df['ref_C_rpm'] = p_df.apply(lambda row: row['ref_C']*row['rpm factor'], axis =1)

                mean_p_df = p_df.groupby('timepoint').mean(numeric_only=True)
                if mean_p_df[plot_column].fillna(value = 0).max() <cutoff and wildcards.mismatch_type != 'positionCoverage':
                    continue
                if wildcards.mismatch_type == 'm5C' and mean_p_df['ref_C_rpm'].fillna(value = 0).max() < config['min_C_fraction']:
                    continue

                mean_p_df.reset_index(inplace = True)
                mean_p_df['time'] = mean_p_df.apply(lambda row: sample_name_dict[row['timepoint']]['timepoint_name'].replace(',', ',\n'), axis =1)
                mean_p_df.plot.line(x = 'timepoint', y = plot_column, ax = axs, color = colors[i], label = pos, linewidth = 0.75)

                for marker, rep_df in p_df.groupby('replicate_symbol'):
                    rep_df.plot.scatter(x = 'timepoint', y = plot_column, ax = axs, color = colors[i], label = pos, alpha = 0.6,  s=7,  marker =marker ,linewidths=0.6 ) #edgecolors='black',  facecolors= None,
                i +=1

            #axs.set_title(cluster+ ' '+wildcards.treatment+' '+plot_column)
            axs.set_title(cluster+ ' '+wildcards.treatment)
            axs.set_xlabel('')

            if plot_column == 'RPM at position':
                axs[1].set_ylim(bottom = 0)
            else:
                axs.set_ylim(0,1)
                axs.set_yticks([0,0.25,0.5,0.75,1.0])
                axs.set_yticklabels( ["{:.2f}".format(l) for l in  axs.get_yticks()])

            axs.set_xlim(0,8)
            axs.set_xticks(list(range(1,8)))
            axs.set_xticklabels([ sample_name_dict[l]['timepoint_name'].replace(', ', ',\n')  for l in list(range(1,8))], rotation=90)

            handles, labels = axs.get_legend_handles_labels()
            new_labels= []
            new_handles = []
            for i,h in enumerate(handles):
                if isinstance(h, matplotlib.lines.Line2D):
                    new_labels.append(labels[i])
                    new_handles.append(h)
            axs.legend(new_handles, new_labels, bbox_to_anchor=(1.,1), loc="upper left")

            fig.savefig(pdf)
            fig.clf()
            plt.close(fig)
            plt.close("all")

        pdfs.sort()
        call_args =['pdfunite'] + pdfs + [output.pdf]
        results = subprocess.run(call_args, capture_output=True)

rule per_pos_cluster_mismatch_coverage_plot:
    input:
        mismatch_tsvs = expand('resources/coverage_counts/pre-filter_{reads_filter}/{ref_set}/clusters-ed-{e_cutoff}-mm-{m_cutoff}_DM/{sample}_per_cluster.tsv',reads_filter = '{reads_filter}',ref_set='{ref_set}', sample = SAMPLES, e_cutoff='{e_cutoff}', m_cutoff='{m_cutoff}' ),
        clusters = 'resources/cluster/pre-filter_{reads_filter}/{ref_set}/clusters-ed-{e_cutoff}-mm-{m_cutoff}_DM.yaml',
        cluster_names = 'resources/cluster/pre-filter_{reads_filter}/{ref_set}/anticodonbased_clusternames-ed-{e_cutoff}-mm-{m_cutoff}_DM.yaml'
    output:
        pdf = 'results/modifications/pre-filter_{reads_filter}/{ref_set}/clusters-ed-{e_cutoff}-mm-{m_cutoff}_DM/per-cluster_per-position-mismatch-fraction-and-coverage/{treatment}/min_cov_{min_cov}/{mismatch_type}/summary.pdf'
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        import pandas as pd
        import yaml
        import matplotlib.pyplot as plt
        import matplotlib
        import seaborn as sn
        import os
        import subprocess

        plot_column = ''
        if wildcards.mismatch_type == 'mismatch':
            plot_column = 'mismatch fraction'
        elif wildcards.mismatch_type== 'm5C':
            plot_column = 'm5C fraction'
        elif wildcards.mismatch_type== 'stop':
            plot_column = 'read start fraction'
        elif wildcards.mismatch_type== 'stopfraction':
            plot_column = 'stop fraction'
        elif wildcards.mismatch_type== 'mismatchStop':
            plot_column = 'stop+mismatch fraction'
        elif wildcards.mismatch_type == 'refCfraction':
            plot_column = 'ref C / all nts'
        else:
            print("can't plot this column")


        data = []
        for tsv_file in input.mismatch_tsvs:
            sample = tsv_file.split('/')[-1].replace('_per_cluster.tsv', '')
            #sample = tsv_file.split('/')[-1].replace('.tsv', '').split('_')[0]
            if sample_dict[sample]['treatment'] != wildcards.treatment:
                continue
            df = pd.read_csv(tsv_file, sep = '\t')
            df['sample'] = sample
            df['timepoint'] = sample_dict[sample]['timepoint']
            df['time point name']  = sample_dict[sample]['timepoint_name']
            df['treatment'] = sample_dict[sample]['treatment']
            #df = df[df['RPM']>= int(wildcards.min_cov)]
            data.append(df)
        df =  pd.concat(data)

        df = df[df['cluster']!= "Nan"]

        with open(input.clusters) as file:
            cluster_dict = yaml.safe_load(file)
        with open(input.cluster_names) as file:
            cluster_name_dict = yaml.safe_load(file)


        df['cluster_name'] = df.apply(lambda row: cluster_name_dict[int(row['cluster'])] if int(row['cluster']) in  cluster_name_dict.keys() else 'todo'  , axis =1)

        plot_dir = '/'.join(output.pdf.split('/')[0:-1])

        pdfs=[]

        colors = [list(plt.cm.tab10(index)) for index in range(0,10)]
        colors = colors *8

        plot_symbols = config['replicate_markes']
        df['replicate_symbol'] = df.apply(lambda row: plot_symbols[int(sample_dict[row['sample']]['replicate'])-1], axis = 1)

        plt.rc('font', size=7) #controls default text size
        plt.rc('axes', titlesize=7) #fontsize of the title
        plt.rc('axes', labelsize=7) #fontsize of the x and y labels
        plt.rc('xtick', labelsize=7) #fontsize of the x tick labels
        plt.rc('ytick', labelsize=7) #fontsize of the y tick labels
        plt.rc('legend', fontsize=7) #fontsize of the legend

        for cluster, cdf in df.groupby('cluster_name'):
            print(cluster)
            print(cdf['RPM at position'].max())
            if cdf['RPM at position'].max() < int(wildcards.min_cov):
                continue
            pdf = os.path.join(plot_dir, cluster+'.pdf')
            pdfs.append(pdf)
            fig, axs = plt.subplots(figsize=(8.5*CM,6*CM), nrows=2, ncols=1, layout = 'tight')


            plt.rc('font', size=7) #controls default text size
            plt.rc('axes', titlesize=7) #fontsize of the title
            plt.rc('axes', labelsize=7) #fontsize of the x and y labels
            plt.rc('xtick', labelsize=7) #fontsize of the x tick labels
            plt.rc('ytick', labelsize=7) #fontsize of the y tick labels
            plt.rc('legend', fontsize=7) #fontsize of the legend

            cdf['position'] = cdf.apply(lambda row: str(row['canonical_pos']).zfill(3) + ' ('+str(row['align_pos'])+ ')', axis = 1)
            cdf['time'] = cdf.apply(lambda row: ' '.join([str(row['timepoint']), row['time point name']]), axis=1)

            i = 0
            for pos, p_df in  cdf.groupby('position'):
                cutoff = config['modification_lineplots_cutoff']['other'] #mismatch and others
                if wildcards.mismatch_type == 'stop':
                    cutoff = config['modification_lineplots_cutoff']['stop']
                elif wildcards.mismatch_type == 'm5C':
                    cutoff = config['modification_lineplots_cutoff']['m5C']
                if p_df[plot_column].fillna(value = 0).max() <cutoff and wildcards.mismatch_type != 'positionCoverage':
                    continue
                if p_df['RPM at position'].max() < int(wildcards.min_cov):
                    continue


                p_df['C and T on refC'] = p_df['C on refC'] + p_df['T on refC']
                p_df['rpm factor'] = p_df['RPM']/p_df['all reads']
                p_df['ref_C_rpm'] = p_df.apply(lambda row: row['ref_C']*row['rpm factor'], axis =1)
                p_df['C and T on refC rpm'] = p_df.apply(lambda row: row['C and T on refC']*row['rpm factor'], axis =1)
                p_df['C on refC rpm'] = p_df.apply(lambda row: row['C on refC']*row['rpm factor'], axis =1)


                mean_p_df = p_df.groupby('timepoint').mean(numeric_only=True)
                if mean_p_df[plot_column].fillna(value = 0).max() <cutoff and wildcards.mismatch_type != 'positionCoverage':
                    continue
                if wildcards.mismatch_type == 'm5C' and mean_p_df['ref_C_rpm'].fillna(value = 0).max() < config['min_C_fraction']:
                    continue


                mean_p_df.reset_index(inplace = True)
                #mean_p_df['timepoint'] = mean_p_df.apply(lambda row: sample_name_dict[row['timepoint']]['timepoint_name'].replace(',', ',\n'), axis =1)
                mean_p_df.plot.line(x = 'timepoint', y = plot_column, ax = axs[0], color = colors[i], label = pos, linewidth = 0.75)
                #axs[0].get_legend().remove()

                for marker, rep_df in p_df.groupby('replicate_symbol'):
                    rep_df.plot.scatter(x = 'timepoint', y = plot_column, ax = axs[0], color = colors[i], label = pos, alpha = 0.6,  s=7,  marker =marker ,linewidths=0.6 ) #edgecolors='black',  facecolors= None,

                if wildcards.mismatch_type== 'm5C':
                    mean_p_df.plot.line(x = 'timepoint', y = 'RPM at position', ax = axs[1], color = colors[i], label = pos, linewidth = 0.75, linestyle = 'dashed', alpha = 0.2)
                    mean_p_df.plot.line(x = 'timepoint', y = 'ref_C_rpm' , ax = axs[1], color = colors[i], label = pos, linewidth = 0.75, linestyle = 'dotted', alpha=0.2)
                    mean_p_df.plot.line(x = 'timepoint', y = 'C and T on refC rpm' , ax = axs[1], color = colors[i], label = pos, linewidth = 0.75, linestyle = 'dashdot', alpha=0.5)
                    mean_p_df.plot.line(x = 'timepoint', y = 'C on refC rpm' , ax = axs[1], color = colors[i], label = pos, linewidth = 0.75, alpha = 0.7)
                else:
                    mean_p_df.plot.line(x = 'timepoint', y = 'RPM at position', ax = axs[1], color = colors[i], label = pos, linewidth = 0.75, linestyle = 'dashed')
                    mean_p_df['mismatch count'] = mean_p_df.apply(lambda row: row[plot_column]*row['RPM at position'], axis =1)
                    mean_p_df.plot.line(x = 'timepoint', y = 'mismatch count' , ax = axs[1], color = colors[i], label = pos, linewidth = 0.75)
                axs[1].get_legend().remove()
                i +=1

            #axs[0].set_xlim(-1,7)
            #axs[0].set_xticks(list(range(0,7)))
            #axs[0].set_xticklabels(axs[1].get_xticklabels(), rotation=90)
            axs[0].set_xlabel('')
            axs[0].set_ylim(0,1)
            axs[0].set_yticks([0.0,0.5,1.0])

            axs[0].set_yticklabels( ["{:.2f}".format(l) for l in  axs[0].get_yticks()])
            axs[0].set_title(cluster+ ' '+wildcards.treatment+' '+plot_column)
            axs[0].set_xlim(0,8)
            axs[0].set_xticks(list(range(1,8)))

            axs[1].set_ylim(bottom = 0)
            axs[1].set_xlabel('')
            axs[1].set_xlim(0,8)
            axs[1].set_xticks(list(range(1,8)))
            axs[1].set_xticklabels([ sample_name_dict[l]['timepoint_name'].replace(', ', ',\n')  for l in list(range(1,8))], rotation=90)

            handles, labels = axs[0].get_legend_handles_labels()
            new_labels= []
            new_handles = []
            for i,h in enumerate(handles):
                if isinstance(h, matplotlib.lines.Line2D):
                    new_labels.append(labels[i])
                    new_handles.append(h)
            axs[0].legend(new_handles, new_labels, bbox_to_anchor=(1.,1), loc="upper left", fontsize =5, markerscale=0.3)


            fig.savefig(pdf)
            fig.clf()
            plt.close(fig)
            plt.close("all")

        pdfs.sort()
        call_args =['pdfunite'] + pdfs + [output.pdf]
        results = subprocess.run(call_args, capture_output=True)



rule per_pos_ref_mismatch_line_plot:
    input:
        mismatch_tsvs = expand('resources/coverage_counts/pre-filter_{reads_filter}/{ref_set}/clusters-ed-{e_cutoff}-mm-{m_cutoff}_DM/{sample}.tsv', reads_filter = '{reads_filter}',ref_set = '{ref_set}', sample = SAMPLES, e_cutoff='{e_cutoff}', m_cutoff='{m_cutoff}' ),
    output:
        pdf = 'results/modifications/pre-filter_{reads_filter}/{ref_set}/clusters-ed-{e_cutoff}-mm-{m_cutoff}_DM/per-high-coverage-ref_modified-position_line-plot/{treatment}/min_cov_{min_cov}/{mismatch_type}/summary.pdf'
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        import pandas as pd
        import yaml
        import matplotlib.pyplot as plt
        import seaborn as sn
        import os
        import subprocess

        plot_column = ''
        if wildcards.mismatch_type == 'mismatch':
            plot_column = 'mismatch fraction'
        elif wildcards.mismatch_type== 'm5C':
            plot_column = 'm5C fraction'
        elif wildcards.mismatch_type== 'stop':
            plot_column = 'read start fraction'
        elif wildcards.mismatch_type== 'stopfraction':
            plot_column = 'stop fraction'
        elif wildcards.mismatch_type== 'mismatchStop':
            plot_column = 'stop+mismatch fraction'
        elif wildcards.mismatch_type == 'positionCoverage':
            plot_column = 'RPM at position'
        elif wildcards.mismatch_type == 'refCfraction':
            plot_column = 'ref C / all nts'
        else:
            print("can't plolt this column")


        data = []
        for tsv_file in input.mismatch_tsvs:
            sample = tsv_file.split('/')[-1].replace('.tsv', '')
            #sample = tsv_file.split('/')[-1].replace('.tsv', '').split('_')[0]
            if sample_dict[sample]['treatment'] != wildcards.treatment:
                continue
            df = pd.read_csv(tsv_file, sep = '\t')
            df['sample'] = sample
            df['timepoint'] = sample_dict[sample]['timepoint']
            df['time point name']  = sample_dict[sample]['timepoint_name']
            df['treatment'] = sample_dict[sample]['treatment']
            #df = df[df['RPM']>= int(wildcards.min_cov)]
            data.append(df)
        df =  pd.concat(data)


        plot_dir = '/'.join(output.pdf.split('/')[0:-1])

        pdfs=[]

        colors = [list(plt.cm.tab10(index)) for index in range(0,10)]
        colors = colors *10
        symbols = ['x', 'o', '^', '1', 's', '+', 'D', 'H']
        symbols = [[symb]*10 for symb in symbols]
        symbols = [item for sublist in symbols for item in sublist]

        for rname, cdf in df.groupby('RNAME'):
            if cdf['RPM'].max() <  int(wildcards.min_cov):
                continue

            pdf = os.path.join(plot_dir, rname+'.pdf')
            pdfs.append(pdf)
            fig, axs = plt.subplots(figsize=(8*CM,4*CM), nrows=1, ncols=1)
            fig.subplots_adjust(hspace=0.2)

            cdf['position'] = cdf.apply(lambda row: str(row['canonical_pos']).zfill(2) + ' ('+str(row['align_pos'])+ ')', axis = 1)
            cdf['time'] = cdf.apply(lambda row: ' '.join([str(row['timepoint']), row['time point name']]), axis=1)


            i = 0

            for pos, p_df in cdf.groupby('position'):
                if p_df[plot_column].fillna(value = 0).max() < 0.25:
                    continue
                mean_p_df = p_df.groupby('timepoint').mean(numeric_only=True)
                mean_p_df.reset_index(inplace = True)
                mean_p_df.plot.line(x = 'timepoint', y = plot_column, ax = axs, color = colors[i], label = pos)
                p_df.plot.scatter(x = 'timepoint', y = plot_column, ax = axs, color = colors[i], label = pos, alpha = 0.6, marker = symbols[i] )
                i +=1


            axs.set_xlim(0,8)
            axs.set_ylim(0,1)
            axs.set_title(rname+'\n'+wildcards.treatment+' '+plot_column, size = 10)

            plt.rc('axes', titlesize=10)
            plt.rc('axes', labelsize=10)

            fig.savefig(pdf, bbox_inches="tight")
            fig.clf()
            plt.close(fig)
            plt.close("all")

        pdfs.sort()
        call_args =['pdfunite'] + pdfs + [output.pdf]
        results = subprocess.run(call_args, capture_output=True)



rule get_all_m5C_plots:
    input:
        #'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-cluster-and-timepoint_position-reference-heatmap/BS/m5C/all_clusters.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-cluster_timepoint-position_heatmap/BS/m5C/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-cluster_timepoint-position_heatmap/BS/refCfraction/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-cluster_sample-position_heatmap/BS/m5C/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-cluster_sample-position_heatmap/BS/refCfraction/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-cluster_per-position-mismatch-fraction-smaller-legend/BS/min_cov_'+str(config['min_coverage'])+'/mismatch/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-cluster_per-position-mismatch-fraction-smaller-legend/BS/min_cov_0/mismatch/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-cluster_per-position-mismatch-fraction/BS/min_cov_'+str(config['min_coverage'])+'/m5C/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-cluster_per-position-mismatch-fraction/BS/min_cov_500/m5C/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-cluster_per-position-mismatch-fraction-and-coverage/BS/min_cov_'+str(config['min_coverage'])+'/m5C/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-cluster_per-position-mismatch-fraction-and-coverage/BS/min_cov_500/m5C/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-cluster_per-position-mismatch-fraction-and-coverage/BS/min_cov_0/m5C/summary.pdf',
        #'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-timepoint_position-cluster_heatmap/BS/m5C/summary.pdf',
        #'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-timepoint_position-cluster_heatmap/BS/refCfraction/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-high-coverage-ref_modified-position_line-plot/BS/min_cov_'+str(config['min_coverage'])+'/m5C/summary.pdf',
        #'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-cluster-and-timepoint_position-reference-heatmap/BS/m5C/all_clusters.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-cluster_timepoint-position_heatmap/BS/m5C/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-cluster_timepoint-position_heatmap/BS/refCfraction/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-cluster_sample-position_heatmap/BS/m5C/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-cluster_sample-position_heatmap/BS/refCfraction/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-cluster_per-position-mismatch-fraction-smaller-legend/BS/min_cov_'+str(config['min_coverage'])+'/mismatch/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-cluster_per-position-mismatch-fraction-smaller-legend/BS/min_cov_0/mismatch/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-cluster_per-position-mismatch-fraction/BS/min_cov_'+str(config['min_coverage'])+'/m5C/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-cluster_per-position-mismatch-fraction/BS/min_cov_500/m5C/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-cluster_per-position-mismatch-fraction-and-coverage/BS/min_cov_'+str(config['min_coverage'])+'/m5C/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-cluster_per-position-mismatch-fraction-and-coverage/BS/min_cov_500/m5C/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-cluster_per-position-mismatch-fraction-and-coverage/BS/min_cov_0/m5C/summary.pdf',
        #'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-timepoint_position-cluster_heatmap/BS/m5C/summary.pdf',
        #'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-timepoint_position-cluster_heatmap/BS/refCfraction/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-high-coverage-ref_modified-position_line-plot/BS/min_cov_'+str(config['min_coverage'])+'/m5C/summary.pdf',
    output:
        'results/modifications/get_all_m5C_calls_plots'
    run:
        with open(str(output),'w') as file:
            for f in input:
                file.write(f+'\n')

rule get_all_mismatch_plots:
    input:
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-cluster-and-timepoint_position-reference-heatmap/{treatment}/mismatch/all_clusters.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-cluster-and-timepoint_position-reference-heatmap/{treatment}/positionCoverage/all_clusters.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-cluster-and-timepoint_position-reference-heatmap/{treatment}/stop/all_clusters.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-cluster-and-timepoint_position-reference-heatmap/{treatment}/stopfraction/all_clusters.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-cluster_timepoint-position_heatmap/{treatment}/positionCoverage/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-cluster_timepoint-position_heatmap/{treatment}/stop/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-cluster_timepoint-position_heatmap/{treatment}/mismatch/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-cluster_timepoint-position_heatmap/{treatment}/stopfraction/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-cluster_timepoint-position_heatmap/{treatment}/mismatchStop/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-cluster_sample-position_heatmap/{treatment}/positionCoverage/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-cluster_sample-position_heatmap/{treatment}/stop/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-cluster_sample-position_heatmap/{treatment}/mismatch/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-cluster_sample-position_heatmap/{treatment}/stopfraction/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-cluster_sample-position_heatmap/{treatment}/mismatchStop/summary.pdf',
        #'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-cluster_per-position-mismatch-fraction-smaller-legend/{treatment}/min_cov_'+str(config['min_coverage'])+'/mismatch/summary.pdf',
        #'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-cluster_per-position-mismatch-fraction-smaller-legend/{treatment}/min_cov_'+str(config['min_coverage'])+'/stop/summary.pdf',
        #'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-cluster_per-position-mismatch-fraction-smaller-legend/{treatment}/min_cov_'+str(config['min_coverage'])+'/mismatchStop/summary.pdf',
        #'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-cluster_per-position-mismatch-fraction-smaller-legend/{treatment}/min_cov_'+str(config['min_coverage'])+'/stopfraction/summary.pdf',
        #'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-cluster_per-position-mismatch-fraction-smaller-legend/{treatment}/min_cov_0/mismatch/summary.pdf',
        #'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-cluster_per-position-mismatch-fraction-smaller-legend/{treatment}/min_cov_0/stop/summary.pdf',
        #'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-cluster_per-position-mismatch-fraction-smaller-legend/{treatment}/min_cov_0/mismatchStop/summary.pdf',
        #'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-cluster_per-position-mismatch-fraction-smaller-legend/{treatment}/min_cov_0/stopfraction/summary.pdf',
        #'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-cluster_per-position-mismatch-fraction/{treatment}/min_cov_'+str(config['min_coverage'])+'/mismatch/summary.pdf',
        #'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-cluster_per-position-mismatch-fraction/{treatment}/min_cov_'+str(config['min_coverage'])+'/stop/summary.pdf',
        #'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-cluster_per-position-mismatch-fraction/{treatment}/min_cov_0/mismatch/summary.pdf',
        #'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-cluster_per-position-mismatch-fraction/{treatment}/min_cov_0/stop/summary.pdf',
        #'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-cluster_per-position-mismatch-fraction-and-coverage/{treatment}/min_cov_'+str(config['min_coverage'])+'/mismatch/summary.pdf',
        #'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-cluster_per-position-mismatch-fraction-and-coverage/{treatment}/min_cov_'+str(config['min_coverage'])+'/stop/summary.pdf',
        #'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-cluster_per-position-mismatch-fraction-and-coverage/{treatment}/min_cov_0/mismatch/summary.pdf',
        #'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-cluster_per-position-mismatch-fraction-and-coverage/{treatment}/min_cov_0/mismatchStop/summary.pdf',
        #'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-cluster_per-position-mismatch-fraction-and-coverage/{treatment}/min_cov_0/stopfraction/summary.pdf',
        #'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-cluster_per-position-mismatch-fraction-and-coverage/{treatment}/min_cov_0/stop/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-timepoint_position-cluster_heatmap/{treatment}/mismatch/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-timepoint_position-cluster_heatmap/{treatment}/stop/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-timepoint_position-cluster_heatmap/{treatment}/stopfraction/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-timepoint_position-cluster_heatmap/{treatment}/mismatchStop/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-timepoint_position-cluster_heatmap/{treatment}/positionCoverage/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-high-coverage-ref_modified-position_line-plot/{treatment}/min_cov_'+str(config['min_coverage'])+'/mismatch/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-high-coverage-ref_modified-position_line-plot/{treatment}/min_cov_'+str(config['min_coverage'])+'/mismatchStop/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-high-coverage-ref_modified-position_line-plot/{treatment}/min_cov_'+str(config['min_coverage'])+'/stop/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-3-mm-50_DM/per-high-coverage-ref_modified-position_line-plot/{treatment}/min_cov_'+str(config['min_coverage'])+'/stopfraction/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-cluster-and-timepoint_position-reference-heatmap/{treatment}/mismatch/all_clusters.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-cluster-and-timepoint_position-reference-heatmap/{treatment}/positionCoverage/all_clusters.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-cluster-and-timepoint_position-reference-heatmap/{treatment}/stop/all_clusters.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-cluster-and-timepoint_position-reference-heatmap/{treatment}/stopfraction/all_clusters.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-cluster_timepoint-position_heatmap/{treatment}/positionCoverage/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-cluster_timepoint-position_heatmap/{treatment}/stop/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-cluster_timepoint-position_heatmap/{treatment}/mismatch/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-cluster_timepoint-position_heatmap/{treatment}/stopfraction/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-cluster_timepoint-position_heatmap/{treatment}/mismatchStop/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-cluster_sample-position_heatmap/{treatment}/positionCoverage/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-cluster_sample-position_heatmap/{treatment}/stop/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-cluster_sample-position_heatmap/{treatment}/mismatch/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-cluster_sample-position_heatmap/{treatment}/stopfraction/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-cluster_sample-position_heatmap/{treatment}/mismatchStop/summary.pdf',
        #'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-cluster_per-position-mismatch-fraction-smaller-legend/{treatment}/min_cov_'+str(config['min_coverage'])+'/mismatch/summary.pdf',
        #'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-cluster_per-position-mismatch-fraction-smaller-legend/{treatment}/min_cov_'+str(config['min_coverage'])+'/stop/summary.pdf',
        #'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-cluster_per-position-mismatch-fraction-smaller-legend/{treatment}/min_cov_'+str(config['min_coverage'])+'/mismatchStop/summary.pdf',
        #'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-cluster_per-position-mismatch-fraction-smaller-legend/{treatment}/min_cov_'+str(config['min_coverage'])+'/stopfraction/summary.pdf',
        #'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-cluster_per-position-mismatch-fraction-smaller-legend/{treatment}/min_cov_0/mismatch/summary.pdf',
        #'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-cluster_per-position-mismatch-fraction-smaller-legend/{treatment}/min_cov_0/stop/summary.pdf',
        #'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-cluster_per-position-mismatch-fraction-smaller-legend/{treatment}/min_cov_0/mismatchStop/summary.pdf',
        #'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-cluster_per-position-mismatch-fraction-smaller-legend/{treatment}/min_cov_0/stopfraction/summary.pdf',
        #'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-cluster_per-position-mismatch-fraction/{treatment}/min_cov_'+str(config['min_coverage'])+'/mismatch/summary.pdf',
        #'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-cluster_per-position-mismatch-fraction/{treatment}/min_cov_'+str(config['min_coverage'])+'/stop/summary.pdf',
        #'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-cluster_per-position-mismatch-fraction/{treatment}/min_cov_0/mismatch/summary.pdf',
        #'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-cluster_per-position-mismatch-fraction/{treatment}/min_cov_0/stop/summary.pdf',
        #'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-cluster_per-position-mismatch-fraction-and-coverage/{treatment}/min_cov_'+str(config['min_coverage'])+'/mismatch/summary.pdf',
        #'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-cluster_per-position-mismatch-fraction-and-coverage/{treatment}/min_cov_'+str(config['min_coverage'])+'/stop/summary.pdf',
        #'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-cluster_per-position-mismatch-fraction-and-coverage/{treatment}/min_cov_0/mismatch/summary.pdf',
        #'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-cluster_per-position-mismatch-fraction-and-coverage/{treatment}/min_cov_0/mismatchStop/summary.pdf',
        #'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-cluster_per-position-mismatch-fraction-and-coverage/{treatment}/min_cov_0/stopfraction/summary.pdf',
        #'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-cluster_per-position-mismatch-fraction-and-coverage/{treatment}/min_cov_0/stop/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-timepoint_position-cluster_heatmap/{treatment}/mismatch/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-timepoint_position-cluster_heatmap/{treatment}/stop/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-timepoint_position-cluster_heatmap/{treatment}/stopfraction/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-timepoint_position-cluster_heatmap/{treatment}/mismatchStop/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-timepoint_position-cluster_heatmap/{treatment}/positionCoverage/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-high-coverage-ref_modified-position_line-plot/{treatment}/min_cov_'+str(config['min_coverage'])+'/mismatch/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-high-coverage-ref_modified-position_line-plot/{treatment}/min_cov_'+str(config['min_coverage'])+'/mismatchStop/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-high-coverage-ref_modified-position_line-plot/{treatment}/min_cov_'+str(config['min_coverage'])+'/stop/summary.pdf',
        'results/modifications/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-2-mm-50_DM/per-high-coverage-ref_modified-position_line-plot/{treatment}/min_cov_'+str(config['min_coverage'])+'/stopfraction/summary.pdf',
    output:
        'results/modifications/get_{treatment}_all_plots'
    run:
        print('test', wildcards.treatment)
        with open(str(output),'w') as file:
            for f in input:
                file.write(f+'\n')
