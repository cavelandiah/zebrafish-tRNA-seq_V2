
rule plot_peranticodon_time_abundance_anticodon_lineperref:
    input:
        tsv = 'resources/coverage/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/coverage_summary_DM.tsv'
    output:
        anticodon_pdf = 'results/abundance/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/anticodon/all_refs/summary.pdf'
    run:
        import pandas as pd
        import matplotlib.pyplot as plt

        gdf = pd.read_csv(input.tsv, sep="\t",)
        gdf.drop(columns = [c for c in gdf.columns if config['abundance_score'] not in c and c not in ['anticodon', 'RNAME']], inplace = True)
        #df.reset_index(inplace = True)
        print(gdf.head())

        pdfs = []
        plot_dir = '/'.join(output.anticodon_pdf.split('/')[0:-1])

        plt.rc('font', size=7) #controls default text size
        plt.rc('axes', titlesize=7) #fontsize of the title
        plt.rc('axes', labelsize=7) #fontsize of the x and y labels
        plt.rc('xtick', labelsize=7) #fontsize of the x tick labels
        plt.rc('ytick', labelsize=7) #fontsize of the y tick labels
        plt.rc('legend', fontsize=7) #fontsize of the legend

        groups  = gdf.groupby('anticodon')
        for anticodon, df in groups:

            fig, axs = plt.subplots(figsize=(8.5*CM, 3.3*CM), nrows= 1, ncols=1, layout = 'tight')
            #fig.subplots_adjust(left=0.5)
            #fig.subplots_adjust(hspace=1)
            #print(plt.rcParams['font.family'])
            #plt.rcParams.update({'font.sans-serif':'Helvetica'})
            #plt.rcParams.update({'font.size':7})

            plt.rc('font', size=7) #controls default text size
            plt.rc('axes', titlesize=7) #fontsize of the title
            plt.rc('axes', labelsize=7) #fontsize of the x and y labels
            plt.rc('xtick', labelsize=7) #fontsize of the x tick labels
            plt.rc('ytick', labelsize=7) #fontsize of the y tick labels
            plt.rc('legend', fontsize=7) #fontsize of the legend


            print(anticodon)
            df.drop(columns=['anticodon'], inplace = True)
            df['RNAME'] = df.apply(lambda row: row['RNAME'].split('(')[0], axis =1)
            #df.reset_index(inplace = True)
            df.set_index('RNAME',inplace=True)
            #df = df*1000000
            df = df.transpose()
            df.reset_index(inplace = True)

            df['raw_fastq'] = df.apply(lambda row: row['index'].split(' ')[0], axis = 1)
            df['time'] = df.apply(lambda row: int(sample_dict[row['raw_fastq']]['timepoint']), axis = 1)
            df['name'] = df.apply(lambda row: sample_dict[row['raw_fastq']]['timepoint_name'], axis = 1)
            #df.drop(columns = ['index', 'time'], inplace = True)
            df.sort_values('time', inplace = True)

            labels = df['name'].to_list()
            labels = sorted(set(labels), key=labels.index)
            group_df = df.groupby('time').mean()
            group_df.reset_index(inplace = True)
            group_df.sort_values('time', inplace = True)
            times = group_df['time'].to_list()


            group_df.set_index('time', inplace = True)
            group_df.plot( ax = axs, c = '0.8',legend = False)


            for col in group_df.columns:
                if group_df[col].max()<0.002:#*1000000:
                    group_df.drop(columns = [col], inplace=True)


            axs.set_title(anticodon)

            if len(group_df.columns) > 0:
                group_df.plot( ax = axs, legend = False)

            axs.set_ylim(0,None)
            #axs.set_yticks(axs.get_yticks())
            axs.set_yticklabels(["{:.3f}".format(l) for l in axs.get_yticks()],)#, max_ylabel_len)
            #axs.tick_params(axis='y', which='major', pad=10)

            axs.set_xticks(times)
            axs.set_xticklabels([l.replace(',', ',\n') for l in labels], rotation = 90)
            axs.set_ylabel('fraction of reads')
            axs.set_xlabel(None)
            handles, labels = axs.get_legend_handles_labels()

            new_labels= []
            new_handles = []
            for i,h in enumerate(handles):
                if h.get_color() != '0.8':
                    new_labels.append(labels[i])
                    new_handles.append(h)
            print(anticodon, len(new_handles))
            if len(new_handles) > 6:
                axs.legend(new_handles, new_labels, bbox_to_anchor=(1.,1), loc="upper left")
            else:
                axs.legend(new_handles, new_labels, bbox_to_anchor=(1.,1), loc="upper left")
            #axs.ticklabel_format(useOffset=False)
            #print(axs.get_legend())
            #axs.legend(loc='upper left')
            #print(axs.get_position())
            #print(axs.get_yaxis())


            fig_path = os.path.join(plot_dir, anticodon + '.pdf')
            pdfs.append(fig_path)

            fig.savefig(fig_path)
            plt.close("all")


        pdfs.sort()
        call_args =['pdfunite'] + pdfs + [output.anticodon_pdf]
        results = subprocess.run(call_args, capture_output=True)


rule plot_percluster_time_abundance_lineperreference:
    input:
        tsv = 'resources/coverage/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/coverage_summary_DM.tsv',
        cluster = 'resources/cluster/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-{e_cutoff}-mm-{m_cutoff}_{treatment}.yaml',
        cluster_name = 'resources/cluster/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/anticodonbased_clusternames-ed-{e_cutoff}-mm-{m_cutoff}_{treatment}.yaml'
    output:
        anticodon_pdf = 'results/abundance/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/cluster_ed-{e_cutoff}-mm-{m_cutoff}/{treatment}_per-cluster_time_abundance_line-per-ref/summary.pdf'
    run:
        import pandas as pd
        import matplotlib.pyplot as plt
        import yaml

        gdf = pd.read_csv(input.tsv, sep="\t",)
        gdf.drop(columns = [c for c in gdf.columns if config['abundance_score'] not in c and c not in [ 'RNAME']], inplace = True)
        #df.reset_index(inplace = True)
        print(gdf.head())

        with open(input.cluster) as file:
            cluster_dict = yaml.safe_load(file)
        with open(input.cluster_name) as file:
            cluster_name_dict = yaml.safe_load(file)
        gdf['cluster'] = gdf.apply(lambda row: cluster_dict[row['RNAME']] if row['RNAME'] in cluster_dict.keys() else 'unknown', axis =1)
        gdf['cluster_name'] = gdf.apply(lambda row: cluster_name_dict[row['cluster']] if row['cluster'] in cluster_name_dict.keys() else 'low coverage', axis =1)

        print(gdf.head())
        gdf.drop(columns = ['cluster'], inplace = True)

        pdfs = []
        plot_dir = '/'.join(output.anticodon_pdf.split('/')[0:-1])

        plt.rc('font', size=7) #controls default text size
        plt.rc('axes', titlesize=7) #fontsize of the title
        plt.rc('axes', labelsize=7) #fontsize of the x and y labels
        plt.rc('xtick', labelsize=7) #fontsize of the x tick labels
        plt.rc('ytick', labelsize=7) #fontsize of the y tick labels
        plt.rc('legend', fontsize=7) #fontsize of the legend

        groups  = gdf.groupby('cluster_name')
        for cluster, df in groups:

            fig, axs = plt.subplots(figsize=(8*CM, 3.3*CM), nrows= 1, ncols=1, layout = 'tight')
            #fig.subplots_adjust(left=0.5)
            #fig.subplots_adjust(hspace=1)
            #print(plt.rcParams['font.family'])
            #plt.rcParams.update({'font.sans-serif':'Helvetica'})
            #plt.rcParams.update({'font.size':7})

            plt.rc('font', size=7) #controls default text size
            plt.rc('axes', titlesize=7) #fontsize of the title
            plt.rc('axes', labelsize=7) #fontsize of the x and y labels
            plt.rc('xtick', labelsize=7) #fontsize of the x tick labels
            plt.rc('ytick', labelsize=7) #fontsize of the y tick labels
            plt.rc('legend', fontsize=7) #fontsize of the legend


            print(cluster)
            df.drop(columns=['cluster_name'], inplace = True)
            df['RNAME'] = df.apply(lambda row: row['RNAME'].split('(')[0], axis =1)
            #df.reset_index(inplace = True)
            df.set_index('RNAME',inplace=True)
            #df = df*1000000
            df = df.transpose()
            df.reset_index(inplace = True)

            df['raw_fastq'] = df.apply(lambda row: row['index'].split(' ')[0], axis = 1)
            df['time'] = df.apply(lambda row: int(sample_dict[row['raw_fastq']]['timepoint']), axis = 1)
            df['name'] = df.apply(lambda row: sample_dict[row['raw_fastq']]['timepoint_name'], axis = 1)
            #df.drop(columns = ['index', 'time'], inplace = True)
            df.sort_values('time', inplace = True)

            labels = df['name'].to_list()
            labels = sorted(set(labels), key=labels.index)
            group_df = df.groupby('time').mean()
            group_df.reset_index(inplace = True)
            group_df.sort_values('time', inplace = True)
            times = group_df['time'].to_list()


            group_df.set_index('time', inplace = True)
            group_df.plot( ax = axs, c = '0.8',legend = False, linewidth=0.75)


            for col in group_df.columns:
                if group_df[col].max()<0.002:#*1000000:
                    group_df.drop(columns = [col], inplace=True)


            axs.set_title(cluster.split('[')[0].replace('_', '/'))

            if len(group_df.columns) > 0:
                group_df.plot( ax = axs, legend = False,  linewidth=0.75)

            axs.set_ylim(0,None)
            axs.set_yticks(axs.get_yticks())
            axs.set_yticklabels(["{:.3f}".format(l) for l in axs.get_yticks()],)#, max_ylabel_len)
            #axs.tick_params(axis='y', which='major', pad=10)

            axs.set_xticks(times)
            axs.set_xticklabels([l.replace(',', ',\n') for l in labels], rotation = 90)
            axs.set_ylabel('fraction of reads')
            axs.set_xlabel(None)
            handles, labels = axs.get_legend_handles_labels()

            new_labels= []
            new_handles = []
            for i,h in enumerate(handles):
                if h.get_color() != '0.8':
                    new_labels.append(labels[i])
                    new_handles.append(h)
            print(cluster, len(new_handles))
            if len(new_handles) > 6:
                axs.legend(new_handles, new_labels, bbox_to_anchor=(1.,1), loc="upper left")
            else:
                axs.legend(new_handles, new_labels, bbox_to_anchor=(1.,1), loc="upper left")
            #axs.ticklabel_format(useOffset=False)
            #print(axs.get_legend())
            #axs.legend(loc='upper left')
            #print(axs.get_position())
            #print(axs.get_yaxis())


            fig_path = os.path.join(plot_dir, cluster + '.pdf')
            pdfs.append(fig_path)

            fig.savefig(fig_path)
            plt.close("all")


        pdfs.sort()
        call_args =['pdfunite'] + pdfs + [output.anticodon_pdf]
        results = subprocess.run(call_args, capture_output=True)

rule plot_percluster_time_abundance_lineperreference_snall_legend:
    input:
        tsv = 'resources/coverage/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/coverage_summary_DM.tsv',
        cluster = 'resources/cluster/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-{e_cutoff}-mm-{m_cutoff}_{treatment}.yaml',
        cluster_name = 'resources/cluster/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/anticodonbased_clusternames-ed-{e_cutoff}-mm-{m_cutoff}_{treatment}.yaml'
    output:
        anticodon_pdf = 'results/abundance/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/cluster_ed-{e_cutoff}-mm-{m_cutoff}/{treatment}_per-cluster_time_abundance_line-per-ref-smaller-legend/summary.pdf'
    run:
        import pandas as pd
        import matplotlib.pyplot as plt
        import yaml

        gdf = pd.read_csv(input.tsv, sep="\t",)
        gdf.drop(columns = [c for c in gdf.columns if config['abundance_score'] not in c and c not in [ 'RNAME']], inplace = True)
        #df.reset_index(inplace = True)
        print(gdf.head())

        with open(input.cluster) as file:
            cluster_dict = yaml.safe_load(file)
        with open(input.cluster_name) as file:
            cluster_name_dict = yaml.safe_load(file)
        gdf['cluster'] = gdf.apply(lambda row: cluster_dict[row['RNAME']] if row['RNAME'] in cluster_dict.keys() else 'unknown', axis =1)
        gdf['cluster_name'] = gdf.apply(lambda row: cluster_name_dict[row['cluster']] if row['cluster'] in cluster_name_dict.keys() else 'low coverage', axis =1)

        print(gdf.head())
        gdf.drop(columns = ['cluster'], inplace = True)

        pdfs = []
        plot_dir = '/'.join(output.anticodon_pdf.split('/')[0:-1])

        plt.rc('font', size=7) #controls default text size
        plt.rc('axes', titlesize=7) #fontsize of the title
        plt.rc('axes', labelsize=7) #fontsize of the x and y labels
        plt.rc('xtick', labelsize=7) #fontsize of the x tick labels
        plt.rc('ytick', labelsize=7) #fontsize of the y tick labels
        plt.rc('legend', fontsize=4) #fontsize of the legend

        groups  = gdf.groupby('cluster_name')
        for cluster, df in groups:

            fig, axs = plt.subplots(figsize=(8*CM, 3.3*CM), nrows= 1, ncols=1, layout = 'tight')
            #fig.subplots_adjust(left=0.5)
            #fig.subplots_adjust(hspace=1)
            #print(plt.rcParams['font.family'])
            #plt.rcParams.update({'font.sans-serif':'Helvetica'})
            #plt.rcParams.update({'font.size':7})

            plt.rc('font', size=7) #controls default text size
            plt.rc('axes', titlesize=7) #fontsize of the title
            plt.rc('axes', labelsize=7) #fontsize of the x and y labels
            plt.rc('xtick', labelsize=7) #fontsize of the x tick labels
            plt.rc('ytick', labelsize=7) #fontsize of the y tick labels
            plt.rc('legend', fontsize=4) #fontsize of the legend


            print(cluster)
            df.drop(columns=['cluster_name'], inplace = True)
            df['RNAME'] = df.apply(lambda row: row['RNAME'].split('(')[0], axis =1)
            #df.reset_index(inplace = True)
            df.set_index('RNAME',inplace=True)
            #df = df*1000000
            df = df.transpose()
            df.reset_index(inplace = True)

            df['raw_fastq'] = df.apply(lambda row: row['index'].split(' ')[0], axis = 1)
            df['time'] = df.apply(lambda row: int(sample_dict[row['raw_fastq']]['timepoint']), axis = 1)
            df['name'] = df.apply(lambda row: sample_dict[row['raw_fastq']]['timepoint_name'], axis = 1)
            #df.drop(columns = ['index', 'time'], inplace = True)
            df.sort_values('time', inplace = True)

            labels = df['name'].to_list()
            labels = sorted(set(labels), key=labels.index)
            group_df = df.groupby('time').mean()
            group_df.reset_index(inplace = True)
            group_df.sort_values('time', inplace = True)
            times = group_df['time'].to_list()


            group_df.set_index('time', inplace = True)
            group_df.plot( ax = axs, c = '0.8',legend = False, linewidth=0.75)


            for col in group_df.columns:
                if group_df[col].max()<0.002:#*1000000:
                    group_df.drop(columns = [col], inplace=True)


            axs.set_title(cluster.split('[')[0].replace('_', '/'))

            if len(group_df.columns) > 0:
                group_df.plot( ax = axs, legend = False,  linewidth=0.75)

            axs.set_ylim(0,None)
            axs.set_yticks(axs.get_yticks())
            axs.set_yticklabels(["{:.3f}".format(l) for l in axs.get_yticks()],)#, max_ylabel_len)
            #axs.tick_params(axis='y', which='major', pad=10)

            axs.set_xticks(times)
            axs.set_xticklabels([l.replace(',', ',\n') for l in labels], rotation = 90)
            axs.set_ylabel('fraction of reads')
            axs.set_xlabel(None)
            handles, labels = axs.get_legend_handles_labels()

            new_labels= []
            new_handles = []
            for i,h in enumerate(handles):
                if h.get_color() != '0.8':
                    new_labels.append(labels[i])
                    new_handles.append(h)
            print(cluster, len(new_handles))
            if len(new_handles) > 6:
                axs.legend(new_handles, new_labels, bbox_to_anchor=(1.,1), loc="upper left")
            else:
                axs.legend(new_handles, new_labels, bbox_to_anchor=(1.,1), loc="upper left")
            #axs.ticklabel_format(useOffset=False)
            #print(axs.get_legend())
            #axs.legend(loc='upper left')
            #print(axs.get_position())
            #print(axs.get_yaxis())


            fig_path = os.path.join(plot_dir, cluster + '.pdf')
            pdfs.append(fig_path)

            fig.savefig(fig_path)
            plt.close("all")


        pdfs.sort()
        call_args =['pdfunite'] + pdfs + [output.anticodon_pdf]
        results = subprocess.run(call_args, capture_output=True)


rule plot_abundance_boxplots_anticodon:
    input:
        tsv = 'resources/coverage/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/coverage_summary_DM.tsv'
    output:
        anticodon_pdf = 'results/abundance/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/anticodon/boxplot.pdf'
    run:
        import pandas as pd
        import matplotlib.pyplot as plt

        df = pd.read_csv(input.tsv, sep="\t",)
        df.drop(columns = [c for c in df.columns if config['abundance_score'] not in c and c != 'anticodon'], inplace = True)
        #df.reset_index(inplace = True)
        df = df.groupby('anticodon').sum()
        df.reset_index(inplace = True)
        df.set_index('anticodon',inplace=True)
        df = df.transpose()
        df.reset_index(inplace = True)
        df['raw_fastq'] = df.apply(lambda row: row['index'].split(' ')[0], axis = 1)
        df['time'] = df.apply(lambda row: int(sample_dict[row['raw_fastq']]['timepoint']), axis = 1)
        df['name'] = df.apply(lambda row: sample_dict[row['raw_fastq']]['timepoint_name'], axis = 1)
        #df.drop(columns = ['index', 'time'], inplace = True)
        df.sort_values('time', inplace = True)
        print(df.head(10))
        labels = df['name'].to_list()
        labels = sorted(set(labels), key=labels.index)
        group_df = df.groupby('time').mean()
        group_df.reset_index(inplace = True)
        group_df.sort_values('time', inplace = True)
        group_df.set_index('time')

        fig, axs = plt.subplots(figsize=(8, 3*len(group_df.columns)), nrows=len(group_df.columns)-1, ncols=1)
        fig.subplots_adjust(hspace=0.5)

        i = 0
        for col in list(group_df.columns):
            if col == 'time':
                #group_df.plot(c='grey', legend = False, ax = axs[i])
                continue
            axs[i].set_xticks(group_df['time'].to_list())
            axs[i].set_xticklabels(labels)


            group_df.plot(x = 'time', y=col,  ax = axs[i])
            #plot_df = df.set_index('name')
            df.boxplot(by = 'time', column = col, ax = axs[i], positions = group_df['time'].to_list(), labels = labels )
            axs[i].set_ylim((0, group_df[col].max()*1.5))
            axs[i].set_ylabel('fraction of reads')
            axs[i].set_xlabel('time points')
            #print(axs[i].get_xticklabels())
            #print(labels)
            axs[i].set_xticklabels(labels+labels)

            i += 1
        fig.savefig(output.anticodon_pdf, bbox_inches="tight")


rule plot_abundance_scatterplots_anticodon:
    input:
        tsv = 'resources/coverage/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/coverage_summary_DM.tsv'
    output:
        anticodon_pdf = 'results/abundance/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/anticodon/scatterplot.pdf'
    run:
        import pandas as pd
        import matplotlib.pyplot as plt

        df = pd.read_csv(input.tsv, sep="\t",)
        df.drop(columns = [c for c in df.columns if config['abundance_score'] not in c and c != 'anticodon'], inplace = True)
        #df.reset_index(inplace = True)
        df = df.groupby('anticodon').sum()
        df.reset_index(inplace = True)
        df.set_index('anticodon',inplace=True)
        df = df.transpose()
        df.reset_index(inplace = True)
        df['raw_fastq'] = df.apply(lambda row: row['index'].split(' ')[0], axis = 1)
        df['time'] = df.apply(lambda row: int(sample_dict[row['raw_fastq']]['timepoint']), axis = 1)
        df['name'] = df.apply(lambda row: sample_dict[row['raw_fastq']]['timepoint_name'], axis = 1)

        #df.drop(columns = ['index', 'time'], inplace = True)
        df.sort_values('time', inplace = True)
        labels = df['name'].to_list()
        labels = sorted(set(labels), key=labels.index)
        print(labels)
        group_df = df.groupby('time').mean()
        group_df.reset_index(inplace = True)
        group_df.sort_values('time', inplace = True)
        group_df.set_index('time')



        fig, axs = plt.subplots(figsize=(8, 3*len(group_df.columns)), nrows=len(group_df.columns)-1, ncols=1)
        fig.subplots_adjust(hspace=0.5)

        plot_symbols = config['replicate_markes']
        df['replicate_symbol'] = df.apply(lambda row: plot_symbols[int(sample_dict[row['raw_fastq']]['replicate'])-1], axis = 1)

        i = 0
        for col in list(group_df.columns):
            if col in ['time', 'replicate_symbol']:
                continue
            axs[i].set_xticks(group_df['time'].to_list())
            axs[i].set_xticklabels(labels)


            group_df.plot(x = 'time', y=col,  ax = axs[i], c = 'grey')
            for marker, rep_df in df.groupby('replicate_symbol'):
                rep_df.plot.scatter(x = 'time', y = col,  ax = axs[i], alpha = 0.7, marker = marker)
            #plot_df = df.set_index('name')
            #df.boxplot(by = 'time', column = col, ax = axs[i], positions = group_df['time'].to_list(), labels = labels )
            axs[i].set_ylim((0, group_df[col].max()*1.5))
            axs[i].set_ylabel('fraction of reads')
            axs[i].set_xlabel('time points')
            #print(axs[i].get_xticklabels())
            #print(labels)
            #axs[i].set_xticklabels(labels+labels)

            i += 1
        fig.savefig(output.anticodon_pdf, bbox_inches="tight")

rule plot_heatmap_anticodon:
    input:
        tsv = 'resources/coverage/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/coverage_summary_DM.tsv'
    output:
        anticodon_pdf = 'results/abundance/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/anticodon/heatmap.pdf'
    run:
        import pandas as pd
        import matplotlib.pyplot as plt
        import seaborn as sn
        df = pd.read_csv(input.tsv, sep="\t",)
        df.drop(columns = [c for c in df.columns if config['abundance_score'] not in c and c != 'anticodon'], inplace = True)
        df = df.groupby('anticodon').sum()

        df = df.transpose()
        df.reset_index(inplace = True)
        df['raw_fastq'] = df.apply(lambda row: row['index'].split(' ')[0], axis = 1)

        df['time'] = df.apply(lambda row: int(sample_dict[row['raw_fastq']]['timepoint']), axis = 1)
        df['name'] = df.apply(lambda row: str(sample_dict[row['raw_fastq']]['timepoint']) +' '+ sample_dict[row['raw_fastq']]['timepoint_name'], axis = 1)
        #df.drop(columns = ['index', 'time'], inplace = True)
        df.sort_values('time', inplace = True)
        labels = df['name'].to_list()
        labels = sorted(set(labels), key=labels.index)
        print(labels)

        group_df = df.groupby('name').mean()
        #group_df.sort_values('time', inplace = True)
        print(group_df.head())
        group_df.reset_index(inplace = True)
        group_df.sort_values('name', inplace = True)
        print(group_df.head())
        group_df.set_index('name', inplace = True)

        group_df.drop(columns = [ 'time'], inplace = True)
        print(group_df.head())
        #plot_df = pd.pivot_table(group_df, index='anticodon', columns='time', values='m5C fraction')
        fig, axs = plt.subplots(figsize=(5, 25), nrows=1, ncols=1)
        if 'low coverage' in  group_df.columns:
            group_df.drop(columns = ['low coverage'], inplace =True)
        for col in group_df.columns:
            if group_df[col].max()<config['abundance_heatmap_cutoff']:
                group_df.drop(columns = [col], inplace =True)
        group_df = group_df.transpose()

        sn.heatmap(group_df, ax=axs, vmin = 0, vmax =0.1, square = True, linewidths= 0.0, cmap = 'Blues')
        #fig.subplots_adjust(hspace=0.5)
        print(group_df.head())

        fig.savefig(output.anticodon_pdf, bbox_inches="tight")


rule plot_heatmap_anticodon_all_points:
    input:
        tsv = 'resources/coverage/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/coverage_summary_DM.tsv'
    output:
        anticodon_pdf = 'results/abundance/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/anticodon/heatmap-singlesample.pdf'
    run:
        import pandas as pd
        import matplotlib.pyplot as plt
        import seaborn as sn
        df = pd.read_csv(input.tsv, sep="\t",)
        df.drop(columns = [c for c in df.columns if config['abundance_score'] not in c and c != 'anticodon'], inplace = True)
        df = df.groupby('anticodon').sum()

        df = df.transpose()
        df.reset_index(inplace = True)
        df['raw_fastq'] = df.apply(lambda row: row['index'].split(' ')[0], axis = 1)
        print(df.head())
        df['time'] = df.apply(lambda row: int(sample_dict[row['raw_fastq']]['timepoint']), axis = 1)
        df['name'] = df.apply(lambda row: str(sample_dict[row['raw_fastq']]['timepoint']) + '_'+sample_dict[row['raw_fastq']]['timepoint_name']+'_'+row['raw_fastq'], axis = 1)


        #df.drop(columns = ['index', 'time'], inplace = True)
        df.sort_values('time', inplace = True)
        labels = df['name'].to_list()
        labels = sorted(set(labels), key=labels.index)
        print(labels)

        group_df = df

        group_df.drop(columns = [ 'time', 'index', 'raw_fastq' ], inplace = True)
        group_df.set_index('name', inplace = True)
        #plot_df = pd.pivot_table(group_df, index='anticodon', columns='time', values='m5C fraction')
        fig, axs = plt.subplots(figsize=(7, 25.5), nrows=1, ncols=1)
        if 'low coverage' in  group_df.columns:
            group_df.drop(columns = ['low coverage'], inplace =True)
        for col in group_df.columns:
            if group_df[col].max()<config['abundance_heatmap_cutoff']:
                group_df.drop(columns = [col], inplace =True)
        group_df = group_df.transpose()
        sn.heatmap(group_df, ax=axs, vmin = 0, vmax =0.1, square = True, linewidths= 0.0, cmap = 'Blues')
        #fig.subplots_adjust(hspace=0.5)
        print(group_df.head())

        fig.savefig(output.anticodon_pdf, bbox_inches="tight")

rule plot_abundance_boxplots_cluster:
    input:
        tsv = 'resources/coverage/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/coverage_summary_DM.tsv',
        cluster = 'resources/cluster/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-{e_cutoff}-mm-{m_cutoff}_{treatment}.yaml',
        cluster_name = 'resources/cluster/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/anticodonbased_clusternames-ed-{e_cutoff}-mm-{m_cutoff}_{treatment}.yaml'
    output:
        cluster_pdf = 'results/abundance/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/cluster/ed-{e_cutoff}-mm-{m_cutoff}_{treatment}-boxplot.pdf'
    run:
        import pandas as pd
        import yaml
        import matplotlib.pyplot as plt

        df = pd.read_csv(input.tsv, sep="\t",)
        #print(df.head())
        #df.set_index('RNAME')
        df.drop(columns = [c for c in df.columns if config['abundance_score'] not in c and c != 'RNAME'], inplace = True)
        with open(input.cluster) as file:
            cluster_dict = yaml.safe_load(file)
        with open(input.cluster_name) as file:
            cluster_name_dict = yaml.safe_load(file)
        df['cluster'] = df.apply(lambda row: cluster_dict[row['RNAME']] if row['RNAME'] in cluster_dict.keys() else 'unknown', axis =1)
        df['cluster_name'] = df.apply(lambda row: cluster_name_dict[row['cluster']] if row['cluster'] in cluster_name_dict.keys() else 'unknowN', axis =1)
        print(df.head(60))
        df.drop(columns = ['cluster', 'RNAME'], inplace = True)
        #df.reset_index(inplace = True)
        df = df.groupby('cluster_name').sum()
        df.reset_index(inplace = True)
        df.set_index('cluster_name',inplace=True)
        df = df.transpose()
        df.reset_index(inplace = True)
        df['raw_fastq'] = df.apply(lambda row: row['index'].split(' ')[0], axis = 1)
        df['time'] = df.apply(lambda row: int(sample_dict[row['raw_fastq']]['timepoint']), axis = 1)
        df['name'] = df.apply(lambda row: sample_dict[row['raw_fastq']]['timepoint_name'], axis = 1)
        #df.drop(columns = ['index', 'time'], inplace = True)
        df.sort_values('time', inplace = True)
        print(df.head(10))
        labels = df['name'].to_list()
        labels = sorted(set(labels), key=labels.index)
        print(labels)
        group_df = df.groupby('time').mean()
        group_df.reset_index(inplace = True)
        group_df.sort_values('time', inplace = True)
        group_df.set_index('time')

        for col in list(group_df.columns):
            if group_df[col].max() < 0.0001:
                group_df.drop(columns = [col], inplace = True)

        fig, axs = plt.subplots(figsize=(8, 3*len(group_df.columns)), nrows=len(group_df.columns)-1, ncols=1)
        fig.subplots_adjust(hspace=0.5)

        i = 0
        for col in list(group_df.columns):
            if col == 'time':
                #group_df.plot(c='grey', legend = False, ax = axs[i])
                continue

            axs[i].set_xticks(group_df['time'].to_list())
            axs[i].set_xticklabels(labels)


            group_df.plot(x = 'time', y=col,  ax = axs[i])
            #plot_df = df.set_index('name')
            df.boxplot(by = 'time', column = col, ax = axs[i], positions = group_df['time'].to_list(), labels = labels )
            axs[i].set_ylim((0, group_df[col].max()*1.5))
            axs[i].set_ylabel('fraction of reads')
            axs[i].set_xlabel('time points')
            #print(axs[i].get_xticklabels())
            #print(labels)
            axs[i].set_xticklabels(labels+labels)
            axs[i].get_legend().remove()

            i += 1
        fig.savefig(output.cluster_pdf, bbox_inches="tight")

rule plot_abundance_scatterplots_cluster:
    input:
        tsv = 'resources/coverage/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/coverage_summary_DM.tsv',
        cluster = 'resources/cluster/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-{e_cutoff}-mm-{m_cutoff}_{treatment}.yaml',
        cluster_name = 'resources/cluster/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/anticodonbased_clusternames-ed-{e_cutoff}-mm-{m_cutoff}_{treatment}.yaml'
    output:
        cluster_pdf = 'results/abundance/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/cluster_ed-{e_cutoff}-mm-{m_cutoff}/{treatment}_per-cluster_time_abundance_scatter/summary.pdf'
    run:
        import pandas as pd
        import yaml
        import matplotlib.pyplot as plt

        df = pd.read_csv(input.tsv, sep="\t",)
        #print(df.head())
        #df.set_index('RNAME')
        df.drop(columns = [c for c in df.columns if config['abundance_score'] not in c and c != 'RNAME'], inplace = True)
        with open(input.cluster) as file:
            cluster_dict = yaml.safe_load(file)
        with open(input.cluster_name) as file:
            cluster_name_dict = yaml.safe_load(file)
        df['cluster'] = df.apply(lambda row: cluster_dict[row['RNAME']] if row['RNAME'] in cluster_dict.keys() else 'unknown', axis =1)
        df['cluster_name'] = df.apply(lambda row: cluster_name_dict[row['cluster']] if row['cluster'] in cluster_name_dict.keys() else 'low coverage', axis =1)

        print(df.head())
        df.drop(columns = ['cluster', 'RNAME'], inplace = True)
        #df.reset_index(inplace = True)
        df = df.groupby('cluster_name').sum()
        df.reset_index(inplace = True)
        df.set_index('cluster_name',inplace=True)
        df = df.transpose()
        df.reset_index(inplace = True)
        df['raw_fastq'] = df.apply(lambda row: row['index'].split(' ')[0], axis = 1)
        df['time'] = df.apply(lambda row: int(sample_dict[row['raw_fastq']]['timepoint']), axis = 1)
        df['name'] = df.apply(lambda row: sample_dict[row['raw_fastq']]['timepoint_name'], axis = 1)
        #df.drop(columns = ['index', 'time'], inplace = True)
        df.sort_values('time', inplace = True)
        print(df.head(10))
        labels = df['name'].to_list()
        labels = sorted(set(labels), key=labels.index)
        print(labels)
        group_df = df.groupby('time').mean()
        group_df.reset_index(inplace = True)
        group_df.sort_values('time', inplace = True)
        group_df.set_index('time')



        plot_symbols = config['replicate_markes']
        df['replicate_symbol'] = df.apply(lambda row: plot_symbols[int(sample_dict[row['raw_fastq']]['replicate'])-1], axis = 1)

        pdfs = []
        plot_dir = '/'.join(output.cluster_pdf.split('/')[0:-1])

        plt.rc('font', size=7) #controls default text size
        plt.rc('axes', titlesize=7) #fontsize of the title
        plt.rc('axes', labelsize=7) #fontsize of the x and y labels
        plt.rc('xtick', labelsize=7) #fontsize of the x tick labels
        plt.rc('ytick', labelsize=7) #fontsize of the y tick labels
        plt.rc('legend', fontsize=7) #fontsize of the legend


        for col in list(group_df.columns):
            if col in ['time', 'replicate_symbol']:
                continue


            fig, axs = plt.subplots(figsize=(6*CM, 3.3*CM), nrows= 1, ncols=1, layout = 'tight')
            #group_df[col] = group_df[col]*1000000
            #df[col] = df[col]*1000000


            plt.rc('font', size=7) #controls default text size
            plt.rc('axes', titlesize=7) #fontsize of the title
            plt.rc('axes', labelsize=7) #fontsize of the x and y labels
            plt.rc('xtick', labelsize=7) #fontsize of the x tick labels
            plt.rc('ytick', labelsize=7) #fontsize of the y tick labels
            plt.rc('legend', fontsize=7) #fontsize of the legend

            group_df.plot(x = 'time', y=col,  ax = axs, c = '0.35', linewidth = 0.75)
            print(axs.get_yticks())
            ylim = group_df[col].max()*1.2
            for marker, rep_df in df.groupby('replicate_symbol'):
                if rep_df[col].max()*1.1>ylim:
                    ylim = rep_df[col].max()*1.15
                rep_df.plot.scatter(x = 'time', y = col,  ax = axs, alpha = 0.6, marker = marker, s=7, c='0.9', edgecolors='black',  facecolors= None, linewidths=0.6)

            axs.set_ylim(0, ylim)

            print(axs.get_yticks())
            axs.set_yticks(axs.get_yticks())

            axs.set_yticklabels(["{:.3f}".format(l) for l in  axs.get_yticks()])#, max_ylabel_len)
            #axs.tick_params(axis='y', which='major', pad=28)
            axs.set_ylabel('fraction of reads')
            axs.set_xlabel(None)
            axs.set_xticks(group_df['time'].to_list())
            axs.set_xticklabels([l.replace(',', ',\n') for l in labels], rotation = 90)
            axs.get_legend().remove()
            axs.set_title(col.split('[')[0].replace('_', '/'))

            fig_path = os.path.join(plot_dir, col + '.pdf')
            pdfs.append(fig_path)

            fig.savefig(fig_path)
            plt.close("all")

        pdfs.sort()
        call_args =['pdfunite'] + pdfs + [output.cluster_pdf]
        results = subprocess.run(call_args, capture_output=True)




rule plot_summary_time_cluster_abundance_heatmap:
    input:
        tsv = 'resources/coverage/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/coverage_summary_DM.tsv',
        cluster = 'resources/cluster/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-{e_cutoff}-mm-{m_cutoff}_{treatment}.yaml',
        cluster_name = 'resources/cluster/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/anticodonbased_clusternames-ed-{e_cutoff}-mm-{m_cutoff}_{treatment}.yaml'
    output:
        cluster_pdf = 'results/abundance/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/cluster_ed-{e_cutoff}-mm-{m_cutoff}/{treatment}_summary_time-cluster-abundance_heatmap.pdf'
    run:
        import pandas as pd
        import matplotlib.pyplot as plt
        import seaborn as sn
        import yaml

        df = pd.read_csv(input.tsv, sep="\t",)

        df.drop(columns = [c for c in df.columns if config['abundance_score'] not in c and c != 'RNAME'], inplace = True)
        with open(input.cluster) as file:
            cluster_dict = yaml.safe_load(file)
        with open(input.cluster_name) as file:
            cluster_name_dict = yaml.safe_load(file)
        df['cluster'] = df.apply(lambda row: cluster_dict[row['RNAME']] if row['RNAME'] in cluster_dict.keys() else 'unknown', axis =1)
        for cluster in config['clusters_to_remove']:
            df = df[df['cluster']!= cluster]
        df['cluster_name'] = df.apply(lambda row: cluster_name_dict[row['cluster']] if row['cluster'] in cluster_name_dict.keys() else 'low coverage', axis =1)

        #print(df.head())
        df.drop(columns = ['cluster', 'RNAME'], inplace = True)
        #df.reset_index(inplace = True)
        df = df.groupby('cluster_name').sum()
        #print(df.head())
        #df.reset_index(inplace = True)
        #df['cluster_name'] = df.apply(lambda row: row['cluster_name'].split('[')[0].replace('_', '/') , axis=1)
        #df.set_index('cluster_name',inplace=True)
        df = df.transpose()
        df.reset_index(inplace = True)
        df['raw_fastq'] = df.apply(lambda row: row['index'].split(' ')[0], axis = 1)
        df['time'] = df.apply(lambda row: int(sample_dict[row['raw_fastq']]['timepoint']), axis = 1)
        df['name'] = df.apply(lambda row: sample_dict[row['raw_fastq']]['timepoint_name'], axis = 1)
        #df.drop(columns = ['index', 'time'], inplace = True)
        df.sort_values('time', inplace = True)
        #print(df.head(10))
        labels = df['name'].to_list()
        labels = sorted(set(labels), key=labels.index)
        #print(labels)
        group_df = df.groupby('time').mean()
        group_df.reset_index(inplace = True)
        group_df.sort_values('time', inplace = True)
        group_df.set_index('time')

        group_df.drop(columns = [ 'time'], inplace = True)


        fig, axs = plt.subplots(figsize=(15*CM, 23*CM), nrows=1, ncols=1, layout='tight')
        if 'low coverage' in  group_df.columns:
            group_df.drop(columns = ['low coverage'], inplace =True)

        '''
        for col in group_df.columns:

            if group_df[col].max()< config['abundance_heatmap_cutoff']:
                if col in config['manual_cluster_recovery']:
                    continue
                print(col)
                group_df.drop(columns = [col], inplace =True)
        '''
        group_df = group_df.transpose()
        axs=sn.heatmap(group_df, ax=axs, vmin = 0, vmax =0.1, square = True,
                        linewidths= 0.0, cmap = 'Blues', linecolor='lightgrey',
                        xticklabels=labels, cbar_kws={'location': 'right','label': 'fraction of reads', 'fraction':0.02, 'aspect':15})
        axs.set_xticklabels(axs.get_xticklabels(),  fontsize=7)
        axs.set_ylabel(None)
        axs.set_xlabel(None)
        y_ticks = [n+0.5 for n in range(0,len(group_df))]
        axs.set_yticks(y_ticks)
        #print(group_df.index)
        axs.set_yticklabels([l.split('[')[0] for l in group_df.index],  fontsize=7)

        fig.savefig(output.cluster_pdf)


rule plot_summary_sample_cluster_abundance_heatmap:
    input:
        tsv = 'resources/coverage/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/coverage_summary_{treatment}.tsv',
        cluster = 'resources/cluster/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/clusters-ed-{e_cutoff}-mm-{m_cutoff}_DM.yaml',
        cluster_name = 'resources/cluster/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/anticodonbased_clusternames-ed-{e_cutoff}-mm-{m_cutoff}_DM.yaml'
    output:
        cluster_pdf = 'results/abundance/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/cluster_ed-{e_cutoff}-mm-{m_cutoff}/{treatment}_summary_sample-cluster-abundance_heatmap.pdf'
    run:
        import pandas as pd
        import matplotlib.pyplot as plt
        import seaborn as sn
        import yaml

        df = pd.read_csv(input.tsv, sep="\t",)
        df.drop(columns = [c for c in df.columns if config['abundance_score'] not in c and c != 'RNAME'], inplace = True)
        with open(input.cluster) as file:
            cluster_dict = yaml.safe_load(file)
        with open(input.cluster_name) as file:
            cluster_name_dict = yaml.safe_load(file)
        df['cluster'] = df.apply(lambda row: cluster_dict[row['RNAME']] if row['RNAME'] in cluster_dict.keys() else 'unknown', axis =1)
        for cluster in config['clusters_to_remove']:
            df = df[df['cluster']!= cluster]
        df['cluster_name'] = df.apply(lambda row: cluster_name_dict[row['cluster']] if row['cluster'] in cluster_name_dict.keys() else 'low coverage', axis =1)

        #print(df.head())
        df.drop(columns = ['cluster', 'RNAME'], inplace = True)
        #df.reset_index(inplace = True)
        df = df.groupby('cluster_name').sum()
        df.reset_index(inplace = True)
        df.set_index('cluster_name',inplace=True)

        df = df.transpose()
        df.reset_index(inplace = True)
        df['raw_fastq'] = df.apply(lambda row: row['index'].split(' ')[0], axis = 1)
        df['time'] = df.apply(lambda row: int(sample_dict[row['raw_fastq']]['timepoint']), axis = 1)
        df['name'] = df.apply(lambda row: str(sample_dict[row['raw_fastq']]['timepoint']) + '_'+sample_dict[row['raw_fastq']]['timepoint_name']+'_'+row['raw_fastq'], axis = 1)
        #df.drop(columns = ['index', 'time'], inplace = True)
        df.sort_values('time', inplace = True)
        labels = df['name'].to_list()
        labels = sorted(set(labels), key=labels.index)
        group_df = df
        group_df.reset_index(inplace = True)
        group_df.sort_values('time', inplace = True)
        group_df.set_index('time')

        group_df.drop(columns = [ 'time', 'raw_fastq', 'index', 'level_0'], inplace = True)
        group_df.set_index('name', inplace = True)
        #print(group_df.head())
        #plot_df = pd.pivot_table(group_df, index='anticodon', columns='time', values='m5C fraction')
        fig, axs = plt.subplots(figsize=(15*CM, 190*CM), nrows=1, ncols=1)
        if 'low coverage' in  group_df.columns:
            group_df.drop(columns = ['low coverage'], inplace =True)
        '''
        for col in group_df.columns:

            if group_df[col].max()<config['abundance_heatmap_cutoff']:
                if col in config['manual_cluster_recovery']:
                    continue
                print(col)
                group_df.drop(columns = [col], inplace =True)
        '''

        group_df = group_df.transpose()

        sn.heatmap(group_df, ax=axs, vmin = 0, vmax =0.1, square = True, linewidths= 0.0, cmap = 'Blues')
        axs.set_xticklabels(axs.get_xticklabels(),  fontsize=7)
        axs.set_ylabel(None)
        axs.set_xlabel(None)
        y_ticks = [n+0.5 for n in range(0,len(group_df))]
        axs.set_yticks(y_ticks)
        #axs.set_yticklabels([l.split('[')[0] for l in group_df.index],  fontsize=7)
        axs.set_yticklabels([l for l in group_df.index],  fontsize=7)

        fig.savefig(output.cluster_pdf, bbox_inches="tight")

rule get_all_abundance_plots:
    input:
        #'results/abundance/pre-filter_umi/selected/cluster_ed-3-mm-50/MOCK_summary_time-cluster-abundance_heatmap.pdf',
        #'results/abundance/pre-filter_umi/selected/cluster_ed-3-mm-50/BS_summary_time-cluster-abundance_heatmap.pdf',
        #'results/abundance/pre-filter_umi/selected/cluster_ed-3-mm-50/MOCK_summary_sample-cluster-abundance_heatmap.pdf',
        #'results/abundance/pre-filter_umi/selected/cluster_ed-3-mm-50/BS_summary_sample-cluster-abundance_heatmap.pdf',
        #'results/abundance/pre-filter_umi/selected/cluster_ed-3-mm-50/MOCK_per-cluster_time_abundance_scatter/summary.pdf',
        #'results/abundance/pre-filter_umi/selected/cluster_ed-3-mm-50/MOCK_per-cluster_time_abundance_line-per-ref/summary.pdf',
        #'results/abundance/pre-filter_umi/selected/cluster_ed-3-mm-50/BS_per-cluster_time_abundance_line-per-ref/summary.pdf',
        #'results/abundance/pre-filter_umi/selected/cluster_ed-3-mm-50/BS_per-cluster_time_abundance_scatter/summary.pdf',
        'results/abundance/pre-filter_umi/selected/cluster_ed-3-mm-50/DM_summary_time-cluster-abundance_heatmap.pdf',
        #'results/abundance/pre-filter_umi/selected/cluster_ed-4-mm-50/DM_summary_time-cluster-abundance_heatmap.pdf',
        'results/abundance/pre-filter_umi/selected/cluster_ed-3-mm-50/DM_summary_sample-cluster-abundance_heatmap.pdf',
        #'results/abundance/pre-filter_umi/selected/cluster_ed-4-mm-50/DM_summary_sample-cluster-abundance_heatmap.pdf',
        #'results/abundance/pre-filter_umi/selected/cluster_ed-2-mm-50/DM_summary_sample-cluster-abundance_heatmap.pdf',
        'results/abundance/pre-filter_umi/selected/cluster_ed-3-mm-50/DM_per-cluster_time_abundance_scatter/summary.pdf',
        #'results/abundance/pre-filter_umi/selected/cluster_ed-4-mm-50/DM_per-cluster_time_abundance_scatter/summary.pdf',
        'results/abundance/pre-filter_umi/selected/cluster_ed-3-mm-50/DM_per-cluster_time_abundance_line-per-ref/summary.pdf',
        'results/abundance/pre-filter_umi/selected/cluster_ed-3-mm-50/DM_per-cluster_time_abundance_line-per-ref-smaller-legend/summary.pdf',
        #'results/abundance/pre-filter_umi/selected/cluster_ed-4-mm-50/DM_per-cluster_time_abundance_line-per-ref/summary.pdf',
        'results/abundance/pre-filter_umi/selected/anticodon/all_refs/summary.pdf',
        'results/abundance/pre-filter_umi/selected/anticodon/heatmap-singlesample.pdf',
        'results/abundance/pre-filter_umi/selected/anticodon/heatmap.pdf',
        'results/abundance/pre-filter_umi/selected/anticodon/scatterplot.pdf'

    output:
        'results/abundance/done.txt'
    run:
        with open(str(output),'w') as file:
            for f in input:
                file.write(f+'\n')
