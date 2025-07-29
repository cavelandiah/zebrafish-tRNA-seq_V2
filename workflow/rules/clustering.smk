# Author: Maria Waldl • code@waldl.org
# Version: 2024-01-24

rule get_ref_edit_dist:
    input:
        fasta = 'resources/references/alignment/{ref_set}_tRNAs.fa'
    output:
        edit_dist_tsv = 'resources/references/cluster/{ref_set}/edit_dist.tsv'
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        import pandas as pd
        from Bio import SeqIO

        records = list(SeqIO.parse(input.fasta, "fasta"))
        out_handle = open(output.edit_dist_tsv, 'w')
        i = 1
        for r1 in records:
            if '(' not in r1.id:
                continue
            s1 = str(r1.seq)
            for r2 in records:
                if '(' not in r2.id:
                    continue

                s2 = str(r2.seq)
                dist = 0
                for p in range(0, len(s1)):
                    if s1[p] != s2[p]:
                        dist += 1
                out_handle.write('\t'.join([r1.id, r2.id, str(dist)])+'\n')

        out_handle.close()

rule plot_dist_hist:
    input:
        edit_dist_tsv = 'resources/references/cluster/{ref_set}/edit_dist.tsv'
    output:
        edit_dist_hist = 'qc/cluster/{ref_set}/editdist_hist.pdf',
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        import pandas as pd
        import matplotlib.pyplot as plt
        df = pd.read_csv(input.edit_dist_tsv, header = None, sep = '\t', names = ['ref1', 'ref2', 'dist'])
        fig, axs = plt.subplots(figsize=(8, 9), nrows=2, ncols=1)
        df.plot.hist(column=["dist"], bins=list(range(0,105,5)), ax = axs[0])
        df.plot.hist(column=["dist"], bins=list(range(0,50)), ax = axs[1])
        axs[1].set_xlim(0,25)
        axs[0].set_xlabel('editdist')
        axs[1].set_xlabel('editdist')
        fig.savefig(output.edit_dist_hist, bbox_inches="tight")


rule multimapper_vs_dist:
    input:
        edit_dist_tsv = 'resources/references/cluster/{ref_set}/edit_dist.tsv',
        multimapper_summary_tsv = 'resources/multimappers/pre-filter_{reads_filter}/{ref_set}/multimapper_between_refs.tsv',
    output:
        plot = 'qc/cluster/pre-filter_{reads_filter}/{ref_set}/dist_multimapper.pdf',
        summary_tsv = 'resources/cluster/pre-filter_{reads_filter}/{ref_set}/dist_multimapper.tsv'
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        import pandas as pd
        import matplotlib.pyplot as plt
        import numpy as np
        df = pd.read_csv(input.edit_dist_tsv, header = None, sep = '\t', names = ['ref1', 'ref2', 'dist'])
        df['name'] = df['ref1']+ df['ref2']
        #df.drop(['r1', 'r2'], axis=1, inplace=True)
        m_df =  pd.read_csv(input.multimapper_summary_tsv, sep = '\t', usecols= ['max_norm_shared', 'ref1', 'ref2'])
        m_df['name'] = m_df['ref1']+ m_df['ref2']
        m_df.drop(['ref1', 'ref2'], axis=1, inplace=True)


        m_df = df.merge(m_df, left_on='name', right_on= 'name', how = 'left')


        fig, axs = plt.subplots(figsize=(8, 9), nrows=4, ncols=1)
        df.plot.hist(column=["dist"], bins=list(range(0,70)), ax = axs[0])
        axs[0].set_xlim(0,70)
        axs[0].set_xticks(list(range(0,75,5)))

        plot_df = m_df.groupby('dist')['max_norm_shared'].agg([np.mean, np.max, np.sum])

        plot_df.plot(y = 'mean',ax = axs[1])
        axs[1].set_xlim(0,70)
        axs[1].set_xticks(list(range(0,75,5)))


        plot_df.plot(y = 'amax',ax = axs[2])
        axs[2].set_xlim(0,70)
        axs[2].set_xticks(list(range(0,75,5)))


        plot_df.plot(y = 'sum',ax = axs[3])
        axs[3].set_xlim(0,70)
        axs[3].set_xticks(list(range(0,75,5)))


        fig.savefig(output.plot, bbox_inches="tight")



        m_df.drop(['name'], axis=1, inplace=True)
        m_df.to_csv(output.summary_tsv, sep = '\t', index= False)



rule cluster_by_editdist:
    input:
        summary_tsv = 'resources/references/cluster/{ref_set}/edit_dist.tsv',
        abundance_tsv = 'resources/coverage/pre-filter_{reads_filter}/{ref_set}/min_coverage_summary_DM.tsv'
    output:
        cluster_map = 'resources/cluster/pre-filter_{reads_filter}/{ref_set}/clusters-editdist-{cutoff}.yaml',
        cluster_names = 'resources/cluster/pre-filter_{reads_filter}/{ref_set}/clusternames-editdist-{cutoff}.yaml',
        cluster_info = 'resources/cluster/pre-filter_{reads_filter}/{ref_set}/clusternames-editdist-{cutoff}.txt',
        cluster_abundance = 'resources/cluster/pre-filter_{reads_filter}/{ref_set}/cluster-abundanceinfo-editdist-{cutoff}.txt',
    wildcard_constraints:
        cutoff="[0-5]",
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        import pandas as pd
        import yaml
        df = pd.read_csv(input.abundance_tsv, sep="\t")
        df.set_index('RNAME', inplace = True)
        for col in df.columns:
            if col != 'max random fraction':
                df.drop(columns = [col], inplace = True)
        abundance_dict = df.to_dict(orient='index')
        df = pd.read_csv(input.summary_tsv, header = None, sep = '\t', names = ['ref1', 'ref2', 'dist'])
        #print(df.head())
        references = list(set(df['ref1'].to_list()))
        references = [[ref] for ref in references]

        df.sort_values('dist', inplace=True, ascending=True, ignore_index=True)

        for i in range(0,len(df)):
            if df.at[i,'dist'] > int(wildcards.cutoff):
                continue
            ref1=df.at[i,'ref1']
            ref2=df.at[i,'ref2']
            for j,group in enumerate(references):
                if ref1 in group:
                    id1 = j
                    group1 = group
                if ref2 in group:
                    id2 =j
                    group2 = group
            if id1==id2:
                continue
            del references[max(id1,id2)]
            del references[min(id1,id2)]
            references.append(group1+group2)

        print(len(df))
        cluster_dict = {}
        cluster_names = {}
        txt_file = open(output.cluster_info, 'w')
        print(len(references))


        i = 0
        abundance_file_handle = open(output.cluster_abundance, 'w')
        for cluster in references:

            #if len(cluster) ==1:
                #print('single')
                #if cluster[0] not in abundance_dict.keys():
                #    continue
                #elif abundance_dict[cluster[0]]['max random fraction'] < 0.00005:
                #    print('low abundance')
                #    continue
            abundance_file_handle.write('=== Cluster: ' + str(i)+ ' ===\n')
            name_info = [[r.split('(')[-1].split(')')[0], r] for r in cluster]
            name_df = pd.DataFrame(name_info, columns = ['anticodon', 'ref'])
            name_df['abundance'] = name_df.apply(lambda row:  abundance_dict[row['ref']]['max random fraction'] if row['ref'] in abundance_dict.keys() else 0, axis =1 )

            name_df.sort_values(by = 'abundance', ascending= False, inplace = True)
            abundance_file_handle.write('---single ref abundance---\n')
            name_df.to_csv(abundance_file_handle,sep='\t', mode = 'a', index = False)
            codons = list(set(name_df['anticodon'].to_list()))
            name_df = name_df.groupby('anticodon').sum()
            name_df.sort_values(by = 'abundance', ascending= False, inplace = True)
            #print(name_df.head(21))
            abundance_file_handle.write('---anticodon abundance---\n')
            total_abundance = name_df['abundance'].sum()
            name_df['cluster_fraction'] = name_df['abundance']/total_abundance
            name_df.to_csv(abundance_file_handle,sep='\t', mode = 'a')

            # get cluster name
            name_df = name_df[name_df['cluster_fraction']>=0.0100]
            name_df.sort_values(by='cluster_fraction', inplace=True, ascending= False)
            name_df.reset_index(inplace = True)
            print(name_df.head())
            anticodons_str = '_'.join(name_df['anticodon'].to_list())
            fractions_str = '_'.join(["{fraction:.2f}".format(fraction = f)  for f in name_df['cluster_fraction'].to_list()])
            cluster_name = anticodons_str + '[' +fractions_str +']('+"{abundance:.3f}".format(abundance= total_abundance)+')'

            #cluster_name=str(max_anticodon)+'(' + str(i) +')' + '.'.join(codons)+'[' + str(int(name_df['abundance'].max()*1000000)) + '_'+"{:.2f}".format(name_df['abundance'].max()*1000000/int(name_df['abundance'].sum()*1000000))+']'
            cluster_names[i] = cluster_name
            abundance_file_handle.write('---max anticodon---\n'+cluster_name+'\n\n')
            txt_file.write('Cluster ' + str(i)+':\n')
            for r in cluster:
                txt_file.write('    ' + r +'\n')
                cluster_dict[r] = i
            i +=1
        txt_file.close()
        abundance_file_handle.close()
        with open(output.cluster_names, 'w') as file:
            out = yaml.dump(cluster_names, file)
        with open(output.cluster_map, 'w') as file:
            out = yaml.dump(cluster_dict, file)
        print('ref count: ', len(references), i)

def get_names_seqs(cluster, df_cluster_seq_names):
    return df_cluster_seq_names[df_cluster_seq_names['editdist_cluster'] == int(cluster)]['SeqName'].to_list()

def get_fasta(cluster_id1, cluster_id2, cluster1_names, cluster2_names, fasta):
    cluster_1 = []
    cluster_2 = []
    cluster_combined = []
    for record in SeqIO.parse(fasta, 'fasta'):
        if record.id in cluster1_names:
            cluster_1.append(record)
            cluster_combined.append(record)
        elif record.id in cluster2_names:
            cluster_2.append(record)
            cluster_combined.append(record)
        else:
            sys.exit('Sequence name is not found')
    final_path = "resources/cluster/cluster_subalignments"
    out1 = f"{final_path}/{cluster_id1}.fa"
    out2 = f"{final_path}/{cluster_id2}.fa"
    out3 = f"{final_path}/combined_{cluster_id1}_{cluster_id2}.fa"
    with open(out1, 'w') as out1fasta:
        SeqIO.write(cluster1, out1fasta, 'fasta')
    with open(out2, 'w') as out2fasta:
        SeqIO.write(cluster2, out2fasta, 'fasta')
    with open(out3, 'w') as out3fasta:
        SeqIO.write(cluster_combined, out3fasta, 'fasta')
    fasta_sequences = [out1, out2, out3]
    return fasta_sequences

def calculate_aln_rfam_model(fasta_files, cm):
    aln_files = []
    for seq in fasta_files:
        aln_file = seq.replace('.fa','.sto')
        call_arg = ['cmalign', cm, seq]
        f = open(aln, 'w')
        subprocess.call(call_arg, stdout=f, stderr=subprocess.STDOUT)
        f.close()
        aln_files.append(aln_file)
    return aln_files

def get_stats_alignment(aln_list):
    rinfo_files = []
    for aln in aln_list:
        rinfofile = aln.replace('.sto','.rinfo')
        # Get stats
        call_stats = ['esl-alistat', '--rinfo', rinfofile, aln]
        subprocess.run(call_stats)
        rinfo_files.append(rinfofile)
    return rinfo_files

def parse_rinfo(inputf):
    """Parse rinfo and decide
    """
    rows = []
    with open(inputf, 'r') as input:
        number_sequences = 0
        for line in input:
            line = line.strip()
            if line.startswith('#'):
                match = re.search(r'Number of sequences:', line)
                if match:
                    number_sequences = int(line.split(":")[-1].replace(' ',''))
                continue
            #   rfpos   alnpos      numres   freqres      numgap   freqgap
            data = line.split()
            if len(data) < 6:
                continue
            rfpos,alnpos,numres,freqres,numgap,freqgap = data[-6:]
            rows.append({
                'rfpos': str(rfpos),
                'alnpos' : alnpos,
                'numres': numres,
                'freqres': freqres,
                'numgap': numgap,
                'freqgap': float(freqgap),
                'numSeq': number_sequences,
                })
    df = pd.DataFrame(rows)
    # Number of alignment positions with >0.5 gaps ~ 50% or more sequences have gaps.
    count = sum(df['freqgap'] >= 0.5)
    # If the alignment contains such positions, it is not compatible.
    if count > 0:
        return True
    else:
        return False

def compatible_alignment(cluster1, cluster2, fastaf, cm, df_cluster_seq_names):
    # SeqName  editdist_cluster
    #0 Ala-AGC-1(Ala-AGC)[H]               286
    df_cluster_seq_names = df_cluster_seq_names.reset_index(names=['SeqName'])
    # index->seq names
    cluster1_seq_names = get_names_seqs(cluster1, df_cluster_seq_names)  
    cluster2_seq_names = get_names_seqs(cluster1, df_cluster_seq_names)  
    # seq names->fasta
    fasta_list = get_fasta(cluster1, cluster2, cluster1_seq_names, cluster2_seq_names, fastaf)
    # fasta->sto cmalign
    aln_files = calculate_aln_rfam_model(fasta_list, cm)
    # Aln -> stats
    rinfo_flist = get_stats_alignment(aln_files)
    compatible = []
    for rinfo in rinfo_flist:
        cluster_name = os.path.basename(rinfo).replace('.rinfo','')
        compatible = parse_rinfo(rinfo)
        compatible.append({
            'cluster_name': cluster_name,
            'compatible': compatible,
            })
    dff = pd.DataFrame(compatible)
    # Get ^combined_ to ask if the aln is compatible or not
    result = dff[dff['cluster_name'].startswith('combined_')]['compatible']
    return result


#rule multimapper_cluster_rnahaydn:
    #input:
        #fasta = 'resources/references/{ref_set}_tRNAs.fa'
        #cm = config['rfam_cm'],
        #tsv = 'resources/multimappers/pre-filter_{reads_filter}/{ref_set}/cluster-editdist-{e_cutoff}_multimappers_between_clusters_{treatment}.tsv',
        #cluster_names = 'resources/cluster/pre-filter_{reads_filter}/{ref_set}/clusters-editdist-{e_cutoff}.yaml',
        #abundance_tsv = 'resources/coverage/pre-filter_{reads_filter}/{ref_set}/min_coverage_summary_DM.tsv',
    #output:
        #cluster_map = 'resources/cluster/pre-filter_{reads_filter}/{ref_set}/clusters-ed-{e_cutoff}-mm-{m_cutoff}_{treatment}.yaml',
        #cluster_names = 'resources/cluster/pre-filter_{reads_filter}/{ref_set}/clusternames-ed-{e_cutoff}-mm-{m_cutoff}_{treatment}.yaml',
        #cluster_abundance = 'resources/cluster/pre-filter_{reads_filter}/{ref_set}/cluster-abundanceinfo-ed-{e_cutoff}-mm-{m_cutoff}_{treatment}.txt',
    #wildcard_constraints:
        #e_cutoff="[0-5]",
        #m_cutoff ="[0-9]+",
        #treatment = "[A-Z]+"
    #run:
        #import subprocess
        #import os
        #import pandas as pd
        #import re
        #import sys
        #import yaml
        
        #from Bio import SeqIO

        ## get original clusters
        #with open(input.cluster_names) as file:
            #cluster_name_dict = yaml.safe_load(file)

        #cluster_df = pd.DataFrame.from_dict(cluster_name_dict, orient='index', columns = ['editdist_cluster'])
        #cluster_df_clean = pd.DataFrame.from_dict(cluster_name_dict, orient='index', columns = ['editdist_cluster'])
        #cluster_df = cluster_df.astype({'editdist_cluster': 'int'})

        ## abundance per ref
        #usecol=['RNAME','max random fraction']
        #df = pd.read_csv(input.abundance_tsv, sep="\t",usecols=usecol)
        #df.set_index('RNAME', inplace = True)
        #cluster_df = cluster_df.merge(df, how='outer', left_index=True, right_index = True)
        #cluster_df = cluster_df.dropna(subset=['editdist_cluster','max random fraction'])
        #colums_to_use = ['cluster1','cluster2','mean_ref1_multimaper_fraction']
        #df = pd.read_csv(input.tsv, sep="\t",
                         #keep_default_na=True,
                         #na_values=['Nan','NA'],
                         #skip_blank_lines=True,
                         #dtype={'cluster1':'Int64',
                                #'cluster2':'Int64',
                                #'mean_ref1_multimaper_fraction':'float64'
                                #},
                         #usecols=colums_to_use)

        #cluster_count = len(list(set(df['cluster1'].to_list())))
        #all_clusters = list(df['cluster1'].dropna().unique().tolist())
        ## Clusters in lists [c]
        #clusters = [[str(c)] for c in all_clusters if str(c)!='nan' and str(c)!='Nan']
        ## Selection of data based on the m_cutoff threshold, for an undocumented reason is set to 50
        #df = df[df['mean_ref1_multimaper_fraction']>=int(wildcards.m_cutoff)/100]
        #df = df.dropna(subset=['cluster1','cluster2'])
        #df = df.astype({'cluster1': 'str', 'cluster2': 'str'})
        #df = df[df['cluster1']!=df['cluster2']]
        #df = df.query("cluster1 != cluster2")

        ## Potential multimapper clusters should be evaluated in terms of the alignment
        #df['Compatible'] = df.apply(lambda row: compatible_alignment(row['cluster1', 'cluster2'], input.fasta, input.cm, cluster_df_clean), axis=1)
        ## Selection of compatible merged clusters
        #df = df[df['Compatible'] == True]
        
        ## TODO: Documentation?
        #for c1, gdf in df.groupby('cluster1'):
            #max_shared_cluster = gdf.sort_values('mean_ref1_multimaper_fraction', ascending=False, inplace=False)['cluster2'].to_list()[0]
            #for j,cluster in enumerate(clusters):
                ## Here make a match of indexes
                #if c1 in cluster:
                    #id1 = j
                    #cluster1 = cluster
                #if max_shared_cluster in cluster:
                    #id2 = j
                    #cluster2 = cluster

            #if id1 == id2:
                ## Only one deletion is needed if both indices are the same.
                #if id1 < len(clusters):
                    #del clusters[id1]
                #else:
                    #print(f"Index error: id1 {id1} out of range for clusters of length {len(clusters)}")
            #else:
                ## Remove the element at the higher index first.
                #to_remove_max = max(id1, id2)
                #to_remove_min = min(id1, id2)
                #if to_remove_max < len(clusters) and to_remove_min < len(clusters):
                    #del clusters[to_remove_max]
                    #del clusters[to_remove_min]
                #else:
                    #print(f"Index error: to_remove_max {to_remove_max} or to_remove_min {to_remove_min} out of range for clusters of length {len(clusters)}")

            ## After deletion, merge the clusters
            #together = list(set(cluster1 + cluster2))
            #if together not in clusters:
                #clusters.append(together)
        
        ## Corrected version 
        ## 1. Build a dict: cluster_id_str → cluster_index
        #mapping = {}
        #for idx, cluster in enumerate(clusters):
            #for member in cluster:
                #mapping[member] = idx

        ## 2. Convert & map in two vectorized steps
        #cluster_df['clust_str'] = (
            #cluster_df['editdist_cluster']
            #.dropna()
            #.astype(int)
            #.astype(str)
        #)
        #cluster_df['cluster'] = cluster_df['clust_str'].map(mapping)

        ## write output
        #cluster_names = {}
        #abundance_file_handle = open(output.cluster_abundance, 'w')

        #for cluster, gdf in cluster_df.groupby('cluster'):
            #gdf.reset_index(inplace=True)
            #gdf.sort_values('max random fraction', ascending=False, inplace=True)
            #gdf['anticodonL'] = gdf.apply(lambda row: anticodonfamily_from_rname(row['index']), axis =1)
            #total_abundance = gdf['max random fraction'].sum()
            #gdf['cluster_fraction'] = gdf['max random fraction']/total_abundance
            ##print(gdf.head())

            #abundance_file_handle.write('=== Cluster: ' + str(cluster)+ ' ===\n')

            #a_gdf = gdf.groupby('anticodonL').sum().sort_values('max random fraction', ascending=False).reset_index()
            ##print(a_gdf)
            #anticodons_str= '_'.join(a_gdf['anticodonL'].to_list())
            #fractions_str = '_'.join(["{fraction:.2f}".format(fraction = f)  for f in a_gdf['cluster_fraction'].to_list()])
            #cluster_name = anticodons_str + '[' +fractions_str +']('+"{abundance:.3f}".format(abundance= total_abundance)+')'
            #print(cluster_name)
            #cluster_names[cluster] = cluster_name

            #abundance_file_handle.write('---single ref abundance---\n')
            #gdf.to_csv(abundance_file_handle,sep='\t', mode = 'a', index = False)

            #abundance_file_handle.write('---anticodon abundance---\n')
            #a_gdf.to_csv(abundance_file_handle,sep='\t', mode = 'a')

        #abundance_file_handle.close()

        #cluster_dict = cluster_df.to_dict()['cluster']
        #with open(output.cluster_map, 'w') as file:
            #out = yaml.dump(cluster_dict, file)
        #with open(output.cluster_names, 'w') as file:
            #out = yaml.dump(cluster_names, file)

rule multimapper_cluster:
    input:
        # multimappers.smk: summarize_multimappers_between_editdist_clusters
        tsv = 'resources/multimappers/pre-filter_{reads_filter}/{ref_set}/cluster-editdist-{e_cutoff}_multimappers_between_clusters_{treatment}.tsv',
        cluster_names = 'resources/cluster/pre-filter_{reads_filter}/{ref_set}/clusters-editdist-{e_cutoff}.yaml',
        abundance_tsv = 'resources/coverage/pre-filter_{reads_filter}/{ref_set}/min_coverage_summary_DM.tsv',
    output:
        cluster_map = 'resources/cluster/pre-filter_{reads_filter}/{ref_set}/clusters-ed-{e_cutoff}-mm-{m_cutoff}_{treatment}.yaml',
        cluster_names = 'resources/cluster/pre-filter_{reads_filter}/{ref_set}/clusternames-ed-{e_cutoff}-mm-{m_cutoff}_{treatment}.yaml',
        cluster_abundance = 'resources/cluster/pre-filter_{reads_filter}/{ref_set}/cluster-abundanceinfo-ed-{e_cutoff}-mm-{m_cutoff}_{treatment}.txt',
    wildcard_constraints:
        e_cutoff="[0-5]",
        m_cutoff ="[0-9]+",
        treatment = "[A-Z]+"
    run:
        import pandas as pd
        import yaml


        # get original clusters
        with open(input.cluster_names) as file:
            cluster_name_dict = yaml.safe_load(file)

        cluster_df = pd.DataFrame.from_dict(cluster_name_dict, orient='index', columns = ['editdist_cluster'])
        cluster_df = cluster_df.astype({'editdist_cluster': 'int'})

        # abundance per ref
        usecol=['RNAME','max random fraction']
        df = pd.read_csv(input.abundance_tsv, sep="\t",usecols=usecol)
        df.set_index('RNAME', inplace = True)

        cluster_df = cluster_df.merge(df, how='outer', left_index=True, right_index = True)
        cluster_df = cluster_df.dropna(subset=['editdist_cluster','max random fraction'])

        colums_to_use = ['cluster1','cluster2','mean_ref1_multimaper_fraction']
        df = pd.read_csv(input.tsv, sep="\t",
                         keep_default_na=True,
                         na_values=['Nan','NA'],
                         skip_blank_lines=True,
                         dtype={'cluster1':'Int64',
                                'cluster2':'Int64',
                                'mean_ref1_multimaper_fraction':'float64'
                                },
                         usecols=colums_to_use)
        
        #df = pd.read_csv(input.abundance_tsv, sep="\t")
        #df.set_index('RNAME', inplace = True)
        #for col in df.columns:
            #if col != 'max random fraction':
                #df.drop(columns = [col], inplace = True)
        #cluster_df = cluster_df.merge(df, how='outer', left_index=True, right_index = True )
        
        #df = pd.read_csv(input.tsv, sep="\t")
        cluster_count = len(list(set(df['cluster1'].to_list())))
        all_clusters = list(df['cluster1'].dropna().unique().tolist())
        # Clusters in lists [c]
        clusters = [[str(c)] for c in all_clusters if str(c)!='nan' and str(c)!='Nan']

        # Selection of data
        df = df[df['mean_ref1_multimaper_fraction']>=int(wildcards.m_cutoff)/100]
        df = df.dropna(subset=['cluster1','cluster2'])
        ## Selection of data
        #df = df[df['mean_ref1_multimaper_fraction']>=int(wildcards.m_cutoff)/100]
        ## Drop all Nan in cluster1
        #df = df[df['cluster1']==df['cluster1']]
        #df = df[df['cluster1']!='Nan']

        #df = df[df['cluster2']==df['cluster2']]
        #df = df[df['cluster2']!='Nan']

        df = df.astype({'cluster1': 'str', 'cluster2': 'str'})
        df = df[df['cluster1']!=df['cluster2']]
        df = df.query("cluster1 != cluster2")

        print(df.head(40))

        for c1, gdf in df.groupby('cluster1'):
            #print(c1)
            #print("Available clusters: "+str(clusters))
            #print(gdf)


            max_shared_cluster = gdf.sort_values('mean_ref1_multimaper_fraction', ascending=False, inplace=False)['cluster2'].to_list()[0]
            #print('max', max_shared_cluster)
            for j,cluster in enumerate(clusters):
                # Here make a match of indexes
                if c1 in cluster:
                    id1 = j
                    cluster1 = cluster
                if max_shared_cluster in cluster:
                    id2 = j
                    cluster2 = cluster
            #print("index1:"+str(id1)+" index2:"+str(id2))
            #print(clusters[id1])
            #print(clusters[id2])
            #print("#############\n")

            if id1 == id2:
                # Only one deletion is needed if both indices are the same.
                if id1 < len(clusters):
                    del clusters[id1]
                else:
                    print(f"Index error: id1 {id1} out of range for clusters of length {len(clusters)}")
            else:
                # Remove the element at the higher index first.
                to_remove_max = max(id1, id2)
                to_remove_min = min(id1, id2)
                if to_remove_max < len(clusters) and to_remove_min < len(clusters):
                    del clusters[to_remove_max]
                    del clusters[to_remove_min]
                else:
                    print(f"Index error: to_remove_max {to_remove_max} or to_remove_min {to_remove_min} out of range for clusters of length {len(clusters)}")

            # After deletion, merge the clusters
            together = list(set(cluster1 + cluster2))
            if together not in clusters:
                clusters.append(together)

        ##################get merged clusters
        #df = pd.read_csv(input.tsv, sep="\t")

        #cluster_count = len(list(set(df['cluster1'].to_list())))
        #clusters = list(set(df['cluster1'].to_list()))
        #clusters = [[str(c)] for c in clusters if str(c)!='nan' and str(c)!='Nan']

        #df = df[df['mean_ref1_multimaper_fraction']>=int(wildcards.m_cutoff)/100]
        #df = df[df['cluster1']==df['cluster1']]
        #df = df[df['cluster1']!='Nan']

        #df = df[df['cluster2']==df['cluster2']]
        #df = df[df['cluster2']!='Nan']

        #df = df.astype({'cluster1': 'str', 'cluster2': 'str'})
        #df = df[df['cluster1']!=df['cluster2']]
        #df = df.query("cluster1 != cluster2")

        #print(df.head(40))

        #for c1, gdf in df.groupby('cluster1'):
            #print(c1)
            #print(clusters)
            #print(gdf)


            #max_shared_cluster = gdf.sort_values('mean_ref1_multimaper_fraction', ascending=False, inplace=False)['cluster2'].to_list()[0]
            #print('max', max_shared_cluster)
            #for j,cluster in enumerate(clusters):
                #if c1 in cluster:
                    #id1 = j
                    #cluster1 = cluster
                #if max_shared_cluster in cluster:
                    #id2 = j
                    #cluster2 = cluster
            #print(id1,id2)
            #print(clusters[id1])
            #print(clusters[id2])

            #del clusters[max(id1,id2)]
            #del clusters[min(id1,id2)]
            #clusters.append(list(set(cluster1+cluster2)))

        # annotate new clusters
        #print(cluster_df.head(20))
        #print(clusters)
        
        # TODO:This part is not used
        #for e_clust in cluster_df['editdist_cluster'].to_list():
            #if pd.isna(e_clust):
                #continue
            #idx_str = str(int(e_clust))
            #match = [i for i, c in enumerate(clusters) if idx_str in c]
            #if not match:
                #print('m', e_clust)

        # TODO: It has a problem when there is no reads mapping here: Nan
        #cluster_df['cluster'] = cluster_df.apply(lambda row: [i for i, c in enumerate(clusters) if str(int(row['editdist_cluster'])) in c][0] , axis = 1)
        # Corrected version 
        # 1. Build a dict: cluster_id_str → cluster_index
        mapping = {}
        for idx, cluster in enumerate(clusters):
            for member in cluster:
                mapping[member] = idx

        # 2. Convert & map in two vectorized steps
        cluster_df['clust_str'] = (
            cluster_df['editdist_cluster']
            .dropna()
            .astype(int)
            .astype(str)
        )
        cluster_df['cluster'] = cluster_df['clust_str'].map(mapping)

        # write output

        cluster_names = {}
        abundance_file_handle = open(output.cluster_abundance, 'w')

        for cluster, gdf in cluster_df.groupby('cluster'):
            gdf.reset_index(inplace=True)
            gdf.sort_values('max random fraction', ascending=False, inplace=True)
            gdf['anticodonL'] = gdf.apply(lambda row: anticodonfamily_from_rname(row['index']), axis =1)
            total_abundance = gdf['max random fraction'].sum()
            gdf['cluster_fraction'] = gdf['max random fraction']/total_abundance
            #print(gdf.head())

            abundance_file_handle.write('=== Cluster: ' + str(cluster)+ ' ===\n')

            a_gdf = gdf.groupby('anticodonL').sum().sort_values('max random fraction', ascending=False).reset_index()
            #print(a_gdf)
            anticodons_str= '_'.join(a_gdf['anticodonL'].to_list())
            fractions_str = '_'.join(["{fraction:.2f}".format(fraction = f)  for f in a_gdf['cluster_fraction'].to_list()])
            cluster_name = anticodons_str + '[' +fractions_str +']('+"{abundance:.3f}".format(abundance= total_abundance)+')'
            print(cluster_name)
            cluster_names[cluster] = cluster_name

            abundance_file_handle.write('---single ref abundance---\n')
            gdf.to_csv(abundance_file_handle,sep='\t', mode = 'a', index = False)

            abundance_file_handle.write('---anticodon abundance---\n')
            a_gdf.to_csv(abundance_file_handle,sep='\t', mode = 'a')

        abundance_file_handle.close()

        cluster_dict = cluster_df.to_dict()['cluster']
        with open(output.cluster_map, 'w') as file:
            out = yaml.dump(cluster_dict, file)

        with open(output.cluster_names, 'w') as file:
            out = yaml.dump(cluster_names, file)

rule get_cluster_names:
    input:
        #mismatch_tsvs = expand('resources/coverage_counts/pre-filter_{reads_filter}/{ref_set}/clusters-ed-{e_cutoff}-mm-{m_cutoff}_{{treatment}}/{sample}_per_cluster.tsv', reads_filter = '{reads_filter}', ref_set='{ref_set}', sample = SAMPLES, e_cutoff='{e_cutoff}', m_cutoff='{m_cutoff}'),
        mismatch_tsvs = expand('resources/coverage_counts/pre-filter_{reads_filter}/{ref_set}/clusters-ed-{e_cutoff}-mm-{m_cutoff}_DM/{sample}_per_cluster.tsv', reads_filter = '{reads_filter}', ref_set='{ref_set}', sample = dm_samples, e_cutoff='{e_cutoff}', m_cutoff='{m_cutoff}' ),
        cluster_names = 'resources/cluster/pre-filter_{reads_filter}/{ref_set}/clusternames-ed-{e_cutoff}-mm-{m_cutoff}_{treatment}.yaml',
    output:
        cluster_names = 'resources/cluster/pre-filter_{reads_filter}/{ref_set}/anticodonbased_clusternames-ed-{e_cutoff}-mm-{m_cutoff}_{treatment}.yaml',
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        import pandas as pd
        import yaml

        with open(input.cluster_names) as file:
            cluster_name_dict = yaml.safe_load(file)

        data = []
        for tsv_file in input.mismatch_tsvs:
            print(tsv_file)
            sample = tsv_file.split('/')[-1].replace('_per_cluster.tsv', '')
            #sample = tsv_file.split('/')[-1].replace('.tsv', '').split('_')[0]
            print(sample_dict)
            if sample_dict[sample]['treatment'] != wildcards.treatment:
                continue
            df = pd.read_csv(tsv_file, sep = '\t')
            df['sample'] = sample
            df['timepoint'] = sample_dict[sample]['timepoint']
            df['time point name']  = sample_dict[sample]['timepoint_name']
            keep_cols = ['ref_T', 'ref_G', 'ref_C', 'ref_A', 'canonical_pos', 'timepoint','cluster', 'all nts']
            for c in df.columns:
                if c in keep_cols:
                    continue
                df.drop(columns= [c], inplace = True)

            for pos in ['34', '35', '36']:
                df_p= df[df['canonical_pos'] == pos].copy()
                df_p['ref_T'] = df_p['ref_T']/df_p['all nts']
                df_p['ref_A'] = df_p['ref_A']/df_p['all nts']
                df_p['ref_G'] = df_p['ref_G']/df_p['all nts']
                df_p['ref_C'] = df_p['ref_C']/df_p['all nts']
                df_p.drop(columns= ['all nts'], inplace = True)
                data.append(df_p)
        df =  pd.concat(data)
        df = df.groupby(['canonical_pos','cluster','timepoint']).mean()
        df.reset_index(inplace = True)
        df = df.groupby(['canonical_pos','cluster']).mean()
        df.drop(columns= ['timepoint'], inplace = True)
        df.reset_index(inplace = True)

        cluster_names = {}
        # The problem is with Clusters that were annotated as Nan
        # as a product of per_ref_nt_count in coverage.smk
        # Get rid of rows with 'Nan' in cluster column
        mask = df['cluster'].notna() & (df['cluster'] != 'Nan')
        df = df[mask]
        
        for cluster, gdf in df.groupby('cluster'):
            print(f'{cluster}')

            cluster_name = None
            if str(cluster) == "Nan":
                print(f"{cluster} with NAn")
                cluster_name = 'Nan'
            else:
                if int(cluster) in cluster_name_dict:
                    cluster_name = cluster_name_dict[int(cluster)]
                    print(cluster_name)
                else:
                    print("DONOTEXITS!")
                    sys.exit()
            codon_aa_dict = {}
            codon_info = cluster_name.split('[')[0]
            codon_info = codon_info.split('_')
            for codon in codon_info:
                codon_aa_dict[codon.split('-')[-1]] = codon


            gdf.set_index('canonical_pos', inplace=True)
            anticodons = []
            p = '34'
            pos34 = [['A',round(gdf.at[p,'ref_A'],2)], ['C',round(gdf.at[p,'ref_C'],2)], ['G',round(gdf.at[p,'ref_G'],2)],['T',round(gdf.at[p,'ref_T'],2)]]
            pos34 = [p for p in pos34 if p[1]>=0.01]
            p = '35'
            pos35 = [['A',round(gdf.at[p,'ref_A'],2)], ['C',round(gdf.at[p,'ref_C'],2)], ['G',round(gdf.at[p,'ref_G'],2)],['T',round(gdf.at[p,'ref_T'],2)]]
            pos35 = [p for p in pos35 if p[1]>=0.01]
            p = '36'
            pos36 = [['A',round(gdf.at[p,'ref_A'],2)], ['C',round(gdf.at[p,'ref_C'],2)], ['G',round(gdf.at[p,'ref_G'],2)],['T',round(gdf.at[p,'ref_T'],2)]]
            pos36 = [p for p in pos36 if p[1]>=0.01]
            for nt34 in pos34:
                for nt35 in pos35:
                    for nt36 in pos36:
                        abundance = nt34[1]*nt35[1]*nt36[1]
                        if abundance < 0.01:
                            continue
                        anticodon = nt34[0]+nt35[0]+nt36[0]
                        #print(anticodon)
                        anticodons.append([codon_aa_dict[anticodon], abundance])

            anticodons  =  sorted(anticodons, key = lambda l: l[1], reverse=True)
            name = '_'.join([n[0] for n in anticodons]) + '[' +'_'.join([str(n[1]) for n in anticodons])+'](' + str(cluster) +')'
            cluster_names[cluster] = name
        with open(output.cluster_names, 'w') as file:
            out = yaml.dump(cluster_names, file)
