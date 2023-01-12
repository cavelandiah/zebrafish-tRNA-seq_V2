

rule get_intarna_predictions:
    input:
        genome_fa = 'resources/references/non-redundant/genomic.fa',
        mitochondrial_fa = 'resources/references/non-redundant/mitochondrial.fa'
    output:
        intarna_tsv = 'resources/probes/intarna/{probe}.tsv'
    run:
        import sys
        sys.path.append("../")
        import scripts.rri_utils as rri
        from Bio import SeqIO
        import pandas as pd
        import os

        probe_seq = probes_dict[wildcards.probe]['sequence']
        record_mit_dict = SeqIO.to_dict(SeqIO.parse(input.mitochondrial_fa, "fasta"))
        record_gen_dict = SeqIO.to_dict(SeqIO.parse(input.genome_fa, "fasta"))
        reference_dict =  {**record_mit_dict, **record_gen_dict}
        p_dir = '/scratch/maria/sequencing/Zebrafish/resources/probes'
        if not os.path.exists(p_dir):
            os.makedirs(p_dir)
        intarna_dir = '/scratch/maria/sequencing/Zebrafish/resources/probes/intarna'
        if not os.path.exists(intarna_dir):
            os.makedirs(intarna_dir)
        probe_dir = os.path.join('/scratch/maria/sequencing/Zebrafish/resources/probes/intarna/', wildcards.probe)
        if not os.path.exists(probe_dir):
            os.makedirs(probe_dir)
        probe_df = pd.DataFrame()
        for ref in reference_dict.keys():
            ref_seq = str(reference_dict[ref].seq)

            single_ref_path = os.path.join(probe_dir, ref + '.tsv')

            intarna_df = rri.run_intarna(
            ref_seq,
            probe_seq,
            id1=ref,
            id2=wildcards.probe,
            temperature = 37,
            intarna_args=config['intarna_args'],
            out_file=single_ref_path,
            )
            if len(probe_df) == 0:
                probe_df = intarna_df
            else:
                probe_df = pd.concat([probe_df, intarna_df])
        probe_df.to_csv( output.intarna_tsv, sep= '\t', index = False)


rule get_interaction_E_3_all:
    input:
        coverage_summary = 'resources/coverage/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/min_coverage_summary_DM.tsv',
        intarna_outputs = expand('resources/probes/intarna/{probe}.tsv', probe=PROBES)
    output:
        intarna_tsv = 'resources/probes/probes_intarna_summary_3.tsv'
    run:
        import pandas as pd
        df = pd.read_csv(input.coverage_summary, sep = '\t')
        df['max norm count'] = df['max random fraction']*1000000
        for c in df.columns:
            if c not in ['RNAME', 'anticodon', 'max norm count']:
                df.drop(columns = [c], inplace = True)
        for probe in PROBES:
            probe_df = pd.read_csv('resources/probes/intarna/'+probe+'.tsv', sep ='\t')
            #print(probe_df.columns)
            #probe_df = probe_df[probe_df['suboptimal']==0]
            probe_df.drop_duplicates('id1', keep = 'first', inplace = True)
            probe_df.set_index('id1', inplace = True)
            probe_len = len(probe_df['seq2'].to_list()[0])
            probe_df[probe + " paired 3'"] = probe_df.apply(lambda row: row['end2']==probe_len, axis =1)

            for c in probe_df.columns:
                if c == 'E':
                    #probe_df['full_E'] = probe_df['full_E'].fillna(0)
                    probe_df.rename(columns={'E': probe}, inplace = True)
                elif c ==probe + " paired 3'":
                    continue
                else:
                    probe_df.drop(columns = [c], inplace = True)
            print(probe_df.head())
            df=df.merge(probe_df, how='left', left_on= 'RNAME', right_index = True )
        df.to_csv(output.intarna_tsv, sep= '\t', index = False)
        print(df.head())

rule get_interaction_E_all:
    input:
        coverage_summary = 'resources/coverage/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/min_coverage_summary_DM.tsv',
        intarna_outputs = expand('resources/probes/intarna/{probe}.tsv', probe=PROBES)
    output:
        intarna_tsv = 'resources/probes/probes_intarna_summary.tsv'
    run:
        import pandas as pd
        df = pd.read_csv(input.coverage_summary, sep = '\t')
        df['max norm count'] = df['max random fraction']*1000000
        for c in df.columns:
            if c not in ['RNAME', 'anticodon', 'max norm count']:
                df.drop(columns = [c], inplace = True)
        for probe in PROBES:
            probe_df = pd.read_csv('resources/probes/intarna/'+probe+'.tsv', sep ='\t')
            #print(probe_df.columns)
            #probe_df = probe_df[probe_df['suboptimal']==0]
            probe_df.drop_duplicates('id1', keep = 'first', inplace = True)
            probe_df.set_index('id1', inplace = True)

            for c in probe_df.columns:
                if c == 'E':
                    #probe_df['full_E'] = probe_df['full_E'].fillna(0)
                    probe_df.rename(columns={'E': probe}, inplace = True)

                else:
                    probe_df.drop(columns = [c], inplace = True)
            print(probe_df.head())
            df=df.merge(probe_df, how='left', left_on= 'RNAME', right_index = True )
        df.to_csv(output.intarna_tsv, sep= '\t', index = False)
        print(df.head())

rule plot_probe_targets:
    input:
        intarna_summary_tsv = 'resources/probes/probes_intarna_summary.tsv'
    output:
        probe_targets_pdf = 'results/probes/probes_targets.pdf'
    run:
        import pandas as pd
        import os
        import matplotlib.pyplot as plt
        from matplotlib.cm import ScalarMappable
        import numpy as np
        import math

        def get_rate(deltaE):
            """Compute transition rate between two states from their free energies."""
            R = 1.98720425864083 * math.pow(10, -3)
            T = 273.15 + 37
            rate = math.exp(-deltaE / (R * T))
            return rate

        df = pd.read_csv(input.intarna_summary_tsv, sep = '\t')


        f, axs = plt.subplots(figsize=(200, 3*len(list(probes_dict.keys()))), nrows=len(list(probes_dict.keys())), ncols=1, sharex=True)

        df['width'] = df['max norm count']/1000
        width = df['width'].to_list()

        pos = [0]*len(width)
        for i in range(0, len(width)):
            pos[i] = sum(width[0:i])+width[i]/2

        print(df.head())
        print(df.columns)

        max_e = 15
        i = 0
        for probe in df.columns:
            print(i)

            if probe in ['RNAME', 'anticodon', 'width', 'max norm count']:
                continue
            #colors = [abs(c/max_e) for c in df[probe].replace(np.nan, 0).to_list() ]
            colors = [get_rate(c) for c in df[probe].replace(np.nan, 0).to_list() ]
            max_c = max(colors)
            colors = [c/max_c for c in colors]
            #print(colors)
            my_cmap = plt.cm.get_cmap('Blues')
            cols = my_cmap(colors)

            axs[i].set_ylabel(probe)

            axs[i].bar(pos, height = [1]*len(df) , width=df['width'], edgecolor=None, color = cols )#cmap='Blues')
            i +=1
        plt.xticks(pos, df['RNAME'].to_list(), rotation=90)
        f.savefig(output.probe_targets_pdf, bbox_inches="tight")
        plt.show()

rule plot_probe_targets_lin:
    input:
        intarna_summary_tsv = 'resources/probes/probes_intarna_summary.tsv'
    output:
        probe_targets_pdf = 'results/probes/probes_targets_linear_colorscale.pdf'
    run:
        import pandas as pd
        import os
        import matplotlib.pyplot as plt
        from matplotlib.cm import ScalarMappable
        import numpy as np
        import math

        def get_rate(deltaE):
            """Compute transition rate between two states from their free energies."""
            R = 1.98720425864083 * math.pow(10, -3)
            T = 273.15 + 37
            rate = math.exp(-deltaE / (R * T))
            return rate

        df = pd.read_csv(input.intarna_summary_tsv, sep = '\t')


        f, axs = plt.subplots(figsize=(200, 3*len(list(probes_dict.keys()))), nrows=len(list(probes_dict.keys())), ncols=1, sharex=True)

        df['width'] = df['max norm count']/1000
        width = df['width'].to_list()

        pos = [0]*len(width)
        for i in range(0, len(width)):
            pos[i] = sum(width[0:i])+width[i]/2

        print(df.head())
        print(df.columns)

        max_e = 35
        i = 0
        for probe in df.columns:
            print(i)

            if probe in ['RNAME', 'anticodon', 'width', 'max norm count']:
                continue
            colors = [abs(c/max_e) for c in df[probe].replace(np.nan, 0).to_list() ]
            colors = [c if c <1 else 1 for c in colors ]
            #colors = [get_rate(c) for c in df[probe].replace(np.nan, 0).to_list() ]
            #max_c = max(colors)
            #colors = [c/max_c for c in colors]
            #print(colors)
            my_cmap = plt.cm.get_cmap('Blues')
            cols = my_cmap(colors)

            axs[i].set_ylabel(probe)

            axs[i].bar(pos, height = [1]*len(df) , width=df['width'], edgecolor=None, color = cols )#cmap='Blues')
            i +=1
        plt.xticks(pos, df['RNAME'].to_list(), rotation=90)
        f.savefig(output.probe_targets_pdf, bbox_inches="tight")
        plt.show()


rule plot_probe_single_target:
    input:
        intarna_summary_tsv = 'resources/probes/probes_intarna_summary.tsv'
    output:
        probe_targets_pdf = 'results/probes/probes_targets_{target}.pdf'
    run:
        import pandas as pd
        import os
        import matplotlib.pyplot as plt
        from matplotlib.cm import ScalarMappable
        import numpy as np
        import math

        def get_rate(deltaE):
            """Compute transition rate between two states from their free energies."""
            R = 1.98720425864083 * math.pow(10, -3)
            T = 273.15 + 37
            rate = math.exp(-deltaE / (R * T))
            return rate

        df = pd.read_csv(input.intarna_summary_tsv, sep = '\t')

        rows = 0

        for key in probes_dict.keys():
            if probes_dict[key]['target'] == wildcards.target:
                rows +=1


        f, axs = plt.subplots(figsize=(200, 3*rows), nrows=rows, ncols=1, sharex=True)

        df['width'] = df['max norm count']/1000
        width = df['width'].to_list()

        pos = [0]*len(width)
        for i in range(0, len(width)):
            pos[i] = sum(width[0:i])+width[i]/2

        print(df.head())
        print(df.columns)

        max_e = 15
        i = 0
        for probe in df.columns:
            print(i)



            if probe in ['RNAME', 'anticodon', 'width', 'max norm count']:
                continue

            target = probes_dict[probe]['target']
            if target != wildcards.target:
                continue
            #colors = [abs(c/max_e) for c in df[probe].replace(np.nan, 0).to_list() ]
            colors = [get_rate(c) for c in df[probe].replace(np.nan, 0).to_list() ]
            max_c = max(colors)
            colors = [c/max_c for c in colors]
            #print(colors)
            my_cmap = plt.cm.get_cmap('Blues')
            cols = my_cmap(colors)

            axs[i].set_ylabel(probe)

            axs[i].bar(pos, height = [1]*len(df) , width=df['width'], edgecolor=None, color = cols )#cmap='Blues')
            i +=1
        plt.xticks(pos, df['RNAME'].to_list(), rotation=90)
        f.savefig(output.probe_targets_pdf, bbox_inches="tight")
        plt.show()
