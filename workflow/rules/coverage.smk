# Author: Maria Waldl • code@waldl.org
# Version: 2024-01-24

rule per_ref_nt_count:
    input:
        sam = 'resources/filtered-mappings/pre-filter_{reads_filter}/{ref_set}/random/{sample}.sam',
        ref_fasta = 'resources/references/{ref_set}_refs.fa',
        alignment_fa = 'resources/references/alignment/{ref_set}_tRNAs.fa',
        cluster_map = 'resources/cluster/pre-filter_{reads_filter}/{ref_set}/clusters-ed-{e_cutoff}-mm-{m_cutoff}_{c_treatment}.yaml',
    output:
        tsv =  'resources/coverage_counts/pre-filter_{reads_filter}/{ref_set}/clusters-ed-{e_cutoff}-mm-{m_cutoff}_{c_treatment}/{sample}.tsv'
    params:
        sample = lambda wc: wc.get('sample'),
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        import pandas as pd
        import os
        from Bio import SeqIO
        import yaml
        
        def correct_cigar(cigar):
            """
            context: It assumes that tRNA reference is composed as: tRNA+(N's). When added tail N's, segemehl
            reports cigar strings with mismatches (X) and not insertions (I) at the end of the sequence. Just
            replace last X by I.
            input: CIGAR string
            output: CIGAR modified string
            """
            if cigar.endswith('X'):
                new_cigar = cigar
                new_cigar = re.sub(r'X$', 'I', new_cigar)
                return new_cigar
            else:
                return cigar

        def get_matched_seq(cigar, seq):
            # Correct 'I' cigar after getting 'N' on ref.
            #cigar = correct_cigar(cigar)
            block_types = [i for i in cigar if not i.isdigit()]
            block_sizes = re.split("I|D|M|S|=|X", cigar)[0:-1]
            block_sizes = [int(s) for s in block_sizes]
            match_seq = ""
            l = 0  # length of ungapped matched seq
            for i, t in enumerate(block_types):
                if t == "S":
                    seq = seq[block_sizes[i] :]
                elif t == "M":
                    match_seq += seq[0 : block_sizes[i]].upper()
                    seq = seq[block_sizes[i] :]
                    l += block_sizes[i]
                elif t == "=":
                    match_seq += seq[0 : block_sizes[i]].upper()
                    seq = seq[block_sizes[i] :]
                    l += block_sizes[i]
                elif t == "X":
                    match_seq += seq[0 : block_sizes[i]].lower()
                    seq = seq[block_sizes[i] :]
                    l += block_sizes[i]
                elif t == "I":
                    seq = seq[block_sizes[i] :]
                    l += block_sizes[i]
                elif t == "D":
                    match_seq += "-" * block_sizes[i]
            return match_seq, l

        def get_padded_seq(matched_seq, pos, lead, tail, seq):
            lead_seq = seq[0:tail].lower()
            tail_seq = seq[len(seq) - tail + 1 :].lower()
            # insertions = len(seq)-lead-tail-len(matched_seq)

            padding_len = pos - 1

            diff = padding_len - len(lead_seq)

            if diff <= 0:
                lead_seq = lead_seq[(diff * (-1)) :]
            else:
                lead_seq = " " * (diff - 1) + ">" + lead_seq

            padded_seq = lead_seq + matched_seq + tail_seq + "<"
            return padded_seq

        def match_cm_position(alignment_fa):
            handle = alignment_fa
            pos_dict = {}
            for record in SeqIO.parse(handle, "fasta"):
                rname = record.id.replace('/', '_')
                seq = str(record.seq)
                pos_list = [i for i,s in enumerate(seq) if s not in  '.-' ]
                pos_dict[rname] = (pos_list, len(seq), seq.replace('.', '').replace('-', ''))
            return pos_dict

        def alignment_to_canonical_positions_mapper(alignment_fa):
            a_to_c_dict= {}

            #get subscripted positions
            #X14835.1/6927-7002 contains subscripted canonical positions at given zero based index
            pos_ref = match_cm_position(alignment_fa)['X14835.1_6927-7002'][0]
            a_to_c_dict[pos_ref[17]] = '17a'
            a_to_c_dict[pos_ref[21]] = '20a'
            a_to_c_dict[pos_ref[22]] = '20b'

            #get regular canonical positions
            #K01390.1_442-514 contains all canonical non subscripted positions
            pos_ref = match_cm_position(alignment_fa)['K01390.1_442-514'][0]
            variable_loop_start = 0
            variable_loop_end = 0
            for i, j in enumerate(pos_ref):
                a_to_c_dict[j] = i+1
                if i+1 == 45:
                    variable_loop_start = j+1
                if i+1 == 46:
                    variable_loop_end = j

            # annotate variable loop
            for j in range(variable_loop_start, variable_loop_end):
                a_to_c_dict[j] = 'e'
            return a_to_c_dict

        ref_dict = SeqIO.to_dict(SeqIO.parse(input.ref_fasta, "fasta"))
        #housekeeping_dict = SeqIO.to_dict(SeqIO.parse(input.ref_house_fasta, "fasta"))
        #ref_dict = {**ref_dict, **housekeeping_dict}

        alignment_pos_dict =  match_cm_position(input.alignment_fa)
        canonical_pos_dict = alignment_to_canonical_positions_mapper(input.alignment_fa)

        with open(input.cluster_map) as file:
            cluster_dict = yaml.safe_load(file)


        df = pd.read_csv(input.sam, sep = '\t')
        total_reads = len(df)

        groups = df.groupby("RNAME")
        data = []
        for name, g_df in groups:

            ref_seq = str(ref_dict[name].seq)
            ref_len = len(ref_seq)
            #print(ref_seq)

            g_df["mapped seq"] = g_df.apply(
                lambda row: get_matched_seq(row["CIGAR"], row["SEQ"])[0], axis=1
            )
            g_df["padded mapped seq"] = g_df.apply(
                lambda row: get_padded_seq(
                row["mapped seq"],
                row["POS"],
                0,
                0,
                row["SEQ"],
            ),
            axis=1,
            )
            #print(g_df.head())
            seqs = g_df["padded mapped seq"].to_list()
            seqs = [s + " " * (ref_len+1 - len(s)) for s in seqs]
            seqs = [list(s) for s in seqs]

            alignment_df = pd.DataFrame(seqs)
            nt_count_df = alignment_df.apply(pd.Series.value_counts, axis=0).fillna(0)
            nt_count_df_transposed = nt_count_df.transpose().copy()
            #print(alignment_df.head())
            #print(nt_count_df.head())
            nt_count_df_transposed.reset_index(inplace = True)
            nt_count_df_transposed.rename(columns = {'index':'ref_pos'}, inplace = True)
            nt_count_df_transposed['ref_nt'] = nt_count_df_transposed.apply(lambda row: ref_seq[int(row['ref_pos'])] if int(row['ref_pos']) < len(ref_seq) else 'N',
                                                                                            axis = 1)
            if name in alignment_pos_dict.keys():
                alignment_pos_ref = alignment_pos_dict[name][0]
                alignment_pos_ref += [alignment_pos_ref[-1]+1]
                canonical_pos_ref = [str(canonical_pos_dict[pos])
                                        if pos in canonical_pos_dict.keys()
                                        else '.'
                                        for pos in alignment_pos_ref]
                #print(nt_count_df_transposed.head(60))
                nt_count_df_transposed['align_pos'] = nt_count_df_transposed.apply(lambda row: alignment_pos_ref[int(row['ref_pos'])],
                                                                                    axis = 1)
                nt_count_df_transposed['canonical_pos'] = nt_count_df_transposed.apply(lambda row: canonical_pos_ref[int(row['ref_pos'])],
                                                                                    axis = 1)
            else:
                nt_count_df_transposed['align_pos']= 'Nan'
                nt_count_df_transposed['canonical_pos']= 'Nan'




            nt_count_df_transposed['RNAME'] = name

            anticodon = anticodonfamily_from_rname(name)
            nt_count_df_transposed['anticodon'] = anticodon

            cluster = 'Nan'
            if name in cluster_dict.keys():
                cluster = cluster_dict[name]
            else:
                print(f"{name} not found!")
                #sys.exit()
            nt_count_df_transposed['cluster'] = cluster


            nt_count_df_transposed['all reads'] = len(g_df)
            nt_count_df_transposed['RPM'] = nt_count_df_transposed.apply(lambda row: row['all reads']*1000000/total_reads, axis=1)

            nt_count_df_transposed['all nts'] = 0
            for nt in ['A', 'C', 'G', 'T', 'U', 'N', 'a', 't', 'u', 'g', 'c', 'n']:
                if nt in nt_count_df_transposed.columns:
                    nt_count_df_transposed['all nts'] = nt_count_df_transposed[nt] + nt_count_df_transposed['all nts']
            nt_count_df_transposed['RPM at position'] = nt_count_df_transposed.apply(lambda row: row['all nts']*1000000/total_reads, axis=1)
            nt_count_df_transposed['mismatch'] = 0
            for nt in ['a', 'c', 'g', 't', 'u', '-', 'n']:
                if nt in nt_count_df_transposed.columns:
                    nt_count_df_transposed['mismatch'] = nt_count_df_transposed[nt] + nt_count_df_transposed['mismatch']
            nt_count_df_transposed['match'] = 0
            for nt in ['A', 'C', 'G', 'T', 'U', 'N']:
                if nt in nt_count_df_transposed.columns:
                    nt_count_df_transposed['match'] = nt_count_df_transposed[nt] + nt_count_df_transposed['match']

            for nt in ['>', '<', 'C', 'T']:
                if nt not in nt_count_df_transposed.columns:
                    nt_count_df_transposed[nt] = 0


            nt_count_df_transposed['mismatch fraction'] = nt_count_df_transposed.apply(lambda row: row['mismatch']/(row['mismatch']+row['match']) if row['mismatch'] >0 else 0,
                                                  axis =1)
            nt_count_df_transposed['stop+mismatch fraction'] = nt_count_df_transposed.apply(lambda row: (row['>']+row['mismatch'])/(row['>']+row['mismatch']+row['match']) if row['mismatch'] >0 else 0,
                                                  axis =1)
            nt_count_df_transposed['read start fraction'] = nt_count_df_transposed.apply(lambda row: row['>']/row['all reads'],
                                                  axis =1)
            nt_count_df_transposed['read end fraction'] = nt_count_df_transposed.apply(lambda row: row['<']/row['all reads'],
                                                  axis =1)
            nt_count_df_transposed['m5C fraction'] = nt_count_df_transposed.apply(lambda row: row['C']/(row['C']+row['T']) if row['ref_nt'] in ['C', 'c'] and (row['C']+row['T'])>0 else 0,
                                                        axis =1)
            if '-' not in nt_count_df_transposed.columns:
                nt_count_df_transposed['-']=0
            if '>' not in nt_count_df_transposed.columns:
                nt_count_df_transposed['>']=0
            nt_count_df_transposed['stop fraction'] = nt_count_df_transposed.apply(
                lambda row: row['>']/(row['all nts']+row['-']+row['>']) if (row['all nts']+row['-']+row['>']) >0 else 0, axis=1)
            data.append(nt_count_df_transposed)


        df =  pd.concat(data)
        df['reads in sample'] = total_reads
        #df['Experiment'] = params.sample
        print(params.sample)
        #print(df.head(180))
        df.to_csv(output.tsv, sep = '\t', index=False)

rule per_cluster_nt_count:
    input:
        tsv =  'resources/coverage_counts/pre-filter_{reads_filter}/{ref_set}/clusters-ed-{e_cutoff}-mm-{m_cutoff}_{c_treatment}/{sample}.tsv'
    output:
        cluster_tsv =  'resources/coverage_counts/pre-filter_{reads_filter}/{ref_set}/clusters-ed-{e_cutoff}-mm-{m_cutoff}_{c_treatment}/{sample}_per_{group}.tsv'
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        import pandas as pd
        import numpy as np

        df = pd.read_csv(input.tsv, sep="\t")

        total_reads = df["reads in sample"].to_list()[0]

        df["T on refC"] = df.apply(
            lambda row: row["T"] if row["ref_nt"] == "C" else np.nan, axis=1
        )
        df["C on refC"] = df.apply(
            lambda row: row["C"] if row["ref_nt"] == "C" else np.nan, axis=1
        )

        data = []
        #group by cluster or anticodon
        groups = df.groupby(wildcards.group)
        for group, gdf in groups:

            # convert position descriptors to strings
            gdf = gdf.astype({"align_pos": str, "canonical_pos": str})

            # for each position within the cluster get refernce composition
            ref_count_dict = {}
            for p, x_df in gdf.groupby("align_pos"):
                ref_nts = x_df.groupby("ref_nt")["all nts"].sum().to_dict()
                ref_count_dict[p] = ref_nts

            # sum coverage counts for equivalent positions within the cluster
            p_df = gdf.groupby(["align_pos", "canonical_pos"]).sum()

            # reintroduce cluster or anticodon lable
            p_df[wildcards.group] = group

            p_df.reset_index(inplace=True)

            p_df.drop(columns="ref_pos", inplace=True)

            # recompute features that can not be computed from individual refs as sum
            p_df["mismatch fraction"] = p_df.apply(
                lambda row: row["mismatch"] / (row["mismatch"] + row["match"])
                if row["mismatch"] > 0
                else 0,
                axis=1,
            )
            p_df["stop+mismatch fraction"] = p_df.apply(
                lambda row: (row[">"] + row["mismatch"])
                / (row[">"] + row["mismatch"] + row["match"])
                if row["mismatch"] > 0
                else 0,
                axis=1,
            )
            p_df["read start fraction"] = p_df.apply(
                lambda row: row[">"] / row["all reads"], axis=1
            )
            p_df["read end fraction"] = p_df.apply(
                lambda row: row["<"] / row["all reads"], axis=1
            )
            p_df["C vs C+T"] = p_df.apply(
                lambda row: row["C"] / (row["C"] + row["T"])
                if (row["C"] + row["T"]) > 0
                else np.nan,
                axis=1,
            )

            for nt in list(set(gdf["ref_nt"].to_list())):
                p_df["ref_" + str(nt)] = p_df.apply(
                    lambda row: ref_count_dict[row["align_pos"]][nt]
                    if nt in ref_count_dict[row["align_pos"]].keys()
                    else 0,
                    axis=1,
                )

            p_df["ref C / (C+T)"] = p_df.apply(
                lambda row: row["ref_C"] / (row["ref_C"] + row["ref_T"])
                if (row["ref_C"] + row["ref_T"]) > 0
                else 0,
                axis=1,
            )
            p_df["ref C / all nts"] = p_df.apply(
                lambda row: row["ref_C"] / row["all nts"] if row["all nts"] > 0 else np.nan,
                axis=1,
            )
            p_df["m5C fraction T based"] = p_df.apply(
                lambda row: 1 - row["T on refC"] / (row["ref_C"])
                if (row["ref_C"]) > 0
                else np.nan,
                axis=1,
            )
            p_df["m5C fraction C based"] = p_df.apply(
                lambda row: row["C on refC"] / (row["ref_C"]) if (row["ref_C"]) > 0 else np.nan,
                axis=1,
            )
            p_df["m5C fraction"] = p_df.apply(
                lambda row: row["C on refC"] / (row["C on refC"] + row["T on refC"])
                if (row["C on refC"] + row["T on refC"]) > 0
                else np.nan,
                axis=1,
            )
            if '-' not in p_df.columns:
                p_df['-']=0
            if '>' not in p_df.columns:
                p_df['>']=0
            p_df['stop fraction'] = p_df.apply(
                lambda row: row['>']/(row['all nts']+row['-']+row['>']) if (row['all nts']+row['-']+row['>']) >0 else 0, axis=1)



            p_df["all reads"] = p_df["all reads"].max()
            p_df["RPM"] = p_df["RPM"].max()

            data.append(p_df)
        df = pd.concat(data)
        df["reads in sample"] = total_reads
        print(df.head(100))
        df.to_csv(output.cluster_tsv, sep="\t", index=False)
