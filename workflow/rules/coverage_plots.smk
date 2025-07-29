
# Author: Maria Waldl • code@waldl.org
# Version: 2024-01-24

rule coverage_plots_per_group_sample:
    input:
        mismatch_tsvs = expand('resources/coverage_counts/pre-filter_{reads_filter}/{ref_set}/clusters-ed-{e_cutoff}-mm-{m_cutoff}_{c_treatment}/{sample}_per_cluster.tsv', reads_filter = '{reads_filter}',ref_set='{ref_set}',sample = SAMPLES, e_cutoff='{e_cutoff}', m_cutoff='{m_cutoff}', c_treatment='{c_treatment}' ),
        cluster_names = 'resources/cluster/pre-filter_{reads_filter}/{ref_set}/anticodonbased_clusternames-ed-{e_cutoff}-mm-{m_cutoff}_{c_treatment}.yaml'
    output:
        pdf = 'results/coverage_plots/pre-filter_{reads_filter}/{ref_set}/per_cluster-ed-{e_cutoff}-mm-{m_cutoff}_{c_treatment}/coverage_plots_per_sample_{treatment}/color_scheme_{coloring}/summary.pdf'
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24

        import pandas as pd
        import matplotlib.pyplot as plt
        import os
        import subprocess
        import yaml

        columnid_color = [
            ["A", "#1f77b4", ""],
            ["C", "#ff7f0e", ""],
            ["G", "#2ca02c", ""],
            ["T", "#d62728", ""],
            ["a", "#aec7e8", ""],
            ["c", "#ffbb78", ""],
            ["g", "#98df8a", ""],
            ["t", "#ff9896", ""],
            ["-", "grey"   , ""],
            [">", "white"  , ""],
            ["<", "white"  , ""],
        ]
        nts_to_plot = ["A", "C", "G", "T", "a", "c", "g", "t", "-"]

        if wildcards.coloring == "mismatch":
            columnid_color = [
                ["a", "#1f77b4", ""],
                ["c", "#ff7f0e", ""],
                ["g", "#2ca02c", ""],
                ["t", "#d62728", ""],
                ["-", "grey"   , ""],
                ["A", "#aec7e8", ""],
                ["C", "#ffbb78", ""],
                ["G", "#98df8a", ""],
                ["T", "#ff9896", ""],
            ]
            nts_to_plot = ["a", "c", "g", "t", "A", "C", "G", "T", "-"]
        elif wildcards.coloring == "simplified":
            columnid_color = [
                ["A", "#1f77b4", ""],
                ["C", "#ff7f0e", ""],
                ["G", "#2ca02c", ""],
                ["T", "#d62728", ""],
                ["a", "#1f77b4", ""],
                ["c", "#ff7f0e", ""],
                ["g", "#2ca02c", ""],
                ["t", "#d62728", ""],
                ["-", "grey"   , ""],
            ]
            nts_to_plot = ["a", "A", "c", "C", "g", "G", "t", "T", "-"]
        elif wildcards.coloring == "bisulfite":
            columnid_color = [
                ["C", "black", "xx"],
                ["a", "#1f77b4", "."],
                ["c", "#ff7f0e", "."],
                ["g", "#2ca02c", "."],
                ["t", "#d62728", "."],
                ["-", "grey", "."],
                ["A", "#aec7e8", "."],
                ["G", "#98df8a", "."],
                ["T", "#ff9896", "."],
            ]
            nts_to_plot = ["C", "a", "c", "g", "t", "A", "G", "T", "-"]
        elif wildcards.coloring == "stop":
            columnid_color = [
                [">", "grey"     , "."],
                ["<", "darkgrey" , "."],
            ]
            nts_to_plot = [">", "<"]

        ref_columnid_color = [
            ["ref_A", "#1f77b4"],
            ["ref_C", "#ff7f0e"],
            ["ref_G", "#2ca02c"],
            ["ref_T", "#d62728"],
            ["ref_U", "#d62728"],
            ["ref_a", "#1f77b4"],
            ["ref_c", "#ff7f0e"],
            ["ref_g", "#2ca02c"],
            ["ref_t", "#d62728"],
            ["ref_u", "#d62728"],
        ]
        ref_nts_to_plot = [
            "ref_A",
            "ref_C",
            "ref_G",
            "ref_T",
            "ref_U",
            "ref_a",
            "ref_c",
            "ref_g",
            "ref_t",
            "ref_u",
        ]


        data = []
        for tsv_file in input.mismatch_tsvs:
            sample = tsv_file.split("/")[-1].replace("_per_cluster.tsv", "")
            #sample = tsv_file.split("/")[-1].replace(".tsv", "").split("_")[0]
            if sample_dict[sample]["treatment"] != wildcards.treatment:
                continue
            df = pd.read_csv(tsv_file, sep="\t")
            df["sample"] = sample
            df["timepoint"] = sample_dict[sample]["timepoint"]
            df["time point name"] = sample_dict[sample]["timepoint_name"]
            df["treatment"] = sample_dict[sample]["treatment"]
            df["position"] = df.apply(
                lambda row: str(row["align_pos"]).zfill(3)
                + " ("
                + str(row["canonical_pos"])
                + ")",
                axis=1,
            )
            data.append(df)
        df = pd.concat(data)




        with open(input.cluster_names) as file:
            cluster_name_dict = yaml.safe_load(file)

        pdfs = []

        plot_dir = "/".join(output.pdf.split("/")[0:-1])
        for group, g_df in df.groupby("cluster"):
            g_df["cluster"] = group
            fig = plt.figure(figsize=(36, 16), constrained_layout=True)
            cluster_name = group

            if group != "Nan" and group in cluster_name_dict.keys():
                cluster_name = cluster_name_dict[group]
            plt.suptitle("cluster: " + str(cluster_name) + "\ntreatment: " + wildcards.treatment)

            if wildcards.treatment == "DM":
                subfigures = fig.subfigures(3, 3, hspace=0.0)
            elif wildcards.treatment == "MOCK":
                subfigures = fig.subfigures(3, 3, hspace=0.0)
            fig_path = os.path.join(plot_dir, str(cluster_name) + ".pdf")
            pdfs.append(fig_path)

            for timepoint, t_df in g_df.groupby("timepoint"):
                for sample, sample_df in t_df.groupby("sample"):
                    replicate = sample_dict[sample]['replicate']
                    axs = subfigures[timepoint, replicate-1].subplots(
                    #axs = subfigures[timepoint - 1, replicate-1].subplots(
                        nrows=2,
                        ncols=1,
                        gridspec_kw={"height_ratios": [8, 1], "hspace": 0.0, "wspace": 0},
                        sharex=True,
                    )
                    #columnid_color = [
                        #["a", "#1f77b4", ""],
                        #["c", "#ff7f0e", ""],
                        #["g", "#2ca02c", ""],
                        #["t", "#d62728", ""],
                        #["-", "grey"   , ""],
                        #["A", "#aec7e8", ""],
                        #["C", "#ffbb78", ""],
                        #["G", "#98df8a", ""],
                        #["T", "#ff9896", ""],
                    #]
                    #nts_to_plot = ["a", "c", "g", "t", "A", "C", "G", "T", "-"]
                    y = [
                        m[0]
                        for m in columnid_color
                        if m[0] in sample_df.columns and m[0] in nts_to_plot
                    ]
                    color = [
                        m[1]
                        for m in columnid_color
                        if m[0] in sample_df.columns and m[0] in nts_to_plot
                    ]
                    hatch = [
                        m[2]
                        for m in columnid_color
                        if m[0] in sample_df.columns and m[0] in nts_to_plot
                    ]
                    #Order by ref_pos
                    sample_df.sort_values(by="position", inplace=True, ignore_index=True)
                    sample_df.plot.bar(x="position", y=y, color=color,  stacked=True, ax=axs[0])
                    axs[0].set_ylabel("read count")
                    axs[0].set_title(
                        sample_dict[sample]["timepoint_name"] + " (" + sample + ")"
                    )
                    axs[0].set_ylim(0, sample_df["all reads"].max())

                    for nt in ref_nts_to_plot:
                        if nt in sample_df.columns:
                            sample_df[nt] = sample_df.apply(
                                lambda row: row[nt] / row["all nts"] if row[nt] > 0 else 0,
                                axis=1,
                            )
                    #ref_columnid_color = [
                        #["ref_A", "#1f77b4"],
                        #["ref_C", "#ff7f0e"],
                        #["ref_G", "#2ca02c"],
                        #["ref_T", "#d62728"],
                        #["ref_U", "#d62728"],
                        #["ref_a", "#1f77b4"],
                        #["ref_c", "#ff7f0e"],
                        #["ref_g", "#2ca02c"],
                        #["ref_t", "#d62728"],
                        #["ref_u", "#d62728"],
                    #]
                    y = [
                        m[0]
                        for m in ref_columnid_color
                        if m[0] in sample_df.columns and m[0] in ref_nts_to_plot
                    ]
                    color = [
                        m[1]
                        for m in ref_columnid_color
                        if m[0] in sample_df.columns and m[0] in ref_nts_to_plot
                    ]
                    bars = sample_df.plot.bar(
                        x="position", y=y, color=color, stacked=True, ax=axs[1], legend=None
                    )

                    axs[1].set_ylim(0, 1)
                    axs[1].set_ylabel("references\ncomposition")
                    axs[1].set_xlabel("alignment position (canonical position)")
                    if wildcards.treatment == "DM":
                        if replicate-1 == 2:
                            timepoint += 1
                    elif wildcards.treatment == "MOCK":
                        if replicate-1 == 2:
                            timepoint += 1

            fig.savefig(fig_path)
            plt.close("all")


        pdfs.sort()
        call_args = ["pdfunite"] + pdfs + [output.pdf]
        results = subprocess.run(call_args, capture_output=True)



rule get_all_coverage_plots:
    input:
        'results/coverage_plots/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/per_cluster-ed-3-mm-50_DM/coverage_plots_per_sample_DM/color_scheme_mismatch/summary.pdf',
        'results/coverage_plots/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/per_cluster-ed-3-mm-50_DM/coverage_plots_per_sample_MOCK/color_scheme_mismatch/summary.pdf',
        #'results/coverage_plots/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/per_cluster-ed-3-mm-50_DM/coverage_plots_per_sample_BS/color_scheme_mismatch/summary.pdf',
        #'results/coverage_plots/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/per_cluster-ed-3-mm-50_DM/coverage_plots_per_sample_DM/color_scheme_stop/summary.pdf',
        #'results/coverage_plots/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/per_cluster-ed-3-mm-50_DM/coverage_plots_per_sample_BS/color_scheme_stop/summary.pdf',
        #'results/coverage_plots/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/per_cluster-ed-3-mm-50_DM/coverage_plots_per_sample_BS/color_scheme_all/summary.pdf',
        #'results/coverage_plots/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/per_cluster-ed-3-mm-50_DM/coverage_plots_per_sample_BS/color_scheme_bisulfite/summary.pdf',
        #'results/coverage_plots/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/per_cluster-ed-2-mm-50_DM/coverage_plots_per_sample_DM/color_scheme_mismatch/summary.pdf',
        #'results/coverage_plots/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/per_cluster-ed-2-mm-50_DM/coverage_plots_per_sample_BS/color_scheme_mismatch/summary.pdf',
        #'results/coverage_plots/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/per_cluster-ed-2-mm-50_DM/coverage_plots_per_sample_DM/color_scheme_stop/summary.pdf',
        #'results/coverage_plots/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/per_cluster-ed-2-mm-50_DM/coverage_plots_per_sample_BS/color_scheme_stop/summary.pdf',
        #'results/coverage_plots/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/per_cluster-ed-2-mm-50_DM/coverage_plots_per_sample_BS/color_scheme_all/summary.pdf',
        #'results/coverage_plots/pre-filter_'+config['reads_filter']+'/'+config['ref_set']+'/per_cluster-ed-2-mm-50_DM/coverage_plots_per_sample_BS/color_scheme_bisulfite/summary.pdf',
    output:
        'results/coverage_plots/done.txt'
    run:
        with open(str(output),'w') as file:
            for f in input:
                file.write(f+'\n')
