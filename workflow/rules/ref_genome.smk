
# Author: Maria Waldl • code@waldl.org
# Version: 2024-01-24

import pandas as pd


def aggregate_fasta_files_for_all_accessions(wildcards):
    assembly_report_file = str(checkpoints.get_assembly_report_from_ncbi.get(assembly_id=config['assembly_id']).output)
    accessions = pd.read_csv(assembly_report_file, sep = '\t', comment='#', usecols = [6,9], names = ['accession', 'chromosom'], index_col=False)['accession'].to_list()
    genomic_fasta_files = [f"resources/references/genome/{accession}.fa" for accession in accessions]
    return genomic_fasta_files


rule download_genome_annotations_as_gff_gz:
    input:
    output: 'resources/references/genome_annotations/{assembly_id}_genomic.gff.gz'
    shell:
        'wget -O {output} https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/Danio_rerio/latest_assembly_versions/{wildcards.assembly_id}/{wildcards.assembly_id}_genomic.gff.gz'


rule fasta_files_from_ncbi_by_accession:
    input:
    output:
        fasta_file = "resources/references/genome/{accession}.fa"
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        import random
        import time
        from Bio import SeqIO
        from Bio import Entrez
        Entrez.email = "maria@tbi.univie.ac.at"
        with open(output.fasta_file, "w") as fasta_handle:
            time.sleep(random.random()*10)
            handle = Entrez.efetch(
                db="nucleotide",
                id=wildcards.accession,
                rettype="fasta",
                retmode="text",
            )
            sequence = SeqIO.read(handle, "fasta")
            r = SeqIO.write(sequence, fasta_handle, "fasta")
            if r != 1:
                print("Error while writing sequence:  " + sequence.id)



checkpoint get_assembly_report_from_ncbi:
# assembly report includes information for mapping chromosome ids to accession numbers
    input:
    output:
        fasta_file = "resources/references/genome/{assembly_id}_assembly_report.txt"
    #wildcard_constraints:
    #    assembly_id ="[GCF_000002035.6_GRCz11]+"
    shell:
        'wget -O {output} https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/Danio_rerio/latest_assembly_versions/{wildcards.assembly_id}/{wildcards.assembly_id}_assembly_report.txt'



rule get_mitochondrial_tRNAs_from_gff:
    input:
        gff_file = 'resources/references/genome_annotations/'+config['assembly_id']+'_genomic.gff',
        mitochondrial_fa = "resources/references/genome/"+config['mitochondrial_accession']+".fa"
    output:
        fa = 'resources/references/non-redundant/mitochondrial_gff.fa'
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        from BCBio import GFF
        from Bio import SeqIO

        mit_genome_record = SeqIO.read(input.mitochondrial_fa, "fasta")

        with open(output.fa, "w") as out:
            # limit_info = dict(gff_source=["RefSeq"], gff_type=["tRNA"])
            limit_info = dict(
                gff_id=[config["mitochondrial_accession"]], gff_type=["tRNA"]
            )
            in_handle = open(input.gff_file)
            for rec in GFF.parse(in_handle, limit_info=limit_info):
                for trna in rec.features:
                    out.write(
                        ">mt{}-{}(mt{}-{})\n".format(
                            trna.qualifiers["product"][0], trna.qualifiers["gene"][0],trna.qualifiers["product"][0], trna.qualifiers["gene"][0]
                        )
                    )
                    seq = mit_genome_record[
                        trna.location.start : trna.location.end
                    ].seq
                    if trna.location.strand == -1:
                        seq = seq.reverse_complement()
                    out.write(str(seq) + "CCA\n")
            in_handle.close()


rule get_mitochondrial_tRNAs_from_mt_tRNADB:
# removes structure information from mtRNADB download and constructs standard header format
    input:
        mitochondrial_fa = config['mitochondrial_tRNAs_from_mt_tRNADB']
    output:
        fa = 'resources/references/non-redundant/mitochondrial_tRNADB.fa'
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        from Bio import SeqIO
        with open(output.fa, "w") as out:
            in_handle = open(input.mitochondrial_fa)
            for line in in_handle:
                if line.startswith('>'):
                    id,species,taxid,aminoacid,anticodon = line.strip().lstrip('>').split('|')
                    out.write(
                        ">mt-{}-{}(mt-{}-{})\n".format(aminoacid, anticodon, aminoacid, anticodon)
                    )
                elif line.startswith('.') or line.startswith('('):
                    continue
                else:
                    out.write(line.strip().replace('U','T')+'CCA\n')

            in_handle.close()


rule get_high_confidence_genomic_tRNAs_from_tRNAscanSE:
    input:
        fasta = config['hc_genomic_tRNAs_from_GtRNAdb']
    output:
        non_redundant_fa = 'resources/references/non-redundant/genomic_hc.fa',
        non_redundant_BS_fa = 'resources/references/non-redundant/genomicBS_hc.fa',
        non_redundant_BSX_fa = 'resources/references/non-redundant/genomicBSX_hc.fa',
        bs_group_info = 'resources/references/non-redundant/group_info_hc.yaml',
        #bs_group_seqs = 'resources/references/non-redundant/group_info.txt'
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        import pandas as pd
        from Bio import SeqIO
        import yaml


        def get_merged_bs_ref(seqs):
            seqs = [s.replace("U", "T") for s in seqs]
            if len(seqs) == 1:
                return seqs[0]
            else:
                seq_len = len(seqs[0])
                ref = []
                for p in range(0, seq_len):
                    nts = [s[p] for s in seqs]
                    if "C" in nts:
                        ref.append("C")
                    else:
                        if len(set(nts)) != 1:
                            print("error: diverse sequences")
                        ref.append(nts[0])
                return "".join(ref)


        def get_merged_bs_refx(seqs):
            seqs = [s.replace("U", "T") for s in seqs]
            if len(seqs) == 1:
                return seqs[0]
            else:
                seq_len = len(seqs[0])
                ref = []
                for p in range(0, seq_len):
                    nts = [s[p] for s in seqs]
                    if len(set(nts)) != 1:
                        ref.append("X")
                    else:
                        ref.append(nts[0])
            return "".join(ref)


        # read sequences from GtRNAdb file
        seqs = []
        for record in SeqIO.parse(input.fasta, "fasta"):
            seq = str(record.seq)
            info = {}
            data = record.id.split("-")
            info["name"] = "-".join(data[1:5])+'('+'-'.join(data[1:3])+')'
            info["aminoacid"] = data[1]
            info["anticodon"] = data[2]
            info["seqgroup"] = int(data[3])
            info["gene_nr"] = int(data[4])
            info["seq"] = seq + 'CCA'
            info["BSseq"] = seq.replace("C", "U")
            seqs.append(info)
        df = pd.DataFrame.from_records(seqs)
        print("{} tRNAs in total".format(len(df)))

        # check if one sequence is part of multiple anticodon or seqeunce groups
        groups = df.groupby("seq")
        for s, gdf in groups:
            if len(list(set(gdf["anticodon"].to_list()))) > 1:
                print(s, "not unique in anticodon")
            if len(list(set(gdf["seqgroup"].to_list()))) > 1:
                print(s, "not unique in seqgroup")

        # only keep one version if same sequence occures in multiple genomic regions
        df.sort_values(["aminoacid", "anticodon", "seqgroup", "gene_nr"], inplace=True)
        df.drop_duplicates(
            subset=["aminoacid", "anticodon", "seqgroup"],
            keep="first",
            inplace=True,
            ignore_index=True,
        )

        # check if any sequence is not unique
        groups = df.groupby("seq")
        for s, gdf in groups:
            if len(gdf) > 1:
                print("non unique refernce sequences:")
                print(gdf.head(), len(gdf))
        print("{} unique reference seqeunces".format(len(df)))

        # write unique sequences
        with open(output.non_redundant_fa, "w") as nr_genomic_fa:
            groups = df.groupby("seq")
            for seq, gdf in groups:
                nr_genomic_fa.write(">" + gdf["name"].to_list()[0] + "\n")
                nr_genomic_fa.write(seq.replace("U", "T") + "\n")
                if len(gdf) != 1:
                    print(gdf["name"].to_list(), "not unique by seq")


        # write unique after BS teatment
        group_info = {}
        with open(output.non_redundant_BSX_fa, "w") as bxnr_genomic_fa:
            with open(output.non_redundant_BS_fa, "w") as bnr_genomic_fa:
                groups = df.groupby(by="BSseq")
                print(len(groups))
                for n, g in groups:
                    g.sort_values(
                        ["aminoacid", "anticodon", "seqgroup", "gene_nr"], inplace=True
                    )
                    if len(g) == 1:
                        bnr_genomic_fa.write(">" + g["name"].to_list()[0] + "\n")
                        bxnr_genomic_fa.write(">" + g["name"].to_list()[0] + "\n")
                        group_info[g["name"].to_list()[0]] = list(
                            zip(g["name"].to_list(), g["seq"].to_list())
                        )
                    else:
                        bnr_genomic_fa.write(">" + g["name"].to_list()[0] + "[B]\n")
                        bxnr_genomic_fa.write(">" + g["name"].to_list()[0] + "[B]\n")
                        group_info[g["name"].to_list()[0] + "[B]"] = [
                            ["merged", get_merged_bs_refx(g["seq"].to_list())]
                        ] + list(zip(g["name"].to_list(), g["seq"].to_list()))
                    bnr_genomic_fa.write(get_merged_bs_ref(g["seq"].to_list()) + "\n")
                    bxnr_genomic_fa.write(
                        get_merged_bs_refx(g["seq"].to_list()) + "\n"
                    )

        with open(output.bs_group_info, "w") as file:
            yaml.dump(group_info, file)


rule get_all_genomic_tRNAs_from_tRNAscanSE:
    input:
        scan_tsv = config['downloaded_tRNAscanSE_summary'],
        tRNA_name_map = config['downloaded_tRNAscanSE_name_mapping'],
        chromosom_accession_map = "resources/references/genome/"+config['assembly_id']+"_assembly_report.txt",
        genomic_fasta_files = aggregate_fasta_files_for_all_accessions
    output:
        non_redundant_fa = 'resources/references/non-redundant/genomic_fullScan.fa',
        non_redundant_fa_highconfidence = 'resources/references/non-redundant/genomic_fullScan_highconfidence.fa',
        non_redundant_fa_additional = 'resources/references/non-redundant/genomic_fullScan_additional.fa',
        non_redundant_BS_fa = 'resources/references/non-redundant/genomicBS_fullScanc.fa',
        non_redundant_BSX_fa = 'resources/references/non-redundant/genomicBSX_fullScan.fa',
        bs_group_info = 'resources/references/non-redundant/group_info_fullScan.yaml',
        #bs_group_seqs = 'resources/references/non-redundant/group_info.txt'
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        import pandas as pd
        from Bio import SeqIO
        import yaml

        def get_merged_bs_ref(seqs):
            seqs = [s.replace("U", "T") for s in seqs]
            if len(seqs) == 1:
                return seqs[0]
            else:
                seq_len = len(seqs[0])
                ref = []
                for p in range(0, seq_len):
                    nts = [s[p] for s in seqs]
                    if "C" in nts:
                        ref.append("C")
                    else:
                        if len(set(nts)) != 1:
                            print("error: diverse sequences")
                        ref.append(nts[0])
                return "".join(ref)


        def get_merged_bs_refx(seqs):
            seqs = [s.replace("U", "T") for s in seqs]
            if len(seqs) == 1:
                return seqs[0]
            else:
                seq_len = len(seqs[0])
                ref = []
                for p in range(0, seq_len):
                    nts = [s[p] for s in seqs]
                    if len(set(nts)) != 1:
                        ref.append("X")
                    else:
                        ref.append(nts[0])
            return "".join(ref)

        def get_seq(accession, start, end, intron_start, intron_end):
            fasta = "resources/references/genome/"+accession+".fa"
            record = SeqIO.read(fasta, "fasta")
            seq = record.seq[min(start,end)-1:max(start,end)]
            if intron_start>0 or intron_end>0:
                intron_local_start = min(intron_start,intron_end)-min(start,end)
                intron_local_end = max(intron_start,intron_end)-min(start,end)
                seq = seq[0:intron_local_start]+seq[intron_local_end+1:]
            if start>end:
                seq = seq.reverse_complement()
            return str(seq).replace('U','T')+'CCA'

        # get chromosome to accession multimapper
        chrom_accession_mapper = pd.read_csv(input.chromosom_accession_map, sep = '\t', comment='#', usecols = [6,9], names = ['accession', 'chromosom'], index_col=False).set_index('chromosom').to_dict(orient='index')

        # get tRNAscan-SE_id to GtRNAdb_id mapper
        id_mapper = pd.read_csv(input.tRNA_name_map, sep = '\t', header = 0).set_index('tRNAscan-SE_id') .to_dict(orient="index")

        # read info tsv
        columns = ['chromosom',
        'tRNA#',
        'start',
        'end',
        'type', #aminoacide
        'codon', #anticodon
        'intron start',
        'intron end',
        'HMM score',
        'str score',
        'hit score',
        'isotype origin',
        'isotype cm',
        'typ score',
        'sources',
        'note']
        df = pd.read_csv(input.scan_tsv, sep = '\t', skiprows=3,header = 0, names = columns, index_col=False, )#nrows = 100)

        df['tRNAscan-SE_id'] = df.apply(lambda row: row['chromosom']+'.trna'+str(row['tRNA#']), axis = 1)
        df['name_full'] = df.apply(lambda row: id_mapper[row['tRNAscan-SE_id']]['GtRNAdb_id'], axis = 1)
        df['name'] = df.apply(lambda row: "-".join(row['name_full'].split('-')[1:4])+'('+'-'.join([row['type'],row['codon']])+')', axis =1)
        df['accession'] = df.apply(lambda row: chrom_accession_mapper[row['chromosom']]['accession'], axis = 1)
        df['aminoacid'] = df.apply(lambda row: row['name_full'].split('-')[1], axis =1)
        df['anticodon'] = df.apply(lambda row: row['name_full'].split('-')[2], axis =1)
        df['seqgroup'] = df.apply(lambda row: row['name_full'].split('-')[3], axis =1)
        df['gene_nr'] = df.apply(lambda row: row['name_full'].split('-')[4], axis =1)
        df['seq'] = df.apply(lambda row:get_seq(row['accession'], int(row['start']), int(row['end']), int(row['intron start']), int(row['intron end'])), axis =1)
        df["BSseq"] = df.apply(lambda row: row['seq'].replace("C", "T"),axis =1)
        df['high confidence'] = df.apply(lambda row: 1 if 'high confidence set' in row['note'] else 0 , axis = 1)
        df['not high confidence'] = df.apply(lambda row: 0 if 'high confidence set' in row['note'] else 1 , axis = 1)

        print("{} tRNAs in total".format(len(df)))
        # check if one sequence is part of multiple anticodon or seqeunce groups
        groups = df.groupby("seq")
        for s, gdf in groups:
            if len(list(set(gdf["anticodon"].to_list()))) > 1:
                print(s, "not unique in anticodon")
            if len(list(set(gdf["seqgroup"].to_list()))) > 1:
                print(s, "not unique in seqgroup")

        # only keep one version if same sequence occures in multiple genomic regions
        df.sort_values(["aminoacid", "anticodon", "seqgroup", 'not high confidence',"gene_nr"], inplace=True)
        df.drop_duplicates(
            subset=["aminoacid", "anticodon", "seqgroup"],
            keep="first",
            inplace=True,
            ignore_index=True,)
        print("{} tRNAs in total".format(len(df)))


        # check if any sequence is not unique
        groups = df.groupby("seq")
        for s, gdf in groups:
            if len(gdf) > 1:
                print("non unique refernce sequences:")
                print(gdf.head(), len(gdf))
        print("{} unique reference sequences".format(len(df)))

        # write unique sequences
        with open(output.non_redundant_fa, "w") as nr_genomic_fa:
            groups = df.groupby("seq")
            for seq, gdf in groups:
                hc = ''
                if 1 in gdf['high confidence'].to_list():
                    hc = '[H]'
                nr_genomic_fa.write(">" + gdf["name"].to_list()[0]+hc + "\n")
                nr_genomic_fa.write(seq.replace("U", "T") + "\n")
                if len(gdf) != 1:
                    print(gdf["name"].to_list(), "not unique by seq")

        # write unique high confidence sequences
        with open(output.non_redundant_fa_highconfidence , "w") as nr_genomic_fa:
            groups = df.groupby("seq")
            for seq, gdf in groups:
                hc = ''
                if 1 in gdf['high confidence'].to_list():
                    hc = '[H]'
                else:
                    continue
                nr_genomic_fa.write(">" + gdf["name"].to_list()[0]+hc + "\n")
                nr_genomic_fa.write(seq.replace("U", "T") + "\n")
                if len(gdf) != 1:
                    print(gdf["name"].to_list(), "not unique by seq")


        # write unique non  high confidence sequences
        with open(output.non_redundant_fa_additional, "w") as nr_genomic_fa:
            groups = df.groupby("seq")
            for seq, gdf in groups:
                hc = ''
                if 1 in gdf['high confidence'].to_list():
                    continue
                nr_genomic_fa.write(">" + gdf["name"].to_list()[0]+hc + "\n")
                nr_genomic_fa.write(seq.replace("U", "T") + "\n")
                if len(gdf) != 1:
                    print(gdf["name"].to_list(), "not unique by seq")


        # write unique after BS teatment
        group_info = {}
        with open(output.non_redundant_BSX_fa, "w") as bxnr_genomic_fa:
            with open(output.non_redundant_BS_fa, "w") as bnr_genomic_fa:
                groups = df.groupby(by="BSseq")
                print(len(groups))
                for n, g in groups:
                    g.sort_values(
                        ["aminoacid", "anticodon", "seqgroup", "gene_nr"], inplace=True
                    )
                    if len(g) == 1:
                        bnr_genomic_fa.write(">" + g["name"].to_list()[0] + "\n")
                        bxnr_genomic_fa.write(">" + g["name"].to_list()[0] + "\n")
                        group_info[g["name"].to_list()[0]] = list(
                            zip(g["name"].to_list(), g["seq"].to_list())
                        )
                    else:
                        bnr_genomic_fa.write(">" + g["name"].to_list()[0] + "[B]\n")
                        bxnr_genomic_fa.write(">" + g["name"].to_list()[0] + "[B]\n")
                        group_info[g["name"].to_list()[0] + "[B]"] = [
                            ["merged", get_merged_bs_refx(g["seq"].to_list())]
                        ] + list(zip(g["name"].to_list(), g["seq"].to_list()))
                    bnr_genomic_fa.write(get_merged_bs_ref(g["seq"].to_list()) + "\n")
                    bxnr_genomic_fa.write(
                        get_merged_bs_refx(g["seq"].to_list()) + "\n"
                    )

        with open(output.bs_group_info, "w") as file:
            yaml.dump(group_info, file)



rule get_merged_tRNA_refs:
    input:
        genomic = config['genomic_refs'],
        mitochondrial = config['mitochondrial_refs']
    output:
        'resources/references/all_refs.fa'
    shell:
        "cat {input.genomic} {input.mitochondrial} > {output}"


rule get_selected_refs:
    input:
        min_cov_ids = 'resources/min_coverage_refs/pre-filter_'+config['reads_filter']+'/min_cov_refs.yaml',
        tRNAs_fa = 'resources/references/all_refs.fa'
    output:
        fasta = 'resources/references/selected_refs.fa',
        fastaN = 'resources/references/selected_refs_N.fa',
        tRNAs_fa = 'resources/references/all_refs_N.fa'
    run:
        # Corrected version: 2025-08-01
        import yaml
        from Bio import SeqIO
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq

        with open(input.min_cov_ids) as file:
            refs = yaml.safe_load(file)

        all_seqs_N = []
        filtered_sequences = []
        filtered_sequences_N = []
        for record in SeqIO.parse(input.tRNAs_fa, "fasta"):
            new_record_all = SeqRecord(
                Seq(str(record.seq) + 'N'),
                id=record.id,
                description=record.description
            )
            all_seqs_N.append(new_record_all)
            if record.id in refs:
                filtered_sequences.append(record)
                new_record = SeqRecord(
                    Seq(str(record.seq) + 'N'),
                    id=record.id,
                    description=record.description
                )
                filtered_sequences_N.append(new_record)

        SeqIO.write(all_seqs_N, output.tRNAs_fa, "fasta") # all + N
        SeqIO.write(filtered_sequences, output.fasta, "fasta")
        SeqIO.write(filtered_sequences_N, output.fastaN, "fasta")

rule remove_CCA_from_manual_refs:
    input:
        fasta = config['raw_manual_refs']
    output:
        fasta = 'resources/references/non-redundant/manual-noCCA.fa',
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        from Bio import SeqIO
        from Bio.SeqRecord import SeqRecord
        with open(output.fasta, 'w') as output_fasta:
            for record in SeqIO.parse(input.fasta, "fasta"):
                seq = str(record.seq[0:-3]).upper().replace('U', 'T')
                output_fasta.write('>' + record.id + '\n')
                output_fasta.write(seq + '\n')


rule copy_manual_refs:
    input:
        fasta = config['raw_manual_refs']
    output:
        fasta = 'resources/references/manual_refs.fa',
    run:
        # Author: Maria Waldl • code@waldl.org
        # Version: 2024-01-24
        from Bio import SeqIO
        from Bio.SeqRecord import SeqRecord
        with open(output.fasta, 'w') as output_fasta:
            for record in SeqIO.parse(input.fasta, "fasta"):
                seq = str(record.seq).upper().replace('U', 'T')
                output_fasta.write('>' + record.id + '\n')
                output_fasta.write(seq + '\n')
