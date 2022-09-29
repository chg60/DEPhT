from depht.functions.fasta import write_fasta
from depht.functions.find_homologs import __feature_in_prophage
from depht.functions.subprocess import run_command
from pandas import read_csv


def find_defense_systems(contigs, prophage_coords, tmp_dir, cpus,
                         cache_scores=True):
    for contig, contig_prophage_coords in zip(contigs, prophage_coords):
        map_geneid_to_feature = dict()
        contig_geneids, contig_sequences = list(), list()
        contig_features = list()
        for i, feature in enumerate(contig.genes):
            geneid = contig.gene_ids[i]

            translation = feature.qualifiers["translation"][0]

            if __feature_in_prophage(feature, contig_prophage_coords):
                map_geneid_to_feature[geneid] = feature
                contig_geneids.append(geneid)
                contig_sequences.append(translation)
                contig_features.append(feature)

        if contig_sequences:
            defense_outdir = tmp_dir.joinpath(f"{contig.id}")
            defense_outdir.mkdir()

            contig_fasta = tmp_dir.joinpath(f"{contig.id}.faa")
            write_fasta(contig_geneids, contig_sequences, contig_fasta)

            contig_gff = tmp_dir.joinpath(f"{contig.id}.gff")
            write_gff(contig.id, contig_geneids, contig_features, contig_gff)

            defensefinder(contig_fasta, defense_outdir, cpus)
            defensefinder_products = parse_defensefinder(defense_outdir)

            padloc(contig_fasta, contig_gff, defense_outdir, cpus)
            padloc_products = parse_padloc(contig.id, defense_outdir)

            hhsearch_scores = [0.0] * len(contig.gene_ids)
            for geneid, product in padloc_products.items():
                feature = map_geneid_to_feature[geneid]
                feature.qualifiers["product"] = [product]

                hhsearch_scores[contig.gene_ids.index(geneid)] = float(100)

            for geneid, product in defensefinder_products.items():
                feature = map_geneid_to_feature[geneid]
                feature.qualifiers["product"] = [product]

                hhsearch_scores[contig.gene_ids.index(geneid)] = float(100)

        else:
            hhsearch_scores = [0.0] * len(contig.gene_ids)

        if cache_scores:
            contig.update_hhsearch_scores(hhsearch_scores)


def parse_defensefinder(outdir):
    defensefinder_products = dict()
    systems_tsv = outdir.joinpath("defense_finder_systems.tsv")
    if not systems_tsv.is_file():
        return defensefinder_products

    systems_dataframe = read_csv(systems_tsv, sep="\t")

    for i in range(len(systems_dataframe["sys_id"])):
        gene_ids = systems_dataframe["protein_in_syst"][i].split(",")
        products = systems_dataframe["name_of_profiles_in_sys"][i].split(",")

        for gene_id, product in zip(gene_ids, products):
            defensefinder_products[gene_id] = product

    return defensefinder_products


def defensefinder(fasta_file, output_dir, cpus):
    command = ("defense-finder run "
               f"--out-dir {output_dir} "
               f"--workers {cpus} "
               f"{fasta_file} ")

    run_command(command)


def parse_padloc(contig_name, outdir):
    padloc_products = dict()
    systems_csv = outdir.joinpath(f"{contig_name}_padloc.csv")
    if not systems_csv.is_file():
        return padloc_products

    systems_dataframe = read_csv(systems_csv)

    for i in range(len(systems_dataframe["seqid"])):
        gene_id = systems_dataframe["target.name"][i]
        product = systems_dataframe["protein.name"][i]

        padloc_products[gene_id] = product

    return padloc_products


def padloc(fasta_file, gff_file, output_dir, cpus):
    command = (f"padloc --faa {fasta_file} --gff {gff_file} "
               f"--outdir {output_dir} "
               f"--cpu {cpus} --fix-prodigal")

    run_command(command)


def write_gff(contig_name, geneids, features, filepath):
    gff_writer = open(filepath, "w")

    for i in range(len(geneids)):
        geneid = geneids[i]
        feature = features[i]

        gene_num = geneid.split("_")[-1]

        gff_writer.write(f"{contig_name}\t")
        gff_writer.write("Prodigal\t")
        gff_writer.write("CDS\t")
        gff_writer.write(f"{feature.location.start}\t")
        gff_writer.write(f"{feature.location.start}\t")
        gff_writer.write(".\t")
        gff_writer.write(f"{feature.location.strand}\t")
        gff_writer.write("0\t")
        gff_writer.write(f"ID=1_{gene_num};\n")

    # Close the file handle
    gff_writer.close()
