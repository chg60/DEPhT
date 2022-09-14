from depht.functions.find_homologs import __feature_in_prophage
from depht.functions.fasta import write_fasta
from depht.functions.subprocess import run_command
from pandas import read_csv


def find_defense_systems(contigs, prophage_coords, tmp_dir, cpus,
                         cache_scores=True):
    for contig, contig_prophage_coords in zip(contigs, prophage_coords):
        map_geneid_to_feature = dict()
        contig_geneids, contig_sequences = list(), list()
        for i, feature in enumerate(contig.genes):
            geneid = contig.gene_ids[i]

            translation = feature.qualifiers["translation"][0]

            if __feature_in_prophage(feature, contig_prophage_coords):
                map_geneid_to_feature[geneid] = feature
                contig_geneids.append(geneid)
                contig_sequences.append(translation)

        if contig_sequences:
            contig_fasta = tmp_dir.joinpath(f"{contig.id}.faa")
            defensefinder_outdir = tmp_dir.joinpath(f"{contig.id}")
            defensefinder_outdir.mkdir()
            write_fasta(contig_geneids, contig_sequences, contig_fasta)

            defensefinder(contig_fasta, defensefinder_outdir, cpus)

            defensefinder_products = parse_defensefinder(defensefinder_outdir)

            hhsearch_scores = [0.0] * len(contig.gene_ids)
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
