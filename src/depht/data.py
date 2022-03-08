import datetime

DATE = datetime.date.today().strftime("%d-%b-%Y").upper()

GLOBAL_VARIABLES = {
    "sequences": {
        "name": "genes",
        "default_product": "hypothetical protein",
        "annotations": {
            "molecule_type": "DNA",
            "topology": "linear",
            "data_file_division": "PHG",
            "date": DATE,
            "accessions": [],
            "sequence_version": "1",
            "keywords": [],
            "source": "",
            "organism": "",
            "taxonomy": [],
            "comment": [""]},
        "feature_types": ("CDS", "tRNA", "tmRNA")},
    "bacterial_sequences": {
        "name": "bacterial_genes",
        "dir_name": "bacteria"},
    "phage_sequences": {
        "name": "phage_genes",
        "prophage_prefix": "prophi",
        "prophage_delimiter": "-",
        "dir_name": "phage"},
    "phage_homologs": {
        "essential_name": "essential",
        "extended_name": "extended",
        "dir_name": "hhsearch"},
    "shell_db": {
        "fasta_name": "bacterial_genes.fasta",
        "hex_value_name": "bacterial_genes.pbv",
        "dir_name": "mmseqs"},
    "reference_db": {
        "name": "references",
        "dir_name": "blastn"},
    "classifier": {
        "name": "classifier.pkl"},
    "model_storage": {
        "tmp_dir": "tmp",
        "model_dir": "models",
        "home_dir": ".depht"}}

PARAMETERS = {
    "annotation": {
        "min_length": 20000,
        "meta_length": 100000},
    "att_detection": {
        "att_quality_weight": 1,
        "integrase_proximity_weight": 0.6,
        "model_coverage_weight": 0.9,
        "trna_weight": 0,
        "reference_concurrence_weight": 1.5,
        "kmer_size": 5,
        "min_att_score": 2.2,
        "evalue_threshold": 10000,
        "extention_length": 5000,
        "att_sensitivity": 7,
        "blast_sort_key": "bitscore"},
    "phage_homology_search": {
        "evalue": 1E-04,
        "probability": 90,
        "coverage": 50,
        "min_products": 5,
        "min_products_strict": 10},
    "phameration": {
        "cmode": 0,
        "cstep": 1,
        "sens": 8,
        "ident": 0.5,
        "cover": 0.80,
        "evalue": 0.001,
        "min_gcs": 0.7},
    "prophage_prediction": {
        "window": 55}}
