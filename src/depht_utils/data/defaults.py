HHSUITEDB_DEFAULTS = {"name": "essential",
                      "default_product": "hypothetical protein",
                      "consensus_annotation_cutoff": 0.33,
                      "cpus": 1,
                      "min_size": 10}

REF_DB_DEFAULTS = {"name": "references"}

SHELL_DB_DEFAULTS = {"name": "bacterial_genes"}

MODEL_SCHEMA_DEFAULTS = {"reference_db": "blastn",
                         "phage_homologs_db": "hhsearch",
                         "shell_db": "mmseqs",
                         "classifier": "classifier.pkl"}
