PARAMETERS = {
    "bacterial_sequences": {
        "phameration": {
            "first_iteration": {
                "--cluster-mode": 0,
                "--cluster-steps": 1,
                "-s": 8,
                "--min-seq-id": 0.50,
                "-c": 0.80,
                "-e": 0.001}},
        "pham_clust": {
            "gcs_threshold": 70}},
    "phage_sequences": {
        "phameration": {
            "first_iteration": {
                "--cluster-mode": 0,
                "--cluster-steps": 1,
                "-s": 8,
                "--min-seq-id": 0.45,
                "-c": 0.75,
                "-e": 0.001},
            "second_iteration": {
                "--min-seq-id": 0.3,
                "-c": 0.5,
                "-e": 0.001}},
        "pham_clust": {
            "gcs_threshold": 35}},
    "phage_homologs": {
        "annotation_consensus_threshold": 0.33,
        "min_hmm_count": 5,
        "essential_annotations": {
            "LIKE": ["terminase",
                     "major capsid",
                     "portal",
                     "lysin",
                     "integrase",
                     "immunity repressor",
                     "transposase",
                     "collar"],
            "NOT LIKE": ["hypothetical",
                         "putative",
                         "NKF"]},
        "extended_annotations": {
            "LIKE": [],
            "NOT LIKE": ["hypothetical",
                         "putative",
                         "NKF",
                         "teminase",
                         "lysin",
                         "immunity repressor",
                         "integrase",
                         "major capsid",
                         "portal",
                         "transposase",
                         "collar"]}},
    "shell_db": {"rep_threshold": 0.60},
    "reference_db": {},
    "classifier": {"window": 55}}
