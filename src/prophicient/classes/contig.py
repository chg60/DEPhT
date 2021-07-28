
# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
CODING_FEATURE_TYPES = ["CDS", "tRNA", "tmRNA"]


class Contig:
    """Class to manage Bacterial genome contig sequences and attributes
    """
    def __init__(self, record, seq_id):
        self.record = record
        self.seq = record.seq
        self.id = seq_id

        self.genes = list()
        self.gene_ids = list()

        self.model_scores = list()
        self.mask_bits = list()
        self.hhsearch_scores = list()

        self.update_gene_attributes()

    def update_gene_attributes(self):
        """Catalogs the genes in a bacterial sequence contig's record
        """
        if self.record is None:
            return

        self.genes = list()
        self.gene_ids = list()

        for index, feature in enumerate(self.record.features):
            if feature.type not in CODING_FEATURE_TYPES:
                continue

            gene_id = "_".join([self.id, str(index+1)])

            feature.qualifiers["locus_tag"] = [gene_id]
            feature.qualifiers["gene"] = [str(index+1)]

            if feature.type != "CDS":
                continue

            self.genes.append(feature)
            self.gene_ids.append(gene_id)

    def update_model_scores(self, model_scores):
        if len(model_scores) != len(self.genes):
            raise ValueError("Inputted model scores iterable does not match "
                             "the length of the stored gene attributes.")

        self.model_scores = list()
        for model_score in model_scores:
            self.model_scores.append(model_score)

    def update_mask_bits(self, mask_bits):
        if len(mask_bits) != len(self.genes):
            raise ValueError("Inputted gene bit mask iterable does not match "
                             "the length of the stored gene attributes")

        self.mask_bits = list()
        for bit in mask_bits:
            self.mask_bits.append(bit)

    def update_hhsearch_scores(self, hhsearch_scores):
        if len(hhsearch_scores) != len(self.genes):
            raise ValueError("Inputted HHsearch score iterable does not match "
                             "the length of the stored gene attributes")

        self.hhsearch_scores = list()
        for score in hhsearch_scores:
            self.hhsearch_scores.append(score)

    def fill_hhsearch_scores(self):
        if not (self.genes):
            return

        self.hhsearch_scores = [0.0] * len(self.genes)
