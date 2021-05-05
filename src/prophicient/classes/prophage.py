from datetime import date

from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord


# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
DATE = date.today().strftime("%d-%b-%Y").upper()

ANNOTATIONS = {"molecule_type": "DNA", "topology": "linear",
               "data_file_division": "PHG", "date": DATE,
               "accessions": [], "sequence_version": "1",
               "keywords": [], "source": "",
               "organism": "", "taxonomy": [],
               "comment": [""]}

GENE_FEATURES = ["CDS", "tRNA", "tmRNA"]


def realign_subrecord(record, subrecord, subrecord_start, subrecord_end,
                      rev_orient=False):
    """Adds and start-aligns features in the parent record to a subrecord
    given subrecord coordinates with respect to the parent record.
    
    :param record: Parent record to draw features from.
    :type record: Bio.SeqRecord.SeqRecord
    :param subrecord: Record to place and align features into.
    :type subrecord: Bio.SeqRecord.SeqRecord
    :param subrecord_start: Coordinate start in the record to take features.
    :type subrecord_start: int
    :param subrecord_end: Coordinate end in the record to stop taking features.
    :type subrecord_end: int
    :param rev_orient: Toggles inversion of feature strands in the subrecord.
    :type rev_orient: bool
    """
    for feature in record.features:
        if (feature.location.start < subrecord_start or
                feature.location.end > subrecord_end):
            continue

        if rev_orient:
            feature_start = (subrecord_end - feature.location.end)
            feature_end = (subrecord_end - feature.location.start)
            feature_location = FeatureLocation(feature_start, feature_end)

            if feature.strand == 1:
                strand = -1
            else:
                strand = 1
        else:
            feature_start = (feature.location.start - subrecord_start)
            feature_end = (feature.location.end - subrecord_start)
            feature_location = FeatureLocation(feature_start, feature_end)

            strand = feature.strand

        new_feature = SeqFeature(feature_location, type=feature.type,
                                 qualifiers=feature.qualifiers,
                                 strand=strand)
        subrecord.features.append(new_feature)
        subrecord.features.sort(key=lambda x: x.location.start)


class Prophage:
    """Class to manage Prophage sequences and attributes with respect to
    a parent bacterial sequence.
    """
    def __init__(self, parent_record, seq_id, start=None, end=None, att_len=0,
                 strand=1):
        self.parent_record = parent_record
        self.parent_seq = parent_record.seq
        self.id = seq_id

        self.start = start
        self.end = end 
        self.strand = strand
        
        self.length = 0
        self.feature = None
        self.seq = None
        self.record = None

        self.att_len = att_len
        self.attL = None
        self.attR = None

        self.update_sequence_attributes()

    def set_coordinates(self, start, end):
        """Sets the coordinates of the prophage in the bacterial sequence

        :param start: The start (left) of the prophage on the sequence.
        :type start: int
        :param end: The end (right) of the prophage on the sequence.
        :type end: int
        """
        self.start = start
        self.end = end
        
    def set_strand(self, strand):
        """Sets the strand orientation of the prophage in the sequence

        :param strand: The orientation of the prophage in the sequence.
        :type strand: int
        """
        self.strand = strand

        self.update_sequence_attributes()

    def detect_orientation(self):
        """Tries to detect the strand orientation of the prophage in the
        sequence
        """
        if self.record is None:
            return

        orientation = 0
        for feature in self.record.features:
            if feature.type not in GENE_FEATURES:
                continue

            orientation += int(feature.strand) * len(feature)

        if orientation >= 0:
            self.strand = 1
        else:
            self.strand = -1

    def set_att_len(self, att_len):
        """Sets the att length of the prophage in the sequence

        :param att_len: The length of the att sites of the prophage
        :type att_len: int
        """
        self.att_len = att_len

    def update_sequence_attributes(self):
        """Sets prophage SeqFeature and SeqRecord attributes based on
        the object's set coordinates.
        """
        if self.start is None or self.end is None:
            return

        # Sets the prophage SeqFeature according to the stored start and end
        self.feature = SeqFeature(FeatureLocation(self.start, self.end),
                                  strand=self.strand, type="source")
        # Extracts the prophage sequence with the created prophage
        self.seq = self.feature.extract(self.parent_seq)
        self.length = len(self.seq)
        # Creates and realigns a Prophage SeqRecord object
        self.record = SeqRecord(self.seq, id=self.id, annotations=ANNOTATIONS)

        realign_subrecord(self.parent_record, self.record,
                          self.start, self.end,
                          rev_orient=(self.strand != 1) )

    def update_att_attributes(self):
        """Sets prophage attL and attR attributes based on the object's
        set sequence and att metadata.
        """
        if self.seq is None or self.record is None or self.att_len <= 0:
            return

        left_location = FeatureLocation(0, self.att_len)
        right_location  = FeatureLocation(self.length-self.att_len,
                                          self.length)
        
        if self.strand == 1:
            attL_location = left_location
            attR_location = right_location
            strand = 1
        if self.strand == -1:
            attL_location = right_location
            attR_location = left_location
        
        self.attL = SeqFeature(attL_location,
                               strand=self.strand, type="misc_recomb")
        self.attL.qualifiers["note"] = [" ".join([self.id, "attL"])]

        self.attR = SeqFeature(attR_location,
                               strand=self.strand, type="misc_recomb")
        self.attR.qualifiers["note"] = [" ".join([self.id, "attR"])]

        self.record.features.append(self.attL)
        self.record.features.append(self.attR)

        self.record.features.sort(key=lambda x: x.location.start)

    def update(self):
        """Sets prophage sequence and att features based on stored information
        """
        self.update_sequence_attributes()
        self.update_att_attributes()

    def clean_record(self):
        """Add quality-of-life features and tidy prophage extraction artifacts.
        """
        self.record.features.sort(key=lambda x: x.location.start) 

        gene_counter = 0
        gene_features = []
        for feature in self.record.features:
            if feature.type not in GENE_FEATURES:
                continue

            gene_counter += 1

            locus_tag = "_".join([self.id, str(gene_counter)])   
            feature.qualifiers["locus_tag"] = [locus_tag]

            gene_feature = SeqFeature(feature.location,
                                      type="gene")
            gene_feature.qualifiers["locus_tag"] = [locus_tag]
            gene_features.append(gene_feature)

            if feature.type == "CDS":
                assert(
                 str(feature.translate(self.seq,  to_stop=True, table=11)) == \
                 feature.qualifiers["translation"][0])

        source_feature = SeqFeature(FeatureLocation(0, self.length),
                                    type="source")

        self.record.features = self.record.features + gene_features
        self.record.features.append(source_feature)
        
        self.record.features.sort(key=lambda x: x.location.start) 
