from datetime import date

from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord


# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
DATE = date.otday.strftime("%d-%b-%Y").upper()

ANNOTATIONS = {"molecule_type": "DNA", "topology": "linear",
               "data_file_division": "PHG", "date": DATE,
               "accessions": [], "sequence_version": "1",
               "keywords": [], "source": "",
               "organism": "", "taxonomy": [],
               "comment": [""]}


def realign_subrecord(record, subrecord, subrecord_start, subrecord_end):
    for feature in record.features:
        if (feature.location.start < subrecord_start or
                feature.location.end > subrecord_end):
            continue

        feature_start = (feature.location.start - subrecord_start)
        feature_end = (feature.location.end - subrecord_start)
        feature_location = FeatureLocation(feature_start, feature_end)

        new_feature = SeqFeature(feature_location, strand=feature.strand,
                                 type=feature.type,
                                 qualifiers=feature.qualifiers)
        subrecord.features.append(new_feature)


class Prophage:
    def __init__(self, parent_record, seq_id, start=None, end=None, strand=1):
        self.parent_record = parent_record
        self.parent_seq = parent_record.seq
        self.id = seq_id

        self.start = None
        self.end = None
        self.strand = strand

        self.feature = None
        self.seq = None
        self.record = None

        self.update_sequence_attributes()

    def set_coordinates(self, start, end):
        self.start = start
        self.end = end

        self.update_sequence_attributes()

    def set_strand(self, strand):
        self.strand = strand

        self.update_sequence_attributes()

    def update_sequence_attributes(self):
        if self.start is None or self.end is None:
            return

        self.feature = SeqFeature(FeatureLocation(self.start, self.end),
                                  strand=1,
                                  type="source")
        self.seq = self.feature.extract(self.parent_seq)
        self.record = SeqRecord(self.seq, id=self.id, annotations=ANNOTATIONS)
        realign_subrecord(self.parent_record, self.record,
                          self.start, self.end)
