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


def realign_subrecord(record, subrecord, subrecord_start, subrecord_end):
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
    """
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
    """Class to manage Prophage sequences and attributes with respect to
    a parent bacterial sequence.
    """
    def __init__(self, parent_record, seq_id, start=None, end=None, strand=1):
        self.parent_record = parent_record
        self.parent_seq = parent_record.seq
        self.id = seq_id

        self.start = start
        self.end = end 
        self.strand = strand

        self.feature = None
        self.seq = None
        self.record = None

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
        
        self.update_sequence_attributes()

    def set_strand(self, strand):
        """Sets the strand orientation of the prophage in the sequence

        :param strand: The orientation of the prophage in the sequence.
        :type strand: int
        """
        self.strand = strand

        self.update_sequence_attributes()

    def update_sequence_attributes(self):
        """Sets prophage SeqFeature and SeqRecord attributes based on
        the object's given coordinates.
        """
        if self.start is None or self.end is None:
            return

        # Sets the prophage SeqFeature according to the stored start and end
        self.feature = SeqFeature(FeatureLocation(self.start, self.end),
                                  strand=self.strand,
                                  type="source")
        # Extracts the prophage sequence with the created prophage
        self.seq = self.feature.extract(self.parent_seq)
        # Creates and realigns a Prophage SeqRecord object
        self.record = SeqRecord(self.seq, id=self.id, annotations=ANNOTATIONS)
        realign_subrecord(self.parent_record, self.record,
                          self.start, self.end)
