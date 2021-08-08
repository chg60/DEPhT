import re

from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# GLOBAL VARIABLES
# -----------------------------------------------------------------------------

# Initial HHresult header identification regular expressions
# ==========================================================
HHRESULT_HEADERS = {
                "query_id": re.compile(r"Query         ([\w\W\s\S]+)\n"),
                "match_cols": re.compile(r"Match_columns ([\d]+)\n"),
                "num_seqs": re.compile(r"No_of_seqs    (\d+ out of \d+)\n"),
                "neff": re.compile(r"Neff          ([\d\.]+)\n"),
                "searched_HMMs": re.compile(r"Searched_HMMs ([\d]+)\n"),
                "date": re.compile(r"Date          ([\w\W\s\S]+)\n"),
                "command": re.compile(r"Command       ([\w\W\s\S]+)\n")}

MATCHES_HEADER = (re.compile(
        " No Hit                             Prob E-value P-value  Score    "
        "SS Cols Query HMM  Template HMM\n"))
# ==========================================================

# HHresult data table parsing regular expressions
# ===============================================
# Parses one HHresult hit line start into two capturing groups: No and Hit
TABLE_INDEX = re.compile(
         (r"\s*([\d]+)"  # 'No'/'match_num'
          r"\s*([\.;\S\s\-_\d\w]+)$"))  # 'HitID'/'target_id'

# Parsers one HHresult hit line tail into eleven capturing groups
TABLE_DATA = re.compile(
         (r"\s*(\d+\.\d+)"  # 'Prob'/'probability'
          r"\s*([-+\.Ee\d]+)"  # 'E-value'/'expect'
          r"\s*([-+\.Ee\d]+)"  # 'P-value'/'p_value'
          r"\s*(\d+\.\d+)"  # 'Score'/'score'
          r"\s*([\.\d\-]+)"  # 'SS' secondary structure
          r"\s*(\d+)"  # 'Cols'/'match_cols'
          r"\s*(\d+)-(\d+)"  # 'Query HMM'/('query_start')-('query_end')
          # 'Template HMM'/('hit_start')-('hit_end') ('hit_length)
          r"\s+(\d+)-(\d+)\s*\((\d+)\)\n"))
# ===============================================

# HHresult data body parsing regular expressions
# ==============================================
# Parses one HHresult body hit line start into one capturing group
BODY_INDEX = re.compile(r"No (\d+)")  # 'No'/'match_num'
BODY_NAME = re.compile(r">([\S\s\-_\d\w]+)\n")

# Parses one HHresult body hit line start into eight capturing groups
BODY_DATA = re.compile(
                (r"Probab=([\d\.]+)"  # 'probability'
                 r"\s*E-value=([\d\.Ee\-+]+)"  # 'expect'
                 r"\s*Score=([\d\.]+)"  # 'score'
                 r"\s*Aligned_cols=(\d+)"  # 'match_cols"
                 r"\s*Identities=([\d\.]+)%"  # 'pid'
                 r"\s*Similarity=([-\d\.]+)"  # 'similarity'
                 r"\s*Sum_probs=([\d\.]+|inf)"  # 'sum_probs'
                 r"\s*Template_Neff=([\d\.]+)"))  # 'template_Neff'

ALIGNMENT_BLANK_LINES = 2
# Parses one HHresult body alignment line into four capturing groups
BODY_ALIGNMENT_SEQ = re.compile(r"(\w)\s{1}\S+\s*(\d+)\s{1}([\.\w\-]+)"
                                r"\s*(\d+)\s{1}\(\d+\)")
# Parses one HHresult body alignment line into four capturing groups
BODY_ALIGNMENT_CONS = re.compile(r"(\w)\s{1}Consensus\s*(\d+)\s{1}([\.\w\-~]+)"
                                 r"\s*(\d+)\s{1}\(\d+\)")
# Parses one HHresult body alignment line into three capturing groups
BODY_ALIGNMENT_SS = re.compile(r"(\w)\s{1}([_\w]+)\s*([\-\w])")
# Parses one HHresult body alignment line into one capturing group
BODY_ALIGNMENT_MATCH = re.compile(r"([-+|\.\s]+)")
# Parses one HHresult body alignment line into one capturing group
BODY_ALIGNMENT_CONF = re.compile(r"Confidence\s+([\d\s]+)")

BODY_ALIGNMENT_REGEX_MAP = {"cons": BODY_ALIGNMENT_CONS,
                            "seq": BODY_ALIGNMENT_SEQ,
                            "SS": BODY_ALIGNMENT_SS,
                            "match": BODY_ALIGNMENT_MATCH,
                            "conf": BODY_ALIGNMENT_CONF}
# ==============================================

HHMATCH_TABLE_ATTR = ["match_num", "target_id", "probability", "expect",
                      "p_value", "score", "SS", "match_cols", "query_start",
                      "query_end", "hit_start", "hit_end", "hit_length"]
HHMATCH_BODY_ATTR = ["target_id", "probability", "expect", "score",
                     "match_cols", "pid", "similarity", "sum_probs",
                     "template_Neff"]

HHMATCH_ATTR = list(set(HHMATCH_TABLE_ATTR + HHMATCH_BODY_ATTR))


# ERROR HANDLING
# -----------------------------------------------------------------------------
HHRESULT_BASE_MESSAGE = ("Encountered improper formatting while parsing "
                         "HHResult file.\n")
HHRESULT_ERROR_MESSAGES = {"header": ("HHResult file header "
                                      "could not be recognized"),
                           "table": ("HHResult file data table "
                                     "could not be recognized"),
                           "body": ("HHResult file match result "
                                    "could not be recognized"),
                           "alignment": ("HHResult file match alignment "
                                         "was improperly formatted")}


class InitializationError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class HHResultFormatError(Exception):
    def __init__(self, e, line, line_num):
        line_traceback = "".join([f"at line {line_num}:\n",
                                  f"Line {line_num}> ", line])
        e = " ".join([e, line_traceback])
        super().__init__(e)


# HHRESULT DATA CLASSES
# -----------------------------------------------------------------------------
class HHResult:
    """Class for handling HHsuite result files.
    """
    def __init__(self, filepath):
        self.__filepath__ = filepath
        self.__initialized__ = False
        self.__lcounter = 1

        self.matches = None

        # Sets object attributes in preparatin for HHresult file headers
        for attr_name in HHRESULT_HEADERS.keys():
            setattr(self, attr_name, None)

    def parse_result(self):
        """Parses HHsuite result file.
        """
        filehandle = open(self.__filepath__, "r")

        self.__lcounter = 0
        self._parse_header(filehandle)
        match_index_map = self._parse_table(filehandle)
        self._parse_body(filehandle, match_index_map)

        filehandle.close()

        self.matches = list(match_index_map.values())

    def _parse_header(self, filehandle):
        """Parses HHsuite result file header.

        :param filehandle: An open filehandle for a HHresult-formatted file
        """
        # For each HHresult file header line parse and store information
        for attr_name, reg_expr in HHRESULT_HEADERS.items():
            # Parses HHresult file header line and takes header_group
            header_line = self.__attempt_explicit_read(
                                        filehandle, reg_expr, ltype="header")
            header_split = re.split(reg_expr, header_line)

            # Sets attribute defined from HHRESULT_HEADERS with
            # expression capturing group
            setattr(self, attr_name, header_split[1])

        self.__lcounter += 1
        filehandle.readline()

    def _parse_table(self, filehandle):
        """Parses HHsuite result file table.

        :param filehandle: An open filehandle for a HHresult-formatted file
        :returns: Dictionary mapping HHresult match index to a HHMatch object
        :rtype: dict{}
        """
        self.__attempt_explicit_read(filehandle, MATCHES_HEADER, ltype="table")

        match_index_map = dict()
        # Reads matches until the tail end of a line that matches the
        # TABLE_DATA regular expression cannot be found
        while True:
            table_line, matches = self.__attempt_read_check(
                                            filehandle, TABLE_DATA)
            if not matches:
                break

            match_data_split = re.split(TABLE_DATA, table_line)
            match_id_split = re.split(TABLE_INDEX,
                                      match_data_split.pop(0))

            # Create match and use parsed information while ignoring
            # line end whitespace
            match = HHMatch(self.query_id)
            match.load_from_table_data(match_id_split[1:3] +
                                       match_data_split[0:11])

            match_index_map[match.match_num] = (match)

        return match_index_map

    def _parse_body(self, filehandle, match_index_map):
        """Parses HHsuite result file body.

        :param filehandle: An open filehandle for a HHresult-formatted file
        :param match_index_map:
        :type match_index_map:
                Dictionary mapping HHresult match index to a HHMatch object
        """
        self.__lcounter += 1
        body_line = filehandle.readline()
        while True:
            if re.match(BODY_INDEX, body_line) is None:
                break

            body_index_split = re.split(BODY_INDEX, body_line)

            # Asserts the body line directly after the body index is a
            # name line
            body_line = self.__attempt_explicit_read(
                                    filehandle, BODY_NAME, ltype="body")
            body_name_split = re.split(BODY_NAME, body_line)

            # Asserts the body line directly after the body index is a
            # data line
            body_line = self.__attempt_explicit_read(
                                    filehandle, BODY_DATA, ltype="body")
            body_data_split = re.split(BODY_DATA, body_line)

            hhmatch = match_index_map[body_index_split[1]]
            # Rewrites match data taken from the body of the HHResult file
            # onto the corresponding match, identified by its match index
            hhmatch.load_from_body_data(body_name_split[1:2] +
                                        body_data_split[1:-1])

            # Parses the rest of the body alignment and stores it
            # in a HHAlignment object
            hhalignment, body_line = self._parse_body_alignment(
                                                filehandle, hhmatch.target_id)
            hhmatch.hhalignment = hhalignment

    def _parse_body_alignment(self, filehandle, target_id):
        """Parses HHsuite result file body alignments.

        :param filehandle: An open filehandle for a HHresult-formatted file
        :param target_id: Identifier for the target of the HHR match
        :type target_id: str
        :returns: Returns an HHAlignment object and the last unparsed line
        :rtype hhalignment: pde_utils.hhsuite.HHAlignment
        :rtype aln_line: str
        """
        hhalignment = HHAlignment(self.query_id, target_id)
        blanks = 0
        while True:
            self.__lcounter += 1
            aln_line = filehandle.readline()

            # If there are an equal to or greater amount of blank file lines
            # than specified by ALIGNMENT_BLANK_LINES assume EOF/alignment end
            if aln_line in ["", "\n"]:
                blanks += 1
                if blanks > ALIGNMENT_BLANK_LINES:
                    break
                continue

            # Try to parse a non-blank file line as an alignment line
            # print(self.__lcounter)
            parsed = hhalignment.parse_alignment(aln_line)
            if not parsed:
                break
            else:
                blanks = 0

        hhalignment.compile_alignment()
        return hhalignment, aln_line

    def __attempt_read_check(self, filehandle, regex):
        """Helper function to read in a new line from a file and check
        whether it matches an expected regular expression

        :param filehandle: An open filehandle for a HHresult-formatted file
        :param regex: A regular expression to check a read-line against
        :type regex: re.Pattern
        """
        self.__lcounter += 1
        line = filehandle.readline()

        return (line, (re.search(regex, line) is not None))

    def __attempt_explicit_read(self, filehandle, regex, ltype=None):
        """Helper function to read in a new line from a file and check
        whether it matches an expected regular expression

        :param filehandle: An open filehandle for a HHresult-formatted file
        :param regex: A regular expression to check a read-line against
        :type regex: re.Pattern
        :param ltype: A HHResult file section descriptor
        :type ltype: str
        """
        line = filehandle.readline()
        if re.search(regex, line) is None:
            e = HHRESULT_BASE_MESSAGE + HHRESULT_ERROR_MESSAGES.get(ltype, "")
            raise HHResultFormatError(e, line, self.__lcounter)

        self.__lcounter += 1

        return line

    def check_initialization(self, caller):
        """Safe programming feature - raise an exception if a client
        tries to perform operations on uninitialized object.
        :param caller: name of the method that called this one
        """
        if not self.initialized:
            m = (f"Cannot call method '{caller}' on uninitialized "
                 "HHResult object")

            raise InitializationError(m)


class HHMatch:
    """Class for handling single matches from a HHsuite HHR file.
    """
    def __init__(self, query_id):
        self.query_id = query_id

        self.hhalignment = None

        # Sets object attributes in preparatin for HHresult file matches
        for i in range(len(HHMATCH_ATTR)):
            setattr(self, HHMATCH_ATTR[i], None)

    def load_from_table_data(self, table_data):
        """Loads match data parsed from the table section of a HHResult file.

        :param table_data: Ordered data parsed from the table section of a HHR
        :type table_data: list
        """
        # For each HHresult table match line column parse and store information
        for i in range(len(table_data)):
            setattr(self, HHMATCH_TABLE_ATTR[i], table_data[i])

    def load_from_body_data(self, body_data):
        """Loads match data parsed from the body section of a HHResult file.

        :param table_data: Ordered data parsed from the body section of a HHR
        :type table_data: list
        """
        # For each HHresult body match line column parse and store information
        for i in range(len(body_data)):
            setattr(self, HHMATCH_BODY_ATTR[i], body_data[i])


class HHAlignment:
    """Class for handling alignments for single matches from a HHsuite HHR
    """
    def __init__(self, query_id, target_id):
        self.query_id = query_id
        self.target_id = target_id

        self.query_seq = None
        self._query_record = None
        self.target_seq = None
        self._target_record = None

        self.query_cons_seq = None
        self._query_cons_record = None
        self.target_cons_seq = None
        self._target_cons_record = None

        self.alignment = None

    # Could be improved in the future to handle HHR SS, match qualifiers
    # and match SS per position
    def parse_alignment(self, line):
        """Parses and stores information from a line read from an HHsuite HHR

        :param line: A line from a HHResult file
        :type line: str
        :return:  Returns if given line was recognized as an HHR alignment line
        :rtype: bool
        """
        parsed = False
        for aln_line_type, regex in BODY_ALIGNMENT_REGEX_MAP.items():
            if re.match(regex, line) is not None:
                parsed = True

                aln_split = re.split(regex, line)[1:-1]
                if aln_line_type == "seq":
                    if aln_split[0] == "Q":
                        self._load_seq("query_seq", aln_split[2])
                    elif aln_split[0] == "T":
                        self._load_seq("target_seq", aln_split[2])
                elif aln_line_type == "cons":
                    if aln_split[0] == "Q":
                        self._load_seq("query_cons_seq", aln_split[2])
                    elif aln_split[0] == "T":
                        self._load_seq("target_cons_seq", aln_split[2])
                break

        return parsed

    def _load_seq(self, seq_attr, seq_str):
        """Stores/appends to alignment sequence data from a HHR match

        :param seq_attr: HHAlignment Seq object attribute to load
        :type seq_attr: str
        :param seq_str: Sequence to load/append
        :type seq_str: str
        """
        seq = getattr(self, seq_attr)

        # Initialize Seq object,
        # with gaps if the match start is mid translation
        if seq is None:
            seq_str = seq_str
        else:
            seq_str = str(seq) + seq_str

        setattr(self, seq_attr, Seq(seq_str))

    def compile_alignment(self):
        """Compiles a MultipleSeqAlignment object from a HHR alignment
        """
        if self.query_seq is None or self.target_seq is None:
            return

        records = list()
        self._query_record = SeqRecord(self.query_seq, id=self.query_id)
        records.append(self._query_record)

        if self.query_cons_seq is not None:
            self._query_cons_record = SeqRecord(
                                        self.query_cons_seq,
                                        id=" ".join(
                                            [self.query_id, "Consensus"]))
            records.append(self._query_cons_record)

        if self.target_cons_seq is not None:
            self._target_cons_record = SeqRecord(
                                        self.target_cons_seq,
                                        id=" ".join(
                                            [self.target_id, "Consensus"]))
            records.append(self._target_cons_record)

        self._target_record = SeqRecord(self.target_seq, id=self.target_id)
        records.append(self._target_record)

        self.alignment = MultipleSeqAlignment(records)
