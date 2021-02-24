import numpy as np


class Database:
    """
    A Database is a collection of geneids and translations, with methods
    for conveniently accessing data keyed on either translations or geneids,
    and for removing redundancy for faster processing.
    """
    def __init__(self, geneids=tuple(), translations=tuple()):
        """
        Database constructor - assumes that geneids and translations are
        index-paired: e.g. geneids[i] is the label for the sequence at
        translations[i].

        :param geneids: gene labels
        :type geneids: list of str or tuple of str
        :param translations: gene sequences
        :type translations: list of str or tuple of str
        :return: database
        :rtype: Database
        """
        self.all_genes = dict()
        self.non_redundant_genes = dict()
        self.add_genes(geneids, translations)

    def add_genes(self, geneids, translations):
        """
        Database batch submission tool - assumes that geneids and
        translations are index-paired: e.g. geneids[i] is the label for
        the sequence at translations[i].

        Iterates over paired elements and calls `add_gene()` on them.

        :param geneids: gene labels to add
        :type geneids: list or tuple
        :param translations: gene sequences to add
        :type translations: list or tuple
        :return:
        """
        if len(geneids) != len(translations):
            raise ValueError(f"number of geneids ({len(geneids)}) != number of translations ({len(translations)}")

        for geneid, translation in zip(geneids, translations):
            self.add_gene(geneid, translation)

    def add_gene(self, geneid, translation):
        """
        Database single submission tool - assumes that this pair of
        geneid and translation are not already present in the database.

        :param geneid: gene label to add
        :type geneid: str
        :param translation: gene sequence to add
        :type translation: str
        :return:
        """
        if geneid in self.all_genes.keys():
            raise ValueError(f"{geneid} already in database")

        self.all_genes[geneid] = translation
        non_redundant_group = self.non_redundant_genes.get(translation, list())
        non_redundant_group.append(geneid)
        self.non_redundant_genes[translation] = non_redundant_group

    def get_geneids_from_translation(self, translation):
        """
        Returns the list of geneids associated with a translation.

        :param translation: the gene sequence to get labels for
        :return: geneids
        """
        return self.non_redundant_genes[translation]

    def get_translation_from_geneid(self, geneid):
        """
        Returns the translation associated with a geneid.

        :param geneid: the gene label to get a translation for
        :return: translation
        """
        return self.all_genes[geneid]

    def get_non_redundant_array(self):
        temp_list = list()
        max_len = 0
        for translation in self.non_redundant_genes.keys():
            length = len(translation)
            if length > max_len:
                max_len = length
            temp_list.append((translation, length))

        datatype = np.dtype([("Translation", np.unicode, max_len), ("Length", "u4")])
        return np.sort(np.array(temp_list, dtype=datatype), order=["Length", "Translation"])

    def __repr__(self):
        s = [f">{geneids[0]}\n{translation}\n" for translation, geneids in self.non_redundant_genes.items()]
        return "".join(s)

    def __str__(self):
        s = [f">{geneid}\n{translation}\n" for geneid, translation in self.all_genes.items()]
        return "".join(s)

    def __len__(self):
        return len(self.all_genes)


class Pham(Database):
    """
    A Pham is a slightly extended Database implementation, whose main
    difference is that instances of pham are sortable (on length).
    """
    def __init__(self, geneids, translations):
        """
        Constructor for Pham

        :param geneids: the geneids of the sequences assorted into this Pham
        :type geneids: list of str
        :param translations: the sequences assorted into this Pham
        :type translations: list of str
        """
        super().__init__(geneids, translations)

    def get_translations(self):
        return self.non_redundant_genes.keys()

    def __lt__(self, other):
        return len(self) < len(other)
