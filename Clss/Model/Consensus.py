from html.parser import HTMLParser


class MyHTMLParser(HTMLParser):
    """
    Class for converting HTML consensus into string consensus with mutations
    """

    def __init__(self, conf_levels):
        super().__init__()
        self.conf_levels = conf_levels
        self.consensus = ""
        self.conf_data = False

    def handle_starttag(self, tag, attrs):
        if tag == "span" and len(attrs) == 1:
            if attrs[0][1] in self.conf_levels:
                self.conf_data = True
            else:
                self.conf_data = False

    def handle_data(self, data):
        if data in "ACGT-":
            if self.conf_data:
                self.consensus += data
            else:
                self.consensus += "*"


class Consensus(dict):
    def get_str_consensus(self, ignore_gaps=False, ignore_level=0.9):
        """
        Calculating string consensus
        :param ignore_gaps: boolean, True or False. Ignoring gaps with high level of confidence
        (with confidence >= ignore_level)
        :param ignore_level: float, level of ignoring gaps
        :return: string, consensus
        """
        str_consensus = ""
        for i, symbol in enumerate(self["symbols"]):
            if symbol == "-" and ignore_gaps and self["confidences"][i] > ignore_level:
                pass
            else:
                str_consensus += symbol
        return str_consensus

    @staticmethod
    def get_consensus_with_mut(html_cons, mut_levels):
        """
        Parse html consensus and convert it into string with mutations marked as '*'
        :param html_cons: string, html consensus
        :param mut_levels: list, classes of confidence that we do not consider be mutation ['c90', 'c80']
        :return: string, consensus with mutations
        """
        parser = MyHTMLParser(mut_levels)
        parser.feed(html_cons)
        return parser.consensus
