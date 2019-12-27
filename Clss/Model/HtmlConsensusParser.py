from abc import ABC
from html.parser import HTMLParser


class HtmlConsensusParser(HTMLParser, ABC):
    """
    Class for converting HTML consensus into string consensus with mutations
    """

    def __init__(self):
        super().__init__()
        self.levels_of_confidence = None
        self.consensus_string = ""
        self.confidence_string = ""
        self.is_confidence = False

    def handle_starttag(self, tag, attrs):
        if tag == "span" and len(attrs) == 1:
            if attrs[0][1] in self.levels_of_confidence:
                self.is_confidence = True
            else:
                self.is_confidence = False

    def handle_data(self, data):
        if data in "ACGT-":
            self.consensus_string += data
            if self.is_confidence:
                self.confidence_string += "*"
            else:
                self.confidence_string += "-"

    def parse_html_consensus(self, html_consensus, levels_of_confidence):
        """
        Parse html consensus and convert it into string with mutations marked as '*'
        :param html_consensus: string, html consensus
        :param levels_of_confidence: list, classes of confidence that we do not consider be mutations ['c90', 'c80']
        :return: string, consensus with mutations
        """
        self.levels_of_confidence = levels_of_confidence
        self.consensus_string = ""
        self.confidence_string = ""
        self.feed(html_consensus)
