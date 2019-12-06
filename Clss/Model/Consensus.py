from html.parser import HTMLParser


class MyHTMLParser(HTMLParser):
    def __init__(self, conf_levels=None):
        super().__init__()
        if conf_levels is None:
            conf_levels = ["c90"]
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
        str_consensus = ""
        for i, symbol in enumerate(self["symbols"]):
            if symbol == "-" and ignore_gaps and self["confidences"][i] > ignore_level:
                pass
            else:
                str_consensus += symbol
        return str_consensus

    @staticmethod
    def get_consensus_with_mut(html_cons, mut_levels):
        parser = MyHTMLParser(mut_levels)
        parser.feed(html_cons)
        return parser.consensus
