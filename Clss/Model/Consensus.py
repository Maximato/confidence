class Consensus(dict):
    def get_str_consensus(self, ignore_gaps=False, ignore_level=0.9):
        str_consensus = ""
        for i, symbol in enumerate(self["symbols"]):
            if symbol == "-" and ignore_gaps and self["confidences"][i] > ignore_level:
                pass
            else:
                str_consensus += symbol
        return str_consensus
