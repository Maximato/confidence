from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Clss.Classification import *


class Consensus:
    @staticmethod
    def get_html_body(consensus, coloring="c", ignore_gaps=False, ignore_level=0.9):
        symbols = consensus["symbols"]
        if coloring == "c":
            cls = consensus["ccls"]
        elif coloring == "d":
            cls = consensus["dcls"]
        else:
            raise ValueError("coloring should be 'c' (for confidence) or 'd' (for deeps)")

        n = len(symbols)
        html_body = ""
        br_count = 1
        for i in range(n):
            if symbols[i] == "-" and ignore_gaps and consensus["confidences"][i] > ignore_level:
                pass
            else:
                html_body += CLASSES[cls[i]] + symbols[i] + "</span>"
                # add line break
                if br_count % 121 == 0:
                    html_body += "<br>\n"
                br_count += 1
        html_body += "\n</body>\n</html>"
        return html_body

    @staticmethod
    def get_str_consensus(consensus, ignore_gaps=False, ignore_level=0.9):
        str_consensus = ""
        for i, symbol in enumerate(consensus["symbols"]):
            if symbol == "-" and ignore_gaps and consensus["confidences"][i] > ignore_level:
                pass
            else:
                str_consensus += symbol
        return str_consensus

    @staticmethod
    def get_seq_consensus(consensus, ignore_gaps=False, ignore_level=0.9):
        s = Consensus.get_str_consensus(consensus, ignore_gaps, ignore_level)
        return Seq(s)
