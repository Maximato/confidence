from Clss.Extractor.Extractor import Extractor
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from os.path import join, basename, abspath, dirname, isdir
import os

CLASSES = {
    "c00": '<span class="c00">',
    "c01": '<span class="c01">',
    "c10": '<span class="c10">',
    "c20": '<span class="c20">',
    "c30": '<span class="c30">',
    "c40": '<span class="c40">',
    "c50": '<span class="c50">',
    "c60": '<span class="c60">',
    "c70": '<span class="c70">',
    "c80": '<span class="c80">',
    "c90": '<span class="c90">'
}


class ConsensusWriter:
    def __init__(self, consensus):
        self.consensus = consensus

    def __get_html_body(self, coloring="c", ignore_gaps=False, ignore_level=0.9):
        symbols = self.consensus["symbols"]
        if coloring == "c":
            cls = self.consensus["ccls"]
        elif coloring == "d":
            cls = self.consensus["dcls"]
        else:
            raise ValueError("coloring should be 'c' (for confidence) or 'd' (for deeps)")

        n = len(symbols)
        html_body = ""
        br_count = 1
        for i in range(n):
            if (symbols[i] == "-") and ignore_gaps and (self.consensus["confidences"][i] > ignore_level):
                pass
            else:
                html_body += CLASSES[cls[i]] + symbols[i] + "</span>"
                # add line break
                if br_count % 121 == 0:
                    html_body += "<br>\n"
                br_count += 1
        html_body += "\n</body>\n</html>"
        return html_body

    def __get_str_consensus(self, ignore_gaps=False, ignore_level=0.9):
        str_consensus = ""
        for i, symbol in enumerate(self.consensus["symbols"]):
            if symbol == "-" and ignore_gaps and self.consensus["confidences"][i] > ignore_level:
                pass
            else:
                str_consensus += symbol
        return str_consensus

    def __get_seq_consensus(self, ignore_gaps=False, ignore_level=0.9):
        s = self.__get_str_consensus(ignore_gaps, ignore_level)
        return Seq(s)

    def __write_html_to(self, filename, coloring="c"):
        html_header = Extractor.get_html_header(coloring)
        html_body = self.__get_html_body(coloring)

        dname = dirname(filename)
        if not isdir(dname):
            os.mkdir(dname)

        with open(filename, "w") as f:
            f.write(html_header + html_body)

    def __write_seqrec_to(self, filename):
        seq = self.__get_seq_consensus()
        seq_rec = SeqRecord(seq)
        print(seq_rec)
        SeqIO.write(seq_rec, filename, "fasta")

    def write(self, filename, coloring="c", fmt="html"):
        if fmt not in ["html", "fasta"]:
            fmt = "html"
            print("Wrong format! Should be fasta or html.")

        if fmt == "html":
            self.__write_html_to(filename, coloring)
        elif fmt == "fasta":
            self.__write_seqrec_to(filename)

    def write_all(self, outdir, prefix):
        self.write(join(outdir, f"{prefix}_TC.html"), "c", "html")
        self.write(join(outdir, f"{prefix}_TC.html"), "d", "html")
        self.write(join(outdir, f"{prefix}_TC.html"), fmt="fasta")
