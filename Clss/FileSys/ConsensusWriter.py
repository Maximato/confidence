from Clss.FileSys.Extractor import Extractor
from Clss.FileSys.check_directory import check_dir
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from os.path import join


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
        """
        :param consensus: dict shape consensus {"symbols": ['A', ..], "ccls": ['c90', ..] ... }
        """
        self.consensus = consensus

    def __get_html_body(self, coloring="c", ignore_gaps=False, ignore_level=0.9):
        """
        Private function for constructing html body of consensus

        :param coloring: 'c' or 'd'. First for coloring confidence, second for coloring deeps
        :param ignore_gaps: boolean, True or False. Ignoring gaps with high level of confidence
        (with confidence >= ignore_level)
        :param ignore_level: float, level of ignoring gaps
        :return: str, html body
        """
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
            # check ignoring gaps
            if (symbols[i] == "-") and ignore_gaps and (self.consensus["confidences"][i] > ignore_level):
                pass
            # adding symbol of consensus and coloring by CLASS
            else:
                html_body += CLASSES[cls[i]] + symbols[i] + "</span>"
                # add line break
                if br_count % 121 == 0:
                    html_body += "<br>\n"
                br_count += 1
        # add closing tags
        html_body += "\n</body>\n</html>"
        return html_body

    def __get_seq_consensus(self, ignore_gaps=False, ignore_level=0.9):
        """
        Private function for constructing Seq consensus
        :param ignore_gaps: boolean, True or False. Ignoring gaps with high level of confidence
        (with confidence >= ignore_level)
        :param ignore_level: float, level of ignoring gaps
        :return: Seq, consensus
        """
        str_consensus = ""
        for i, symbol in enumerate(self.consensus["symbols"]):
            if symbol == "-" and ignore_gaps and self.consensus["confidences"][i] > ignore_level:
                pass
            else:
                str_consensus += symbol
        return Seq(str_consensus)

    def __write_html_to(self, filename, coloring="c", ignore_gaps=False, ignore_level=0.9):
        """
        Private function for writing html consensus
        :param filename: filename to write
        :param coloring: 'c' or 'd'. First for coloring confidence, second for coloring deeps
        :param ignore_gaps: boolean, True or False. Ignoring gaps with high level of confidence
        (with confidence >= ignore_level)
        :param ignore_level: float, level of ignoring gaps
        """
        html_header = Extractor.extract_html_header(coloring)
        html_body = self.__get_html_body(coloring, ignore_gaps, ignore_level)

        check_dir(filename)

        with open(filename, "w") as f:
            f.write(html_header + html_body)

    def __write_seqrec_to(self, filename, ignore_gaps=False, ignore_level=0.9):
        """
        Private function for writing Seq consensus
        :param filename: filename to write
        :param coloring: 'c' or 'd'. First for coloring confidence, second for coloring deeps
        :param ignore_gaps: boolean, True or False. Ignoring gaps with high level of confidence
        (with confidence >= ignore_level)
        :param ignore_level: float, level of ignoring gaps
        """
        seq = self.__get_seq_consensus(ignore_gaps, ignore_level)
        seq_rec = SeqRecord(seq)
        check_dir(filename)
        SeqIO.write(seq_rec, filename, "fasta")

    def write(self, filename, coloring="c", fmt="html", ignore_gaps=False, ignore_level=0.9):
        """
        Function for writing consensus
        :param filename: filename to write
        :param coloring: 'c' or 'd'. First for coloring confidence, second for coloring deeps
        :param fmt: 'html' or 'fasta', format of consensus
        :param ignore_gaps: boolean, True or False. Ignoring gaps with high level of confidence
        (with confidence >= ignore_level)
        :param ignore_level: float, level of ignoring gaps
        """
        if fmt not in ["html", "fasta"]:
            fmt = "html"
            print("Wrong format! Should be fasta or html.")

        if fmt == "html":
            self.__write_html_to(filename, coloring, ignore_gaps, ignore_level)
        elif fmt == "fasta":
            self.__write_seqrec_to(filename, ignore_gaps, ignore_level)

    def write_all(self, outdir, prefix, igl):
        """
        Function for writing consensus in all combinations
        :param outdir: out directory
        :param prefix: prefix for out files
        :param igl: float, level of ignoring gaps
        """
        self.write(join(outdir, f"{prefix}_CI.html"), "c", "html", ignore_gaps=True, ignore_level=igl)
        self.write(join(outdir, f"{prefix}_DI.html"), "d", "html", ignore_gaps=True, ignore_level=igl)
        self.write(join(outdir, f"{prefix}_I.fasta"), fmt="fasta", ignore_gaps=True, ignore_level=igl)

        self.write(join(outdir, f"{prefix}_CNI.html"), "c", "html", ignore_gaps=False)
        self.write(join(outdir, f"{prefix}_DNI.html"), "d", "html", ignore_gaps=False)
        self.write(join(outdir, f"{prefix}_NI.fasta"), fmt="fasta", ignore_gaps=False)
