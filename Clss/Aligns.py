from Bio._py3k import basestring
from Clss.Fasta import Fasta
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

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

HTML_HEADER = "Clss/header.txt"


def get_html_header(file_name):
    with open(file_name, "r") as f:
        header = f.read()
    return header


def get_group(confidence):
    if confidence == 0:
        return "c00"
    elif 0 < confidence < 0.1:
        return "c01"
    elif 0.1 <= confidence < 0.2:
        return "c10"
    elif 0.2 <= confidence < 0.3:
        return "c20"
    elif 0.3 <= confidence < 0.4:
        return "c30"
    elif 0.4 <= confidence < 0.5:
        return "c40"
    elif 0.5 <= confidence < 0.6:
        return "c50"
    elif 0.6 <= confidence < 0.7:
        return "c60"
    elif 0.7 <= confidence < 0.8:
        return "c70"
    elif 0.8 <= confidence < 0.9:
        return "c80"
    elif 0.9 <= confidence:
        return "c90"


class Aligns(Fasta):
    def __init__(self, seqs, id="<unknown id>", name="<unknown name>", description="<unknown description>"):
        if id is not None and not isinstance(id, basestring):
            raise TypeError("id argument should be a string")
        if not isinstance(name, basestring):
            raise TypeError("name argument should be a string")
        if not isinstance(description, basestring):
            raise TypeError("description argument should be a string")
        self._seqs = seqs
        self.n = len(seqs[1].seq)
        self.id = id
        self.name = name
        self.description = description

    def content_calc(self):
        content = {"starts": [], "ends": []}
        for record in self._seqs:
            # calculate start
            i = 0
            while record.seq[i] == "-":
                i += 1
            # calculate end
            j = self.n - 1
            while record.seq[j] == "-":
                j -= 1
            content["starts"].append(i)
            content["ends"].append(j)
        return content

    def counts_calc(self):
        content = self.content_calc()
        # calculating of counts of nucleotides in alignment
        counts = {
            "A": [0 for _ in range(self.n)], "C": [0 for _ in range(self.n)],
            "G": [0 for _ in range(self.n)], "T": [0 for _ in range(self.n)],
            "-": [0 for _ in range(self.n)]
        }

        for l, record in enumerate(self._seqs):
            for d in range(content["starts"][l], content["ends"][l]):
                nucl = record.seq[d]
                if nucl not in "ACGT-":
                    if nucl == "M":
                        counts["A"][d] += 1
                        counts["C"][d] += 1
                    if nucl == "K":
                        counts["G"][d] += 1
                        counts["T"][d] += 1
                    if nucl == "R":
                        counts["A"][d] += 1
                        counts["G"][d] += 1
                    if nucl == "Y":
                        counts["C"][d] += 1
                        counts["T"][d] += 1
                else:
                    counts[nucl][d] += 1
        return counts

    def get_consensus(self):
        counts = self.counts_calc()

        # calculating of confidence from counts
        consensus = {"symbols": [], "deeps": [], "confidences": [], "groups": []}

        for i in range(self.n):
            summ = 0
            max_score = 0
            symbol = "-"
            for key in counts:
                count = counts[key][i]
                summ += count
                if count > max_score:
                    max_score = count
                    symbol = key

            if summ == 0:
                confidence = 1
            else:
                confidence = max_score / summ

            consensus["symbols"].append(symbol)
            consensus["deeps"].append(summ)
            consensus["confidences"].append(confidence)
            consensus["groups"].append(get_group(confidence))
        return consensus

    def get_html_consensus(self, consensus):
        html = get_html_header(HTML_HEADER)

        symbols = consensus["symbols"]
        groups = consensus["groups"]
        for i in range(self.n):
            html += CLASSES[groups[i]] + symbols[i] + "</span>"
            # add line break
            if (i+1) % 121 == 0:
                html += "<br>\n"
        html += "\n</body>\n</html>"
        return html

    def get_str_consensus(self, consensus, ignore_gaps=False, conf_of_gap_ignoring=0.9):
        str_consensus = ""
        for i, symbol in enumerate(consensus["symbols"]):
            if symbol == "-" and ignore_gaps and consensus["confidences"][i] > conf_of_gap_ignoring:
                pass
            else:
                str_consensus += symbol
        return str_consensus

    def get_seq_record_consensus(self, consensus):
        srt_consensus = self.get_str_consensus(consensus)
        return SeqRecord(Seq(srt_consensus), id=self.id, name=self.name, description=self.description)
