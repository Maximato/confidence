from Clss.Model.Alignment import Alignment
from Clss.Model.Consensus import Consensus
from Clss.Classification import *
from copy import deepcopy
from Bio import pairwise2


class Counts:
    def __init__(self, counts):
        self.counts = counts
        self.n = len(counts["A"])

        # checking counts length
        for key in counts:
            if len(counts[key]) != self.n:
                raise ValueError("Lists in counts should be the same length")

    def get_consensus(self):
        # calculating of confidence from counts
        consensus = {"symbols": [], "deeps": [], "confidences": [], "ccls": [], "dcls": []}

        for i in range(self.n):
            summ = 0
            max_score = 0
            symbol = "-"
            for key in self.counts:
                count = self.counts[key][i]
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
            consensus["ccls"].append(get_ccls(confidence))
            consensus["dcls"].append(get_dcls(summ))
        return consensus

    def recount_with(self, seqs):
        consensus = self.get_consensus()
        str_consensus = Consensus.get_str_consensus(consensus, ignore_gaps=False)
        # new_consensus = deepcopy(consensus)
        new_counts = deepcopy(self.counts)
        k = 0
        for seq in seqs:
            print(f"calculated {k} from {len(seqs)}")
            k += 1
            align = pairwise2.align.localxs(str_consensus, seq, -2, -1, one_alignment_only=1)
            template = align[0][1]
            start = align[0][3]
            end = align[0][4]
            for i in range(start-1, end):
                nucl = template[i]
                Alignment.count(nucl, new_counts, i)
        return new_counts
