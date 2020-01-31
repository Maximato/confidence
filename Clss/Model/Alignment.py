from Clss.Model.Counts import Counts
from Clss.Model.Classification import *
from Clss.Model.Consensus import Consensus
from Bio import pairwise2


class Alignment:

    def __init__(self, aligned_seqs, name="undefined"):
        self.aligned_seqs = aligned_seqs
        self.name = name
        self.max_deep = len(aligned_seqs)
        self.n = aligned_seqs[0].n

        # alignment check
        for als in aligned_seqs:
            if als.n != self.n:
                raise AttributeError("Not equal sizes of aligned sequences")

    def __get_counts(self, full_length=True):
        """
        Calculating of counts of nucleotides in alignment

        :param full_length: boolean, taking into account 'gap ends' if True
        :return: Counts, counts of nucleotides ('A', 'T', 'G', 'C', '-') in every position of alignment
        """
        counts = Counts(self.n)
        for als in self.aligned_seqs:
            if full_length:
                for i in range(0, self.n):
                    nucl = als.seq[i]
                    counts.count(nucl, i)
            else:
                for i in range(als.start_content, als.end_content + 1):
                    nucl = als.seq[i]
                    counts.count(nucl, i)
        return counts

    def get_consensus(self, counts=None, full_length=True):
        """
        Calculating of consensus alignment.

        :param counts: Counts, counts of nucleotides ('A', 'T', 'G', 'C', '-') in every position of alignment
        :param full_length: boolean, taking into account 'gap ends' if True
        :return: dict shape, consensus that contains symbols, deeps, confidences in every position of alignment
        """
        if counts is None:
            counts = self.__get_counts(full_length)
        # calculating of confidence from counts
        consensus = Consensus({"symbols": [], "deeps": [], "confidences": [], "ccls": [], "dcls": []})

        n = counts.n
        for i in range(n):
            summ = 0
            max_score = 0
            symbol = "-"
            for key in counts.counts:
                count = counts.counts[key][i]
                summ += count
                if count > max_score:
                    max_score = count
                    symbol = key

            confidence = max_score / summ

            consensus["symbols"].append(symbol)
            consensus["deeps"].append(summ)
            consensus["confidences"].append(confidence)
            consensus["ccls"].append(get_ccls(confidence))
            consensus["dcls"].append(get_dcls(summ))
        return consensus

    def consensus_with(self, seqs, full_length=True):
        """
        Recount consensus using additional information from sequences
        :param seqs: list sequences
        :param full_length: boolean, taking into account 'gap ends' if True
        :return: dict shape, new consensus
        """
        counts = self.__get_counts(full_length)
        consensus = self.get_consensus(counts=counts)
        str_consensus = consensus.get_str_consensus(ignore_gaps=False)

        k = 1
        for seq in seqs:
            print(f"calculated {k} from {len(seqs)}")
            k += 1
            align = pairwise2.align.localxs(str_consensus, seq, -2, -1, one_alignment_only=1)
            template = align[0][1]
            start = align[0][3]
            end = align[0][4]
            for i in range(start-1, end):
                nucl = template[i]
                counts.count(nucl, i)
        return self.get_consensus(counts)
