from Bio import pairwise2
from copy import deepcopy


class Alignment:
    @staticmethod
    def count(nucl, counts, i):
        if nucl not in "ACGT-":
            if nucl == "M":
                counts["A"][i] += 0.5
                counts["C"][i] += 0.5
            if nucl == "K":
                counts["G"][i] += 0.5
                counts["T"][i] += 0.5
            if nucl == "R":
                counts["A"][i] += 0.5
                counts["G"][i] += 0.5
            if nucl == "Y":
                counts["C"][i] += 0.5
                counts["T"][i] += 0.5
        else:
            counts[nucl][i] += 1

    def __init__(self, aligned_seqs):
        self.aligned_seqs = aligned_seqs
        self.max_deep = len(aligned_seqs)
        self.n = aligned_seqs[0].n

        # alignment check
        for als in aligned_seqs:
            if als.n != self.n:
                raise AttributeError("Not equal sizes of aligned sequences")

    def get_counts(self, full_length=True):
        # calculating of counts of nucleotides in alignment
        counts = {
            "A": [0 for _ in range(self.n)], "C": [0 for _ in range(self.n)],
            "G": [0 for _ in range(self.n)], "T": [0 for _ in range(self.n)],
            "-": [0 for _ in range(self.n)]
        }
        for als in self.aligned_seqs:
            if full_length:
                for i in range(0, self.n):
                    nucl = als.seq[i]
                    self.count(nucl, counts, i)
            else:
                for i in range(als.start_content, als.end_content + 1):
                    nucl = als.seq[i]
                    self.count(nucl, counts, i)
        return counts

    def local_align(self, counts, seqs):
        consensus = self.get_consensus(counts, full_length=True)
        str_consensus = self.get_str_consensus(consensus, ignore_gaps=False)
        # new_consensus = deepcopy(consensus)
        new_counts = deepcopy(counts)
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
                self.count(nucl, new_counts, i)
        return new_counts
