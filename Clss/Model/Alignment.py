class Alignment:
    @staticmethod
    def count(nucl, counts, i):
        if nucl not in "ACGT-":
            # base represent 2
            if nucl == "W":
                counts["A"][i] += 0.5
                counts["T"][i] += 0.5
            elif nucl == "S":
                counts["G"][i] += 0.5
                counts["C"][i] += 0.5
            elif nucl == "M":
                counts["A"][i] += 0.5
                counts["C"][i] += 0.5
            elif nucl == "K":
                counts["G"][i] += 0.5
                counts["T"][i] += 0.5
            elif nucl == "R":
                counts["A"][i] += 0.5
                counts["G"][i] += 0.5
            elif nucl == "Y":
                counts["C"][i] += 0.5
                counts["T"][i] += 0.5
            # base represent 3
            elif nucl == "B":
                counts["C"][i] += 0.33
                counts["G"][i] += 0.33
                counts["T"][i] += 0.33
            elif nucl == "D":
                counts["A"][i] += 0.33
                counts["G"][i] += 0.33
                counts["T"][i] += 0.33
            elif nucl == "H":
                counts["A"][i] += 0.33
                counts["C"][i] += 0.33
                counts["T"][i] += 0.33
            elif nucl == "V":
                counts["A"][i] += 0.33
                counts["C"][i] += 0.33
                counts["G"][i] += 0.33
            # any nucleotide
            elif nucl == "N":
                counts["A"][i] += 0.25
                counts["C"][i] += 0.25
                counts["G"][i] += 0.25
                counts["T"][i] += 0.25
            # zero nucleotide
            elif nucl == "Z":
                counts["-"][i] += 1
            else:
                print("Warning! Undefined nucleotide")
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
