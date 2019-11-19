from Clss.Fasta import Fasta


class Aligns(Fasta):

    @staticmethod
    def content_calc(seqs):
        n = len(seqs[1].seq)
        content = {
            "starts": [],
            "ends": []
        }

        for record in seqs:
            # calculate start
            i = 0
            while record.seq[i] == "-":
                i += 1

            # calculate end
            j = n - 1
            while record.seq[j] == "-":
                j -= 1

            content["starts"].append(i)
            content["ends"].append(j)

        print(content["starts"])
        return content

    def counts_calc(self, seqs):
        n = len(seqs[1].seq)

        content = self.content_calc(seqs)

        # calculating of counts of nucleotides in alignment
        counts = {
            "A": [0 for _ in range(n)],
            "C": [0 for _ in range(n)],
            "G": [0 for _ in range(n)],
            "T": [0 for _ in range(n)],
            "-": [0 for _ in range(n)]
        }

        for l, record in enumerate(seqs):
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
