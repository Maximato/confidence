class Counts:

    def __init__(self, n):
        self.n = n
        self.counts = {
            "A": [0 for _ in range(n)], "C": [0 for _ in range(n)],
            "G": [0 for _ in range(n)], "T": [0 for _ in range(n)],
            "-": [0 for _ in range(n)]
        }

    def count(self, nucl, i):
        counts = self.counts
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
