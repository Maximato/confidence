# confidence calculation
from Bio import SeqIO
from Consensus import Consensus
from HTML import HTML

file_name1 = "./TBEV_nuc_Aligned.fasta"
file_name2 = "./DQ112_aligned.fasta"
seqs = list(SeqIO.parse(file_name2, "fasta"))

N = len(seqs[1].seq)

# calculating starts and ends of information in sequence
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
    j = N - 1
    while record.seq[j] == "-":
        j -= 1

    content["starts"].append(i)
    content["ends"].append(j)

print(content["starts"])


# calculating of counts of nucleotides in alignment
counts = {
    "A": [0 for _ in range(N)],
    "C": [0 for _ in range(N)],
    "G": [0 for _ in range(N)],
    "T": [0 for _ in range(N)],
    "-": [0 for _ in range(N)]
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

consensus = Consensus(counts)
sdcg = consensus.confidence_calc()
print(sdcg["deeps"])

with open("consensus.html", "w") as f:
    f.write(HTML().create_consensus(sdcg))
