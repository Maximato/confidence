# confidence calculation
from Clss.Consensus import Consensus
from Clss.HTML import HTML
from Clss.Aligns import Aligns

file_name1 = "./TBEV_nuc_Aligned.fasta"
file_name2 = "Data/DQ112_aligned.fasta"


aligns = Aligns()
seqs = aligns.extract_from(file_name2)
counts = aligns.counts_calc(seqs)

consensus = Consensus(counts)
sdcg = consensus.confidence_calc()
print(sdcg["deeps"])

with open("consensus.html", "w") as f:
    f.write(HTML().create_consensus(sdcg))
