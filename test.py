# confidence calculation
from Clss.Consensus import Consensus
from Clss.HTML import HTML
from Clss.Aligns import Aligns
from Clss.Sequences import Sequences
from Clss.GramAlign import GramAlign
import os

aligns = Aligns()
seqs = aligns.extract_from("Data/DQ_aligned.fasta")
counts = aligns.counts_calc(seqs)

consensus = Consensus(counts)
sdcg = consensus.confidence_calc()
print(sdcg["deeps"])

with open("Consensuses.html", "w") as f:
    f.write(HTML().create_consensus(sdcg))
