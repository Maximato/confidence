# confidence calculation
from Clss.Consensus import Consensus
from Clss.HTML import HTML
from Clss.Aligns import Aligns
from Clss.Sequences import Sequences
from Clss.GramAlign import GramAlign
import os

all_nucl = "Data/TBEV_all_nucleotides.fasta"

# grouping nucleotides
sequences = Sequences()
records = sequences.extract_from(all_nucl)
groups = sequences.group_seqs(records)
sequences.write_groups(groups)

# running GramAlign for all groped sequences
ga = GramAlign()
ga.run_gram_align("file.fasta")
ga.run_for_all_in("Groups/Sequences", "Groups/Aligns")


aligns = Aligns()
files = os.listdir("Groups/Aligns")
for file in files:
    seqs = aligns.extract_from("Groups/Aligns/" + file)
    counts = aligns.counts_calc(seqs)

    consensus = Consensus(counts)
    sdcg = consensus.confidence_calc()
    print(sdcg["deeps"])

    with open("Groups/Consensuses/" + file[0:5] + ".html", "w") as f:
        f.write(HTML().create_consensus(sdcg))
