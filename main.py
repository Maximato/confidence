# confidence calculation
from Clss.Aligns import Aligns
from Clss.Fasta import Fasta
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


files = os.listdir("Groups/Aligns")
for file in files:
    seqs = Fasta.extract_from("Groups/Aligns/" + file)
    align = Aligns(seqs, id=file[0:5], name=file[0:5], description="consensus")
    consensus = align.get_consensus()
    print(consensus["deeps"])

    with open("Groups/Consensuses/" + file[0:5] + ".html", "w") as f:
        f.write(align.get_html_consensus(consensus))
