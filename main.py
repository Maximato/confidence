# confidence calculation
from Clss.Aligns import Aligns
from Clss.Fasta import Fasta
from Clss.Sequences import Sequences
from Clss.GramAlign import GramAlign
from Bio import SeqIO
import os
from os.path import join

"""
all_nucl = "Data/TBEV_all_nucleotides.fasta"

# grouping nucleotides
sequences = Sequences()
records = sequences.extract_from(all_nucl)
groups = sequences.group_seqs(records)
sequences.write_groups(groups, "Groups/Sequences")

# running GramAlign for all groped sequences
ga = GramAlign()
ga.run_gram_align("file.fasta")
ga.run_for_all_in("Groups/Sequences", "Groups/Aligns")
"""


def run(directory, odir):
    files = os.listdir(directory)
    for file in files:
        seqs = Fasta.extract_from(join(directory, file))
        align = Aligns(seqs, id=file[0:5], name=file[0:5], description="consensus")
        consensus = align.get_consensus()
        print(consensus["deeps"])

        if not os.path.isdir(odir):
            os.mkdir(odir)

        # writing html
        with open(odir + file[0:5] + ".html", "w") as f:
            f.write(align.get_html_consensus(consensus))

        # writing seq record as fasta
        fn = odir + file[0:5] + ".fasta"
        SeqIO.write(align.get_seq_record_consensus(consensus), fn, "fasta")


run("Groups/Aligns", "Groups/Consensuses/")
