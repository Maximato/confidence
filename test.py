# confidence calculation
from Clss.Aligns import Aligns
from Clss.Sequences import Sequences
from Clss.GramAlign import GramAlign
from Clss.Fasta import Fasta
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.AlignIO import read
import time
import os


sequences = Sequences()
seqs = sequences.extract_from("Data/Groups/TBEV_aligns/cds_aln.fasta")
align = Aligns(seqs)
counts = align.counts_calc(full_length=True)

shseqs = sequences.extract_from("Data/Groups/TBEV_aligns/AF_aln.fasta")
new_counts = align.local_align(counts, shseqs)


#align = pairwise2.align.localxs("TACGCATCGACG-ACTGGGGGAA", "ACGC-ATCG", -2, -1, one_alignment_only=1)
#time.

"""align = read(open("Data/Aligns_compare/TBEV_clustalo_align.phylip", "r"), "phylip")
with open("Data/Aligns_compare/TBEV_clustalo_align.fasta", "w") as f:
    f.write(align.format("fasta"))"""

#align = pairwise2.align.globalxs(seqs[0].seq, seqs[3].seq, -2, -1)
#print(align)
#print(*align[0])
#print(format_alignment(*align[0]))
