# confidence calculation
from Clss.Aligns import Aligns
from Clss.Sequences import Sequences
from Clss.GramAlign import GramAlign
from Clss.Fasta import Fasta
import os

seqs = Fasta.extract_from("Data/DQ_aligned.fasta")
aligns = Aligns(seqs)
consensus = aligns.get_consensus()

sr = aligns.get_seq_record_consensus(consensus)
print(sr)
