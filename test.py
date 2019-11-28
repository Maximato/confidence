# confidence calculation
from Clss.Model.Alignment import Alignment
from Clss.Model.Records import Records

sequences = Records()
cds_seqs = sequences.extract_from("Data/Groups/TBEV_aligns/cds_aln.fasta")
align = Alignment(cds_seqs)
counts = align.__get_counts(full_length=True)
consensus = align.get_consensus(counts)
with open("consensus.html", "w") as f:
    f.write(align.get_html_consensus(consensus))

organism = "tick-borne encephalitis virus"
all_seqs = sequences.extract_from("Data/TBEV.fasta")
short_seqs = sequences.filtr_organism_by_size(all_seqs, organism, 100, 10000)

new_counts = align.local_align(counts, short_seqs)
new_consensus = align.get_consensus(new_counts, full_length=True)
html = align.get_html_consensus(new_consensus)
with open("new_consensus.html", "w") as f:
    f.write(html)


#align = pairwise2.align.localxs("TACGCATCGACG-ACTGGGGGAA", "ACGC-ATCG", -2, -1, one_alignment_only=1)
#time.

"""align = read(open("Data/Aligns_compare/TBEV_clustalo_align.phylip", "r"), "phylip")
with open("Data/Aligns_compare/TBEV_clustalo_align.fasta", "w") as f:
    f.write(align.format("fasta"))"""

#align = pairwise2.align.globalxs(seqs[0].seq, seqs[3].seq, -2, -1)
#print(align)
#print(*align[0])
#print(format_alignment(*align[0]))
