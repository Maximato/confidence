from Clss.Aligns import Aligns
from Clss.Fasta import Fasta
from Clss.Sequences import Sequences
from Clss.GramAlign import GramAlign
import os
from os.path import join, basename, dirname


def grouping(nucl_fasta, outdir):
    sequences = Sequences()
    records = sequences.extract_from(nucl_fasta)
    groups = sequences.group_seqs(records)
    sequences.write_groups(groups, outdir)


def align_groups(groups_dir, align_dir):
    # running GramAlign for all groped sequences
    ga = GramAlign()
    # ga.run_gram_align("file.fasta")
    ga.run_for_all_in(groups_dir, align_dir)


def filtr_nucl_by(nucl_fasta, out_file, organism):
    sequences = Sequences()
    seqs = sequences.extract_from(nucl_fasta)
    fseqs = sequences.filtr_by(seqs, organism)
    sequences.write_to(fseqs, out_file)


def align_to_consensus(align_fasta, outdir=None):
    seqs = Aligns.extract_from(align_fasta)
    name = basename(align_fasta)

    if outdir is None:
        outdir = dirname(align_fasta)

    align = Aligns(seqs, id=name[:-5], name=name[:-5], description="consensus")
    consensus = align.get_consensus()

    # writing html
    with open(join(outdir, name[:-5] + "_cons.html"), "w") as f:
        f.write(align.get_html_consensus(consensus))

    # writing seq record as fasta
    fn = join(outdir, name[:-5] + "_cons.fasta")
    align.write_to(align.get_seq_record_consensus(consensus), fn)


def aligns_to_consensuses(aligns_dir, outdir):
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    files = os.listdir(aligns_dir)
    for file in files:
        align_to_consensus(join(aligns_dir, file), outdir)


"""
# filtrating
organism = "tick-borne encephalitis virus"
filtr_nucl_by("Data/rTBEV.fasta", "Data/TBEV.fasta", organism)
"""

align_to_consensus("Data/TBEV_align.fasta")

#run("Groups/Aligns", "Groups/Consensuses/")
