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


def filtr_nucl_by(nucl_fasta, out_file, organism, minsize=100):
    sequences = Sequences()
    seqs = sequences.extract_from(nucl_fasta)
    fseqs = sequences.filtr_organizm_by_minsize(seqs, organism, minsize)
    sequences.write_to(fseqs, out_file)


def align_to_consensus(align_fasta, outdir=None):
    seqs = Aligns.extract_from(align_fasta)
    name = basename(align_fasta)[:-6]

    if outdir is None:
        outdir = dirname(align_fasta)
    else:
        if not os.path.isdir(outdir):
            os.mkdir(outdir)

    align = Aligns(seqs, id=name, name=name, description="consensus")
    consensus_from_full = align.get_consensus(full_length=True)
    consensus_from_part = align.get_consensus(full_length=False)

    # writing html
    # from full sequences
    with open(join(outdir, name + "_ffcons_conf.html"), "w") as f:
        f.write(align.get_html_consensus(consensus_from_full, "c", ignore_gaps=False))

    with open(join(outdir, name + "_ffcons_deep.html"), "w") as f:
        f.write(align.get_html_consensus(consensus_from_full, "d", ignore_gaps=False))

    with open(join(outdir, name + "_ffcons_conf_ignore_gaps.html"), "w") as f:
        f.write(align.get_html_consensus(consensus_from_full, "c", ignore_gaps=True, ignore_level=0.9))

    # from part sequences
    with open(join(outdir, name + "_fpcons_conf.html"), "w") as f:
        f.write(align.get_html_consensus(consensus_from_part, "c", ignore_gaps=False))

    with open(join(outdir, name + "_fpcons_deep.html"), "w") as f:
        f.write(align.get_html_consensus(consensus_from_part, "d", ignore_gaps=False))

    with open(join(outdir, name + "_fpcons_conf_ignore_gaps.html"), "w") as f:
        f.write(align.get_html_consensus(consensus_from_part, "c", ignore_gaps=True, ignore_level=0.9))

    # writing seq record as fasta
    # from full sequences
    fn1 = join(outdir, name + "_ffcons.fasta")
    align.write_to(align.get_seq_record_consensus(consensus_from_full, ignore_gaps=False), fn1)
    fn2 = join(outdir, name + "_ffcons_ignore_gaps.fasta")
    align.write_to(align.get_seq_record_consensus(consensus_from_full, ignore_gaps=True, ignore_level=0.9), fn2)

    # from part sequences
    fn1 = join(outdir, name + "_fpcons.fasta")
    align.write_to(align.get_seq_record_consensus(consensus_from_part, ignore_gaps=False), fn1)
    fn2 = join(outdir, name + "_fpcons_ignore_gaps.fasta")
    align.write_to(align.get_seq_record_consensus(consensus_from_part, ignore_gaps=True, ignore_level=0.9), fn2)


def aligns_to_consensuses(aligns_dir, outdir):
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    files = os.listdir(aligns_dir)
    for file in files:
        align_to_consensus(join(aligns_dir, file), outdir)


# filtrating
#organism = "west nile virus"
#filtr_nucl_by("Data/rWNV.fasta", "Data/WNV_full_genome.fasta", organism, minsize=1000)

#ga = GramAlign()
#ga.run_gram_align("Data/WNV.fasta")
#ga.run_gram_align("Data/WNV_full_genome.fasta")

#align_to_consensus("Data/Aligns_compare/TBEV_mscl_align.fasta", "Data/Aligns_compare2")
#align_to_consensus("Data/Aligns_compare/TBEV_ga_align.fasta", "Data/Aligns_compare")
#align_to_consensus("Data/Aligns_compare/TBEV_clustalo_align.fasta", "Data/Aligns_compare")

#align_to_consensus("Data/Groups/WNV_aligns/AF_aln.fasta", "Data/clustalo")
#align_to_consensus("Data/WNV_aln.fasta", "Data/WNV")

#run("Groups/Aligns", "Groups/Consensuses/")
