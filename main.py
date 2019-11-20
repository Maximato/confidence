from Clss.Aligns import Aligns
from Clss.Fasta import Fasta
from Clss.Sequences import Sequences
from Clss.GramAlign import GramAlign
import os
from os.path import join


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


def aligns_to_consensus(aligns_dir, outdir):
    files = os.listdir(aligns_dir)
    for file in files:
        seqs = Fasta.extract_from(join(aligns_dir, file))
        align = Aligns(seqs, id=file[0:5], name=file[0:5], description="consensus")
        consensus = align.get_consensus()
        print(consensus["deeps"])

        if not os.path.isdir(outdir):
            os.mkdir(outdir)

        # writing html
        with open(join(outdir, file[0:5] + ".html"), "w") as f:
            f.write(align.get_html_consensus(consensus))

        # writing seq record as fasta
        fn = join(outdir, file[0:5] + ".fasta")
        align.write_to(align.get_seq_record_consensus(consensus), fn)


# filtrating
organism = "tick-borne encephalitis virus"
filtr_nucl_by("Data/TBEV.fasta", "Data/fTBEV.fasta", organism)

#run("Groups/Aligns", "Groups/Consensuses/")
