from Clss.FileSys.ConsensusWriter import ConsensusWriter
from Clss.Model.Counts import Counts
from os.path import join, basename, abspath, dirname, isdir


class AlignWriter:
    def __init__(self, align):
        self.alignment = align

    def write_html_to(self, filename, full_length, coloring):
        counts = Counts(self.alignment.get_counts(full_length))
        consensus = counts.get_consensus()
        cc = ConsensusWriter(consensus)
        cc.write_html_to(filename, coloring)

    def write_seqrec_to(self, filename, full_length):
        counts = Counts(self.alignment.get_counts(full_length))
        consensus = counts.get_consensus()
        cc = ConsensusWriter(consensus)
        cc.write_seqrec_to(filename)

    def write_all_to(self, outdir, prefix):
        # prefix = basename(filename).split('.')[0]
        self.write_html_to(join(outdir, f"{prefix}_TC.html"), True, "c")
        self.write_html_to(join(outdir, f"{prefix}_TD.html"), True, "d")
        self.write_html_to(join(outdir, f"{prefix}_FC.html"), False, "c")
        self.write_html_to(join(outdir, f"{prefix}_FD.html"), False, "d")

        self.write_seqrec_to(join(outdir, f"{prefix}_T.txt"), True)
        self.write_seqrec_to(join(outdir, f"{prefix}_F.txt"), False)
