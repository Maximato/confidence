from Clss.Extractor.Extractor import Extractor
from Clss.FileSys.AlignWriter import AlignWriter
from Clss.Model.Alignment import Alignment
from Clss.Model.AlignedSeq import AlignedSeq
from Clss.FileSys.ConsensusWriter import ConsensusWriter


class AlignController:

    @staticmethod
    def __get_alignment_from(filename):
        records = Extractor.recs_extractor(filename)
        aligned_seqs = [AlignedSeq(rec.seq) for rec in records]
        return Alignment(aligned_seqs)

    @staticmethod
    def convert_to(filename, outfile, full_length, coloring, fmt="html"):

        align = AlignController.__get_alignment_from(filename)
        consensus = align.get_consensus(full_length=full_length)
        ConsensusWriter(consensus).write(outfile, coloring, fmt)
        #AlignWriter(align).write_html_to(outfile, full_length, coloring)

    @staticmethod
    def convert_in_all_combinations(filename, outdir, prefix):
        align = AlignController.__get_alignment_from(filename)
        consensus_full = align.get_consensus(full_length=True)
        consensus_part = align.get_consensus(full_length=False)

        ConsensusWriter(consensus_full).write_all(outdir, prefix)
        ConsensusWriter(consensus_part).write_all(outdir, prefix)

        # AlignWriter(align).write_all_to(outdir, prefix) #outdir, basename(filename).split(".")[0])

    @staticmethod
    def consensus_with(align_file, seqs_file, outfile, fmt="html"):
        align = AlignController.__get_alignment_from(align_file)
        recs = Extractor.recs_extractor(seqs_file)
        seqs = [r.seq for r in recs]
        new_consensus = align.consensus_with(seqs, full_length=True)
        ConsensusWriter(new_consensus).write(outfile, fmt)
