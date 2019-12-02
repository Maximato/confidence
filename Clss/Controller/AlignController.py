from Clss.FileSys.Extractor import Extractor
from Clss.Model.Alignment import Alignment
from Clss.Model.AlignedSeq import AlignedSeq
from Clss.FileSys.ConsensusWriter import ConsensusWriter
from Bio.SeqIO import SeqRecord
from Bio.Seq import Seq
from Clss.FileSys.RecordsWriter import RecordsWriter
from os.path import basename


class AlignController:

    @staticmethod
    def __get_alignment_from(filename):
        records = Extractor.extract_records(filename)
        aligned_seqs = [AlignedSeq(rec.seq) for rec in records]
        return Alignment(aligned_seqs)

    @staticmethod
    def convert_to(filename, outfile, full_length, coloring="c", fmt="html", ignore_gaps=False, ignore_level=0.9):

        align = AlignController.__get_alignment_from(filename)
        consensus = align.get_consensus(full_length=full_length)
        ConsensusWriter(consensus).write(outfile, coloring, fmt, ignore_gaps, ignore_level)

    @staticmethod
    def convert_in_all_combinations(filename, outdir, prefix):
        align = AlignController.__get_alignment_from(filename)
        consensus_full = align.get_consensus(full_length=True)
        consensus_part = align.get_consensus(full_length=False)

        ConsensusWriter(consensus_full).write_all(outdir, prefix, 0.9)
        ConsensusWriter(consensus_part).write_all(outdir, prefix, 0.9)

    @staticmethod
    def consensus_with(align_file, seqs_file, outfile, fmt="html"):
        align = AlignController.__get_alignment_from(align_file)
        recs = Extractor.extract_records(seqs_file)
        seqs = [r.seq for r in recs]
        new_consensus = align.consensus_with(seqs, full_length=True)
        ConsensusWriter(new_consensus).write(outfile, fmt)

    @staticmethod
    def unite_aligns(dirname, outfile, full_length, ignore_gaps=False, ignore_level=0.9):
        filenames = Extractor.extract_filenames(dirname)
        records = []
        for file in filenames:
            align = AlignController.__get_alignment_from(file)
            consensus = align.get_consensus(full_length=full_length)
            str_cons = consensus.get_str_consensus(ignore_gaps=ignore_gaps, ignore_level=ignore_level)
            name = f"{basename(file).split('.')[0]}:"
            description = f"consensus sequence with parameters: full_length={full_length}, ignore_gaps={ignore_gaps}, " \
                          f"gnore_level={ignore_level}"
            records.append(SeqRecord(Seq(str_cons), id=name, name=name, description=description))
        RecordsWriter(records).write_to(outfile)
