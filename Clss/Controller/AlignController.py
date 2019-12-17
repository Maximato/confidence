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
        """
        Private function for extracting alignment from file fasta format. Sequences in file should be the same length
        for correct work

        :param filename: fasta file with sequences the same length, is alignment in fasta format.
        :return: Alignment object
        """
        records = Extractor.extract_records(filename)
        aligned_seqs = [AlignedSeq(rec.seq) for rec in records]
        return Alignment(aligned_seqs)

    @staticmethod
    def convert_to(filename, outfile, full_length, coloring="c", fmt="html", ignore_gaps=False, ignore_level=0.9):
        """
        The function for converting align in fasta format to consensus. All parameters is needed for representation
        of consensus.

        :param filename: fasta file with sequences the same length, is alignment in fasta format.
        :param outfile: output filename, consensus in '.fasta' or '.html' format
        :param full_length: boolean, True or False. If True the confidence level and deeps of position calculating
        will be done using all length of sequences. For False parameter calculating will be done without full gaps
        ends of sequences
        :param coloring: 'c' or 'd'. First for coloring confidence, second for coloring deeps. Coloring used in 'html'
        format outfile only
        :param fmt: format of outfile ('html' or 'fasta')
        :param ignore_gaps: boolean, True or False. Ignoring gaps with high level of confidence
        (with confidence >= ignore_level)
        :param ignore_level: float, level of ignoring gaps
        """
        align = AlignController.__get_alignment_from(filename)
        consensus = align.get_consensus(full_length=full_length)
        ConsensusWriter(consensus).write(outfile, coloring, fmt, ignore_gaps, ignore_level)

    @staticmethod
    def convert_in_all_combinations(filename, outdir, prefix):
        """
        The function for converting align in fasta format to consensus in all combinations of parameters.
        :param filename: fasta file with sequences the same length, is alignment in fasta format
        :param outdir: output directory for consensuses
        :param prefix: prefix of files in outdir
        """
        align = AlignController.__get_alignment_from(filename)
        consensus_full = align.get_consensus(full_length=True)
        consensus_part = align.get_consensus(full_length=False)

        # default level of ignoring gaps
        lig = 0.9

        # writing consensuses
        ConsensusWriter(consensus_full).write_all(outdir, prefix, lig)
        ConsensusWriter(consensus_part).write_all(outdir, prefix, lig)

    @staticmethod
    def consensus_with(align_file, seqs_file, outfile, fmt="html"):
        """
        Get new consensus based on aligning of complete genomes and short sequences from database. This function
        calculate consensus taking into account all data from nucleotide database (genomes and other shorter sequences)
        and save it in file. Aligning of genomes is used as first approximation of consensus that modified in running
        program.

        :param align_file: filename of aligning complete genomes in fasta format
        :param seqs_file: filename with short sequences in fasta
        :param outfile: outfile name
        :param fmt: format of outfile name ('fasta' or 'html')
        """
        align = AlignController.__get_alignment_from(align_file)
        recs = Extractor.extract_records(seqs_file)
        seqs = [r.seq for r in recs]
        new_consensus = align.consensus_with(seqs, full_length=True)
        ConsensusWriter(new_consensus).write(outfile, fmt)

    @staticmethod
    def unite_aligns(dirname, outfile, full_length, ignore_gaps=False, ignore_level=0.9):
        """
        Uniting consensuses of aligning from directory into one fasta file.

        :param dirname: directory name with aligns
        :param outfile: out filename
        :param full_length: boolean, True or False. If True the confidence level will be done using all length of
        sequences. For False parameter calculating will be done without 'full gap' ends of sequences
        :param ignore_gaps: boolean, True or False. Ignoring gaps with high level of confidence
        (with confidence >= ignore_level)
        :param ignore_level: float, level of ignoring gaps
        """
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
