from Clss.FileSys.Extractor import Extractor
from Clss.Model.HtmlConsensusParser import HtmlConsensusParser
from Clss.FileSys.PeMutateConsensusWriter import PeMutateConsensusWriter
from Clss.FileSys.RecordsWriter import RecordsWriter
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class HtmlConsensusController:
    @staticmethod
    def convert_to_mutations(html_file, outfile, levels_of_confidence, cut_from=0, cut_to=None, fmt="fasta"):
        """
        Converting html consensus into '.fasta' or 'primer explorer' format file containing consensus string
        and confidence string with information about mutations. Mutations is positions with low level of confidence
        (in this position high probability to find different nucleotides) marked as '-'

        :param html_file: filename with consensus
        :param outfile: out filename
        :param levels_of_confidence: list with classes, that we considered as 'reliable' position. All other positions
        will mark as mutations. Format: ['c90', 'c80', ... ]
        :param cut_to: cutting consensus from this position
        :param cut_from: cutting consensus to this position
        :param fmt: format of output file: 'fasta' or 'pe' for primer explorer
        """
        html_consensus = Extractor.extract_html_consensus(html_file)

        html_consensus_parser = HtmlConsensusParser()
        html_consensus_parser.parse_html_consensus(html_consensus, levels_of_confidence)

        sequence = html_consensus_parser.consensus_string[cut_from:cut_to]
        consensus = html_consensus_parser.confidence_string[cut_from:cut_to]

        if fmt == "fasta":
            consensus_record = SeqRecord(Seq(sequence), id="sequence", description=f"sequence of {html_file}")
            confidence_record = SeqRecord(Seq(consensus), id="consensus", description=f"consensus with levels:"
                                                                                      f" {levels_of_confidence}")
            RecordsWriter([consensus_record, confidence_record]).write_to(outfile)
        elif fmt == "pe":
            PeMutateConsensusWriter(sequence, consensus).write_in_pe_format(outfile)
        else:
            raise AttributeError("Only 'fasta' or 'pe' output formats are available")
