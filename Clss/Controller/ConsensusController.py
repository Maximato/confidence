from Clss.FileSys.Extractor import Extractor
from Clss.Model.Consensus import Consensus
from Clss.FileSys.RecordsWriter import RecordsWriter
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class ConsensusController:
    @staticmethod
    def convert_to_mutations(html_file, outfile, mut_level):
        """
        Converting html consensus into '.txt' consensus with mutations. Mutations is positions with low level
        of confidence (in this position high probability to find different nucleotides) marked as '*'

        :param html_file: filename with consensus
        :param outfile: out filename
        :param mut_level: list with classes, that we considered as 'reliable' position. All other positions will
        mark as mutations. Format: ['c90', 'c80', ... ]
        """
        html_consensus = Extractor.extract_html_consensus(html_file)
        str_cons = Consensus.get_consensus_with_mut(html_consensus, mut_level)
        record = SeqRecord(Seq(str_cons))
        RecordsWriter(record).write_to(outfile)
