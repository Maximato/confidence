from Clss.Extractor.Extractor import Extractor
from Clss.Model.Consensus import Consensus
from Clss.FileSys.RecordsWriter import RecordsWriter
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class ConsensusController:
    @staticmethod
    def convert_to_mutations(html_file, outfile, mut_level):
        html_consensus = Extractor.extract_html_consensus(html_file)
        str_cons = Consensus.get_consensus_with_mut(html_consensus, mut_level)
        record = SeqRecord(Seq(str_cons))
        RecordsWriter(record).write_to(outfile)


cc = ConsensusController()
cc.convert_to_mutations("..\..\consensus.html", "sdf", ["c80", "c90"])
