from Clss.Controller.AlignController import AlignController
from Clss.Controller.GaController import GaController
from Clss.Controller.RecordsController import RecordsController

from Clss.Extractor.Extractor import Extractor

from Clss.Model.Counts import Counts
from Clss.Model.Alignment import Alignment
from Clss.Model.Consensus import Consensus
from Clss.Model.AlignedSeq import AlignedSeq
from Clss.Model.Records import Records


#ga_controller = GaController()
#ga_controller.align_groups("Data/Groups/TBEV_sequences", "Data/Groups/TBEV_aligns_")

#ac = AlignController()
#ac.convert_in_all_combinations("Data/Groups/TBEV_aligns/cds_aln.fasta", "dir", "cds")

#ac.consensus_with("Data/Groups/TBEV_aligns/cds_aln.fasta", "Data/Groups/TBEV_aligns/AB_aln.fasta", "mm.html")

rc = RecordsController()
rc.filtrating("Data/rWNV.fasta", "ttsd", "west nile virus", 9000, 12000)
