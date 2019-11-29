from Clss.Controller.AlignController import AlignController
from Clss.Controller.GaController import GaController
from Clss.Controller.RecordsController import RecordsController
import argparse


version = "1.0.0"


def create_parser():
    parser = argparse.ArgumentParser(
        prog="fstage",
        description="Program for sequences manipulations, aligning and finding consensus.",
        epilog="(c) Maximato 2019."
    )

    parser.add_argument("--version", action='version', help="print version of program",
                        version='%(prog)s {}'.format(version))

    subparsers = parser.add_subparsers(dest="command", description="Commands for first argument")

    align_parser = subparsers.add_parser("align", help="running GramAlign for aligning",
                          description="Run GramAlign for aligning.")
    align_parser.add_argument("-i", "--input", help="Input file in fasta", metavar="input")
    align_parser.add_argument("-o", "--output", help="Output file align in fasta", metavar="output")
    align_parser.add_argument("-m", "--mode", help="Single or group running", default="s", metavar="mode")

    convert_parser = subparsers.add_parser("convert", help="convert align file into consensus",
                                           description="Convert align file into consensus")
    convert_parser.add_argument("-i", "--input", help="Input align file", metavar="input")
    convert_parser.add_argument("-o", "--output", help="Output dir", metavar="output")
    convert_parser.add_argument("-p", "--prefix", help="Prefix for output files", metavar="prefix")

    consensus_parser = subparsers.add_parser("consensus", help="find consensus for old consensus and sequences",
                                             description="Find consensus for old consensus and sequences")
    consensus_parser.add_argument("-i", "--input", help="Input align file", metavar="input")
    consensus_parser.add_argument("-s", "--seqs", help="Fasta file with sequences", metavar="seqs")
    consensus_parser.add_argument("-o", "--output", help="Output file", metavar="output")

    group_parser = subparsers.add_parser("group", help="grouping of sequences",
                                         description="Grouping of sequences")
    group_parser.add_argument("-i", "--input", help="Input file", metavar="input")
    group_parser.add_argument("-o", "--output", help="Output dir", metavar="output")
    group_parser.add_argument("--ming", help="Minimal size of genome", metavar="ming")
    group_parser.add_argument("--maxg", help="Maximal size of genome", metavar="maxg")

    filtr_parser = subparsers.add_parser("filtr", help="filtrate of sequences",
                                         description="Filtrate of sequences")
    filtr_parser.add_argument("-i", "--input", help="Input file", metavar="input")
    filtr_parser.add_argument("-o", "--output", help="Output file", metavar="output")
    filtr_parser.add_argument("--organism", help="Organism", metavar="organism")
    filtr_parser.add_argument("--mins", help="Minimal size of sequence", metavar="mins")
    filtr_parser.add_argument("--maxs", help="Maximal size of sequence", metavar="maxs")

    return parser


if __name__ == "__main__":
    ga = GaController()
    ac = AlignController()
    rc = RecordsController()

    parser = create_parser()
    namespace = parser.parse_args()
    if namespace.command == "align":
        input = namespace.input
        output = namespace.output
        if namespace.mode == "s":
            ga.align(input, output)
        if namespace.mode == "g":
            ga.align_groups(input, output)
    elif namespace.command == "convert":
        ac.convert_in_all_combinations(namespace.input, namespace.output, namespace.prefix)
    elif namespace.command == "consensus":
        ac.consensus_with(namespace.input, namespace.seqs, namespace.output)
    elif namespace.command == "group":
        rc.grouping(namespace.input, int(namespace.ming), int(namespace.maxg), namespace.output)
    elif namespace.command == "filtr":
        rc.filtrating(namespace.input, namespace.output, namespace.organism, int(namespace.mins), int(namespace.maxs))
    else:
        parser.print_help()

#ga_controller = GaController()
#ga_controller.align_groups("Data/Groups/TBEV_sequences", "Data/Groups/TBEV_aligns_")

#ac = AlignController()
#ac.convert_in_all_combinations("Data/Groups/TBEV_aligns/cds_aln.fasta", "dir", "cds")

#ac.consensus_with("Data/Groups/TBEV_aligns/cds_aln.fasta", "Data/Groups/TBEV_aligns/AB_aln.fasta", "mm.html")

#rc = RecordsController()
#rc.filtrating("Data/rWNV.fasta", "ttsd", "west nile virus", 9000, 12000)
