from Clss.Controller.AlignController import AlignController
from Clss.Controller.GaController import GaController
from Clss.Controller.RecordsController import RecordsController
from Clss.Controller.ConsensusController import ConsensusController
from Clss.Controller.ClusterizationController import ClusterisationController
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

    # add parser of align function of program
    align_parser = subparsers.add_parser("align", help="running GramAlign for aligning",
                                         description="Run GramAlign for aligning.")
    align_parser.add_argument("-i", "--input", help="Input file in fasta", metavar="input")
    align_parser.add_argument("-o", "--output", help="Output file align in fasta", metavar="output")
    align_parser.add_argument("-m", "--mode", help="Single or group running", default="s", metavar="mode")

    # add parser of convert function of program
    convert_parser = subparsers.add_parser("convert", help="convert align file into consensus",
                                           description="Convert align file into consensus")
    convert_parser.add_argument("-i", "--input", help="Input align file", metavar="input")
    convert_parser.add_argument("-o", "--output", help="Output dir", metavar="output")
    convert_parser.add_argument("-p", "--prefix", help="Prefix for output files", metavar="prefix")

    # add parser of consensus function of program
    consensus_parser = subparsers.add_parser("consensus", help="find consensus for old consensus and sequences",
                                             description="Find consensus for old consensus and sequences")
    consensus_parser.add_argument("-i", "--input", help="Input align file", metavar="input")
    consensus_parser.add_argument("-s", "--seqs", help="Fasta file with sequences", metavar="seqs")
    consensus_parser.add_argument("-o", "--output", help="Output file", metavar="output")

    # add parser of uniting aligns into one file
    unite_parser = subparsers.add_parser("unite", help="unite aligns from files into one with consensuses",
                                         description="Unite aligns from files into one with consensuses")
    unite_parser.add_argument("-i", "--input", help="Input dir with files aligning", metavar="input")
    unite_parser.add_argument("-o", "--output", help="Output file with consensuses", metavar="output")
    unite_parser.add_argument("-fl", "--flength", help="Count consensus for full length of aligned sequences or only "
                                                       "for content (in case fl=False)", metavar="flength")
    unite_parser.add_argument("-ig", help="Boolean, ignore gaps in consensus", default="False", metavar="ig")
    unite_parser.add_argument("-il", help="From 0 to 1, level of ignoring gaps in consensus. If "
                                          "confidens of gap in position more or equal ignore level, "
                                          "gap will be ignored. Work only if ignore gaps is True",
                              default="0.9", metavar="il")

    # add parser of grouping function of program
    group_parser = subparsers.add_parser("group", help="grouping of sequences",
                                         description="Grouping of sequences")
    group_parser.add_argument("-i", "--input", help="Input file", metavar="input")
    group_parser.add_argument("-o", "--output", help="Output dir", metavar="output")
    group_parser.add_argument("--ming", help="Minimal size of genome", metavar="ming")
    group_parser.add_argument("--maxg", help="Maximal size of genome", metavar="maxg")

    # add parser of filtrating function of program
    filtr_parser = subparsers.add_parser("filtr", help="filtrate of sequences",
                                         description="Filtrate of sequences")
    filtr_parser.add_argument("-i", "--input", help="Input file", metavar="input")
    filtr_parser.add_argument("-o", "--output", help="Output file", metavar="output")
    filtr_parser.add_argument("--organism", help="Organism", metavar="organism")
    filtr_parser.add_argument("--mins", help="Minimal size of sequence", metavar="mins")
    filtr_parser.add_argument("--maxs", help="Maximal size of sequence", metavar="maxs")

    # add parser of converting to mutations
    mut_parser = subparsers.add_parser("mut", help="convert html consensus to consensus with mutations",
                                       description="Convert html consensus to consensus with mutations")
    mut_parser.add_argument("-i", "--input", help="Input html consensus", metavar="input")
    mut_parser.add_argument("-o", "--output", help="Out file with mutations", metavar="output")
    mut_parser.add_argument("-ml", help="Levels of mutations in string: 'c90 c80 ...'", metavar="ml")

    # add parser of clusterisation
    clust_parser = subparsers.add_parser("clust", help="clusterize sequences", description="clusterize sequences")
    clust_parser.add_argument("-i", "--input", help="Input sequences", metavar="input")
    clust_parser.add_argument("-o", "--output", help="Output dir", metavar="output")
    clust_parser.add_argument("-e", "--eps", help="The maximum distance between two samples", default=1, metavar="eps")
    clust_parser.add_argument("-s", "--minsamples", help="The number of samples (or total weight) in a neighborhood "
                                                         "for a point to be considered as a core point. "
                                                         "This includes the point itself.", default=2,
                              metavar="minsamples")
    clust_parser.add_argument("-d", "--dm", help="Distance matrix for input data records", default=None, metavar="dm")
    return parser


if __name__ == "__main__":
    ga = GaController()
    ac = AlignController()
    rc = RecordsController()
    cc = ConsensusController()
    clc = ClusterisationController()

    parser = create_parser()
    ns = parser.parse_args()
    if ns.command == "align":
        input = ns.input
        output = ns.output
        if ns.mode == "s":
            ga.align(input, output)
        if ns.mode == "g":
            ga.align_groups(input, output)
    elif ns.command == "convert":
        ac.convert_in_all_combinations(ns.input, ns.output, ns.prefix)
    elif ns.command == "consensus":
        ac.consensus_with(ns.input, ns.seqs, ns.output)
    elif ns.command == "unite":
        ac.unite_aligns(ns.input, ns.output, bool(ns.flength), bool(ns.ig),
                        float(ns.il))
    elif ns.command == "group":
        rc.grouping(ns.input, int(ns.ming), int(ns.maxg), ns.output)
    elif ns.command == "filtr":
        rc.filtrating(ns.input, ns.output, ns.organism, int(ns.mins), int(ns.maxs))
    elif ns.command == "mut":
        cc.convert_to_mutations(ns.input, ns.output, ns.ml.split())
    elif ns.command == "clust":
        clc.clusterization(ns.input, ns.output, float(ns.eps), int(ns.minsamples), ns.dm)
    else:
        parser.print_help()
