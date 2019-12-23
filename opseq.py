from Clss.Controller.RecordsController import RecordsController
from Clss.Controller.ClusterizationController import ClusterisationController
import argparse

version = "1.0.0"


def create_parser():
    parser = argparse.ArgumentParser(
        prog="opseq",
        description="One of the parts fstage program. Serves for manipulating with sequences",
        epilog="(c) Maximato 2019."
    )

    parser.add_argument("--version", action='version', help="print version of program",
                        version='%(prog)s {}'.format(version))

    subparsers = parser.add_subparsers(dest="command", description="Commands for first argument")

    # add parser of grouping function
    group_parser = subparsers.add_parser("group", help="grouping sequences by names. All big sequences (genomes) "
                                                       "with minsog <= size <= maxsog form separate group 'cds'",
                                         description="Grouping sequences by names. All big sequences (genomes) "
                                                     "with minsog <= size <= maxsog form separate group 'cds'")
    group_parser.add_argument("-i", "--input", help="Input file with sequences (.fasta)", metavar="input")
    group_parser.add_argument("-o", "--output", help="Output directory to save groups", metavar="output")
    group_parser.add_argument("--minsog", help="Minimal size of genome", type=int, metavar="minsog")
    group_parser.add_argument("--maxsog", help="Maximal size of genome", type=int, metavar="maxsog")

    # add parser of filtrating function
    filtr_parser = subparsers.add_parser("filtr", help="Filtrating sequences by parameters",
                                         description="Filtrate of sequences by parameters")
    filtr_parser.add_argument("-i", "--input", help="Input file with sequences (.fasta)", metavar="input")
    filtr_parser.add_argument("-o", "--output", help="Output file with filtrated sequences", metavar="output")
    filtr_parser.add_argument("--organism", help="Organism to filtrate", metavar="organism")
    filtr_parser.add_argument("--mins", help="Minimal size of sequence", type=int, metavar="mins")
    filtr_parser.add_argument("--maxs", help="Maximal size of sequence", type=int, metavar="maxs")

    # add parser of clusterisation
    clust_parser = subparsers.add_parser("clust", help="clusterization of sequences. "
                                                       "Visualisation and saving of clusters into files.",
                                         description="Clusterization of sequences. "
                                                     "Visualisation and saving of clusters into files.")
    clust_parser.add_argument("-i", "--input", help="Input file with sequences (.fasta)", metavar="input")
    clust_parser.add_argument("-o", "--output", help="Output directory for saving data", metavar="output")
    clust_parser.add_argument("-e", "--eps", help="The maximum distance between two samples", default=0.5, type=float,
                              metavar="eps")
    clust_parser.add_argument("-s", "--minsamples", help="The number of samples (or total weight) in a neighborhood "
                                                         "for a point to be considered as a core point. "
                                                         "This includes the point itself.", default=2, type=int,
                              metavar="minsamples")
    clust_parser.add_argument("-d", "--dm", help="Distance matrix for input data records", default=None, metavar="dm")
    return parser


if __name__ == "__main__":
    rc = RecordsController()

    parser = create_parser()
    ns = parser.parse_args()

    if ns.command == "group":
        rc.grouping(ns.input, ns.output, ns.minsog, ns.maxsog)
    elif ns.command == "filtr":
        rc.filtrating(ns.input, ns.output, ns.organism, ns.mins, ns.maxs)
    elif ns.command == "clust":
        ClusterisationController.clusterization(ns.input, ns.output, ns.eps, ns.minsamples, ns.dm)
    else:
        parser.print_help()
