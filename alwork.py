from Clss.Controller.AlignController import AlignController
import argparse

version = "1.0.0"


def create_parser():
    parser = argparse.ArgumentParser(
        prog="alwork",
        description="One of the parts fstage program. Serves for working with aligned sequences",
        epilog="(c) Maximato 2019."
    )

    parser.add_argument("--version", action='version', help="print version of program",
                        version='%(prog)s {}'.format(version))

    subparsers = parser.add_subparsers(dest="command", description="Commands for first argument")

    # 1. add parser of converting aligned sequences
    convert_parser = subparsers.add_parser("convert", help="convert in different combinations aligned sequences "
                                                           "stored in fasta format into consensus",
                                           description="Convert in different combinations aligned sequences "
                                                       "stored in fasta format into consensus")
    convert_parser.add_argument("-i", "--input", help="Input align file", metavar="input")
    convert_parser.add_argument("-o", "--output", help="Output directory", metavar="output")
    convert_parser.add_argument("-p", "--prefix", help="Prefix for output files", metavar="prefix")

    # 2. add parser of creating consensus based on aligned sequences of genome and short sequences (not complete genome)
    consensus_parser = subparsers.add_parser("consensus",
                                             help="create consensus based on aligned sequences of complete genome "
                                                  "and short sequences (not complete genome)",
                                             description="Create consensus based on aligned sequences of complete "
                                                         "genome and short sequences (not complete genome)")
    consensus_parser.add_argument("-i", "--input", help="Input align file", metavar="input")
    consensus_parser.add_argument("-s", "--seqs", help="Fasta file with sequences", metavar="seqs")
    consensus_parser.add_argument("-o", "--output", help="Output file", metavar="output")
    consensus_parser.add_argument("-f", "--format", help="Format of output file ('html', 'fasta')", default="html",
                                  metavar="format")

    # 3. add parser of uniting aligns into one file
    unite_parser = subparsers.add_parser("unite", help="unite aligns from fasta files into one as consensuses",
                                         description="Unite aligns from fasta files into one as consensuses")
    unite_parser.add_argument("-i", "--input", help="Input directory with aligns", metavar="input")
    unite_parser.add_argument("-o", "--output", help="Output file with consensuses", metavar="output")
    unite_parser.add_argument("-fl", "--flength",
                              help="Count consensus for full length of aligned sequences (fl=True) or "
                                   "only for content (fl=False)", default=True, type=bool, metavar="flength")
    unite_parser.add_argument("-ig", help="Boolean, ignore gaps in consensus", default=False, type=bool, metavar="ig")
    unite_parser.add_argument("-il", help="From 0 to 1, level of ignoring gaps in consensus. If "
                                          "confidence of gap in position more or equal ignore level, "
                                          "gap will be ignored. Work only if ignore gaps is True",
                              default=0.9, type=float, metavar="il")
    return parser


if __name__ == "__main__":
    ac = AlignController()

    parser = create_parser()
    ns = parser.parse_args()

    if ns.command == "convert":
        ac.convert_in_all_combinations(ns.input, ns.output, ns.prefix)
    elif ns.command == "consensus":
        ac.consensus_with(ns.input, ns.seqs, ns.output, ns.format)
    elif ns.command == "unite":
        ac.unite_aligns(ns.input, ns.output, ns.flength, ns.ig, ns.il)
    else:
        parser.print_help()
