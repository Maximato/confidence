"""
-----------------------------------------------------------------------------------
Script mutate.py
-----------------------------------------------------------------------------------

This script is used to convert a consensus sequence in html format to a text format
with consensus string and confidence string in which all unreliable positions
(with frequent mutations) are marked with a special symbol '-'.

    python3 mutate.py –i consensus.html – o output.fasta –ml c90 c80

-i, --input (consensus.html) - html file containing consensus;
-o, --output (output.fasta) - output file name, fasta format;
-ml - levels of nucleotide occurrence, below which nucleotides are noted as mutations.
-cf - start position for cutting
-ct - end position for cutting
-fmt - format of out file: fasta or pe
"""


from Clss.Controller.HtmlConsensusController import HtmlConsensusController
import argparse

version = "1.0.0"


def create_parser():
    parser = argparse.ArgumentParser(
        prog="mutate",
        description="One of the parts fstage program. Serves for converting html consensus to consensus with mutations",
        epilog="(c) Maximato 2019."
    )

    parser.add_argument("--version", action='version', help="print version of program",
                        version='%(prog)s {}'.format(version))

    parser.add_argument("-i", "--input", help="Input file as html consensus",
                        metavar="input")
    parser.add_argument("-o", "--output", help="Out file with mutations", metavar="output")
    parser.add_argument("-ml", help="Levels of mutations in string: c90 c80 ...", type=str,  metavar="ml", nargs='+')
    parser.add_argument("-cf", help="Start position for cutting", type=int, default=0, metavar="cf")
    parser.add_argument("-ct", help="End position for cutting", type=int, default=None, metavar="ct")
    parser.add_argument("-fmt", help="Format of out file: fasta or pe", default="fasta", metavar="fmt")
    return parser


if __name__ == "__main__":

    parser = create_parser()
    ns = parser.parse_args()
    HtmlConsensusController.convert_to_mutations(ns.input, ns.output, ns.ml, ns.cf, ns.ct, ns.fmt)
