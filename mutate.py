from Clss.Controller.ConsensusController import ConsensusController
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
    parser.add_argument("-ml", help="Levels of mutations in string: 'c90 c80 ...'", metavar="ml")
    return parser


if __name__ == "__main__":

    parser = create_parser()
    ns = parser.parse_args()

    ConsensusController.convert_to_mutations(ns.input, ns.output, ns.ml.split())
