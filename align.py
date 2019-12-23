from Clss.Controller.RunAlignController import GaController
import argparse

version = "1.0.0"


def create_parser():
    parser = argparse.ArgumentParser(
        prog="align",
        description="One of the parts fstage program. Serves for running align of sequences",
        epilog="(c) Maximato 2019."
    )

    parser.add_argument("--version", action='version', help="print version of program",
                        version='%(prog)s {}'.format(version))

    parser.add_argument("-i", "--input", help="Input data. File in fasta format for single mode or directory for "
                                              "multiply mode",)
    parser.add_argument("-p", "--program", help="Name of program that will be used for aligning. "
                                                "One of the GramAlign, muscle of clustalo", default="GramAlign")
    parser.add_argument("-o", "--output", help="Output filename to save aligning or output directory fo multiply mode",
                        metavar="output")
    parser.add_argument("-m", "--multiply", action="store_true", help="Turn on multiply mode for aligning all data "
                                                                      "in directory",)
    return parser


if __name__ == "__main__":
    a = GaController()

    parser = create_parser()
    ns = parser.parse_args()

    inp = ns.input
    output = ns.output
    if ns.multiply:
        a.align_groups(inp, ns.program, output)
    else:
        a.align(inp, ns.program, output)
