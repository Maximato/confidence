"""
-----------------------------------------------------------------------------------
Script align.py
-----------------------------------------------------------------------------------

The program is designed to run multiple sequence alignment of one of the three
programs (GramAlign, muscle or clustalo). These programs should be installed
on your computer and added to path variable for correct work. If programs were
not be added to 'path' system variable they can be added to definitions.py.
The launch is as follows:

    python3 align.py -i sequences.fasta –p program_name -o out_file.fasta -m

-i, --input (sequences.fasta) - a set of sequences for which multiple alignment
is necessary (or the name of the directory with the files with sequences
stored in it);
-p, --program (program_name) - (optional parameter, default GramAlign)
name of the program that performs alignment: GramAlign, muscle or clustalo
-o, --output (out_file.fasta) - name of the output file (or directory);
–m, --multiply - an optional parameter for automatically starting alignments
of all files contained in the directory specified in the -i (--input) parameter.
"""

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
                                              "multiply mode", )
    parser.add_argument("-p", "--program", help="Name of program that will be used for aligning. "
                                                "One of the GramAlign, muscle of clustalo", default="GramAlign")
    parser.add_argument("-o", "--output", help="Output filename to save aligning or output directory fo multiply mode",
                        metavar="output")
    parser.add_argument("-m", "--multiply", action="store_true", help="Turn on multiply mode for aligning all data "
                                                                      "in directory", )
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
