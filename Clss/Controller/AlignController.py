from Clss.Extractor.Extractor import Extractor
from Clss.FileSys.AlignWriter import AlignWriter


class AlignController:

    @staticmethod
    def convert_to_html(filename, outfile, full_length, coloring):
        align = Extractor.align_extractor(filename)
        AlignWriter(align).write_html_to(outfile, full_length, coloring)

    @staticmethod
    def convert_to_seqrec(filename, outfile, full_length):
        align = Extractor.align_extractor(filename)
        AlignWriter(align).write_seqrec_to(outfile, full_length)

    @staticmethod
    def convert_in_all_combinations(filename, outdir, prefix):
        align = Extractor.align_extractor(filename)
        AlignWriter(align).write_all_to(outdir, prefix) #outdir, basename(filename).split(".")[0])
