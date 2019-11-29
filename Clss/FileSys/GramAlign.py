from os.path import join, dirname
from Clss.FileSys.check_directory import check_dir
import os

ga_path = "~/GRAMALIGN/src/GramAlign"


class GramAlign:
    @staticmethod
    def run_gram_align(filename, outfile):
        # creating output directory
        check_dir(filename)

        command = ga_path + " -i " + filename + " -o " + outfile + " -f 2"
        print(command)
        os.system(command)

    @staticmethod
    def run_for_all_in(directory, odir=""):
        files = os.listdir(directory)
        for filename in files:
            path = join(directory, filename)
            outfile = join(odir, f"{filename.split('.')[0]}_aln.fasta")
            GramAlign.run_gram_align(path, outfile)
