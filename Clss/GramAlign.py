from os.path import join, basename, abspath
import os

ga_path = "~/GRAMALIGN/src/GramAlign"


class GramAlign:
    @staticmethod
    def run_gram_align(filename, odir="."):
        # creating output directory
        if not os.path.isdir(odir):
            os.mkdir(odir)

        bname = basename(filename)
        opath = join(abspath(odir), bname[:-6] + "_aln.fasta")

        command = ga_path + " -i " + filename + " -o " + opath + " -f 2"
        print(command)
        # os.system(command)

    @staticmethod
    def run_for_all_in(directory, odir=""):
        ga = GramAlign()
        files = os.listdir(directory)
        for file in files:
            apath = join(abspath(directory), file)
            ga.run_gram_align(apath, odir)


"""
# testing
ga = GramAlign()
ga.run_gram_align("sdgf", ".")


path = join(os.path.abspath("."), "edfs")
print(path)
print(os.path.basename(path))
print(os.path.basename("PycharmProjects\confidence\Clss\edfs"))
"""

