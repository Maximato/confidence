from os.path import join, basename, abspath, dirname
import os

ga_path = "~/GRAMALIGN/src/GramAlign"


class GramAlign:
    @staticmethod
    def run_gram_align(filename, outfile):
        dname = dirname(outfile)
        # creating output directory
        if not os.path.isdir(dname):
            os.mkdir(dname)

        command = ga_path + " -i " + filename + " -o " + outfile + " -f 2"
        print(command)
        os.system(command)

    @staticmethod
    def run_for_all_in(directory, odir=""):
        ga = GramAlign()
        files = os.listdir(directory)
        for filename in files:
            path = join(directory, filename)
            outfile = join(odir, f"{filename.split('.')[0]}_aln.fasta")
            ga.run_gram_align(path, outfile)


"""l = "dfv.fff"
print(l.split("."))"""

"""
# testing
ga = GramAlign()
ga.run_gram_align("sdgf", ".")


path = join(os.path.abspath("."), "edfs")
print(path)
print(os.path.basename(path))
print(os.path.basename("PycharmProjects\confidence\Clss\edfs"))
"""

