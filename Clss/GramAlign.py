from os.path import join, basename, abspath
import os


class GramAlign:
    @staticmethod
    def run_gram_align(filename, odir=""):
        bname = basename(filename)
        opath = join(abspath(odir), bname[:-6] + "_aln.fasta")
        command = "GramAlign -i " + filename + " -o " + opath + " -f 2"
        print(command)
        #os.system(command)

    @staticmethod
    def run_for_all_in(directory, odir=""):
        ga = GramAlign()
        files = os.listdir(directory)
        for file in files:
            apath = join(abspath(directory), file)
            ga.run_gram_align(apath, odir)


"""
# testing
path = join(os.path.abspath(""), "edfs")
print(path)
print(os.path.basename(path))
print(os.path.basename("PycharmProjects\confidence\Clss\edfs"))
"""

