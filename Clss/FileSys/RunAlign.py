from Clss.FileSys.check_directory import check_dir
from definitions import *
import os


class RunAlign:
    @staticmethod
    def __run_gram_align(filename, outfile):
        # checking directory of filename
        check_dir(filename)

        if ga_path:
            program = ga_path
        else:
            program = "GramAlign"
        command = program + " -i " + filename + " -o " + outfile + " -f 2"
        print(command)
        os.system(command)

    @staticmethod
    def __run_clustalo(filename, outfile):
        # checking directory of filename
        check_dir(filename)

        if clustalo_path:
            program = clustalo_path
        else:
            program = "clustalo"
        command = program + " -i " + filename + " -o " + outfile
        print(command)
        os.system(command)

    @staticmethod
    def __run_muscle(filename, outfile):
        # checking directory of filename
        check_dir(filename)

        if clustalo_path:
            program = muscle_path
        else:
            program = "muscle"
        command = program + " -in " + filename + " -out " + outfile
        print(command)
        os.system(command)

    @staticmethod
    def run_align(filename, prog, outfile):
        if prog not in ["GramAlign", "muscle", "clustalo"]:
            raise AttributeError(f"'{prog}' is not available program for multiply alignment. "
                                 "Available is: 'GramAlign', 'muscle', 'clustalo'")
        elif prog == "GramAlign":
            RunAlign.__run_gram_align(filename, outfile)
        elif prog == "muscle":
            RunAlign.__run_muscle(filename, outfile)
        elif prog == "clustalo":
            RunAlign.__run_clustalo(filename, outfile)

    @staticmethod
    def run_for_all_in(directory, prog, odir=""):
        files = os.listdir(directory)
        for filename in files:
            path = join(directory, filename)
            outfile = join(odir, f"{filename.split('.')[0]}_aln.fasta")
            RunAlign.run_align(path, prog, outfile)
