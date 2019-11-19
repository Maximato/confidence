import os


class GramAlign:
    @staticmethod
    def run_gram_align(filename):
        command = "GramAlign -i " + filename + " -o " + filename[:-6] + "_aln.fasta -f 2"
        print(command)
        #os.system(command)

    @staticmethod
    def run_for_all_in(directory):
        ga = GramAlign()
        files = os.listdir(directory)
        for file in files:
            ga.run_gram_align(file)


ga = GramAlign()
ga.run_gram_align("file.fasta")
ga.run_for_all_in("Groups")
