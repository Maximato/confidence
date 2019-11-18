import os


class GramAlign:
    @staticmethod
    def run_gram_align(filename):
        command = "GramAlign -i " + filename + " -o " + filename[:-6] + "_aln.fasta -f 2"
        print(command)
        #os.system(command)


ga = GramAlign()
ga.run_gram_align("file.fasta")
