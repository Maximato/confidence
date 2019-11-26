from Clss.FileSys.GramAlign import GramAlign


class GaController:

    @staticmethod
    def align(filename, outfile):
        ga = GramAlign()
        ga.run_gram_align(filename, outfile)

    @staticmethod
    def align_groups(groups_dir, align_dir):
        # running GramAlign for all groped sequences
        ga = GramAlign()
        # ga.run_gram_align("file.fasta")
        ga.run_for_all_in(groups_dir, align_dir)
