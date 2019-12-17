from Clss.FileSys.GramAlign import GramAlign


class GaController:

    @staticmethod
    def align(filename, outfile):
        """
        Run aligning program from command line

        :param filename: filename of sequences for aligning
        :param outfile: out filename with aligning
        """
        ga = GramAlign()
        ga.run_gram_align(filename, outfile)

    @staticmethod
    def align_groups(groups_dir, align_dir):
        """
        Running aligning program for all groped sequences

        :param groups_dir: directory name with files
        :param align_dir: output directory
        """
        ga = GramAlign()
        ga.run_for_all_in(groups_dir, align_dir)
