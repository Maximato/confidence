from Clss.FileSys.RunAlign import RunAlign


class GaController:

    @staticmethod
    def align(filename, prog, outfile):
        """
        Run aligning program from command line

        :param filename: filename of sequences for aligning
        :param prog: string, name of program to run aligning
        :param outfile: out filename with aligning
        """
        ra = RunAlign()
        ra.run_align(filename, prog, outfile)

    @staticmethod
    def align_groups(groups_dir, prog, align_dir):
        """
        Running aligning program for all groped sequences

        :param groups_dir: directory name with files
        :param prog: string, name of program to run aligning
        :param align_dir: output directory
        """
        ra = RunAlign()
        ra.run_for_all_in(groups_dir, prog, align_dir)
