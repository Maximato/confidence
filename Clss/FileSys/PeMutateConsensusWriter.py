from Clss.FileSys.check_directory import check_dir


class PeMutateConsensusWriter:
    def __init__(self, sequence, consensus):
        super()
        self.sequence = sequence
        self.consensus = consensus

    def write_in_pe_format(self, filename):
        check_dir(filename)

        with open(filename, 'w') as f:
            f.write(f"sequence={self.sequence}\n"
                    f"consensus={self.consensus}\n")
            f.write("F3_5pos=-1\nF3_3pos=-1\nF2_5pos=-1\nF2_3pos=-1\nF1c_5pos=-1\nF1c_3pos=-1\nB3_5pos=-1\nB3_3pos=-1\n"
                    "B2_5pos=-1\nB2_3pos=-1\nB1c_5pos=-1\nB1c_3pos=-1\ntarget_range_type=0\ntarget_range_from=\n"
                    "target_range_to=")
