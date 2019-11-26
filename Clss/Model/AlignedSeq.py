class AlignedSeq:
    @staticmethod
    def content_calc(seq):
        # calculate start
        n = len(seq)
        i = 0
        while seq[i] == "-":
            i += 1

        # calculate end
        j = n - 1
        while seq[j] == "-":
            j -= 1
        return i, j

    def __init__(self, seq):
        self.seq = seq
        self.n = len(seq)
        self.start_content, self.end_content = self.content_calc(seq)
