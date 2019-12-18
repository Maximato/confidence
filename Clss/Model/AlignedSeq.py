class AlignedSeq:
    @staticmethod
    def content_calc(seq):
        """
        Calculation of content in sequence. Content is all position in the middle of sequence.
        Start content is the firs not '-' symbol, end position is the last not  '-' symbol

        :param seq: Seq, sequence
        :return: int, start and end of content
        """
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
