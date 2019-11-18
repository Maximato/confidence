class Consensus:
    def __init__(self, counts):
        super()
        self._counts = counts
        self._n = len(counts["A"])

    @staticmethod
    def get_group(c):
        if c == 0:
            return "c00"
        elif 0 < c < 0.1:
            return "c01"
        elif 0.1 <= c < 0.2:
            return "c10"
        elif 0.2 <= c < 0.3:
            return "c20"
        elif 0.3 <= c < 0.4:
            return "c30"
        elif 0.4 <= c < 0.5:
            return "c40"
        elif 0.5 <= c < 0.6:
            return "c50"
        elif 0.6 <= c < 0.7:
            return "c60"
        elif 0.7 <= c < 0.8:
            return "c70"
        elif 0.8 <= c < 0.9:
            return "c80"
        elif 0.9 <= c:
            return "c90"

    def confidence_calc(self):
        # calculating of confidence from counts
        sdcg = {
            "symbols": [],
            "deeps": [],
            "confidences": [],
            "groups": []
        }

        for i in range(self._n):
            summ = 0
            max_score = 0
            symbol = "-"
            for key in self._counts:
                count = self._counts[key][i]
                summ += count
                if count > max_score:
                    max_score = count
                    symbol = key

            if summ == 0:
                confidence = 1
            else:
                confidence = max_score / summ

            sdcg["symbols"].append(symbol)
            sdcg["deeps"].append(summ)
            sdcg["confidences"].append(confidence)
            sdcg["groups"].append(self.get_group(confidence))
        return sdcg
