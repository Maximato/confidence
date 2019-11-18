import os
from GramAlign import GramAlign


class RunnerGramAlign:
    @staticmethod
    def run_all(directory):
        ga = GramAlign()
        files = os.listdir(directory)
        for file in files:
            ga.run_gram_align(file)


RunnerGramAlign.run_all("Groups")
