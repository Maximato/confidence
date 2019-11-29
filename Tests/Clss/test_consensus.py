from Clss.Model.Consensus import Consensus
import unittest


class TestConsensus(unittest.TestCase):

    def test_get_str_consensus(self):
        # setup data
        consensus = Consensus({"symbols": ["A", "C"], "deeps": [101, 5], "confidences": [0.62, 0.91],
                               "ccls": ["c60", "c90"], "dcls": ["c10", "c01"]})

        # test function
        s = Consensus.get_str_consensus(consensus)

        # assertion
        self.assertEqual("AC", s)


"""    def test_get_html_body(self):
        # setup data
        consensus = {"symbols": ["A", "C"], "deeps": [101, 5], "confidences": [0.62, 0.91],
                     "ccls": ["c60", "c90"], "dcls": ["c10", "c01"]}

        # test function
        s = Consensus.get_html_body(consensus)

        # assertion
        expected = '<span class="c60">A</span><span class="c90">C</span>\n</body>\n</html>'
        self.assertEqual(expected, s)"""
