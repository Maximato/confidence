from unittest import TestCase
from Clss.Model.HtmlConsensusParser import HtmlConsensusParser


class TestHtmlConsensusParser(TestCase):
    def test_get_consensus_with_mut(self):
        # setup data
        html_consensus = "<span class='c90'>G</sapan><span class='c60'>C</sapan>"

        # test function
        html_consensus_parser = HtmlConsensusParser()
        html_consensus_parser.parse_html_consensus(html_consensus, "c90 c80")
        consensus_string = html_consensus_parser.consensus_string
        confidence_string = html_consensus_parser.confidence_string

        # assertions
        self.assertEqual("GC", consensus_string)
        self.assertEqual("*-", confidence_string)
