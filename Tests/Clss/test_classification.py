import unittest
from Clss.Model.Classification import *


class TestAligns(unittest.TestCase):
    def test_get_ccls(self):
        self.assertEqual("c10", get_ccls(0.11))

    def test_get_dcls(self):
        self.assertEqual("c50", get_dcls(555))
