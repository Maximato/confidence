from os.path import join, abspath, basename, dirname
import os


PROJECT_PATH = os.path.dirname(__file__)

# for tests
TEST_DATA_PATH = join(PROJECT_PATH, "Tests", "Data")
TEST_IN_DATA = join(TEST_DATA_PATH, "in")
TEST_OUT_DATA = join(TEST_DATA_PATH, "out")
