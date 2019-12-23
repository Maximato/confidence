from os.path import join, abspath, basename, dirname
import os


PROJECT_PATH = os.path.dirname(__file__)

# directory for multiply alignment programs
ga_path = "~/GRAMALIGN/src/GramAlign"
muscle_path = None
clustalo_path = None

# for tests
TEST_DATA_PATH = join(PROJECT_PATH, "Tests", "Data")
TEST_IN_DATA = join(TEST_DATA_PATH, "in")
TEST_OUT_DATA = join(TEST_DATA_PATH, "out")
