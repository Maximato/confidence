from os.path import dirname
import os


def check_dir(filename):
    dname = dirname(filename)
    # creating output directory
    if not os.path.isdir(dname) and dname != "":
        os.mkdir(dname)
