from os.path import dirname
import os


def check_dir(filename):
    """
    Supporting function
    :param filename: filename to check
    """

    # get directory of filename
    dname = dirname(filename)

    # creating output directory if not created
    if not os.path.isdir(dname) and dname != "":
        os.mkdir(dname)
