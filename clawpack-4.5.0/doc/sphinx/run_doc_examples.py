"""
Run all the examples in the $CLAW/doc/sphinx subdirectories that are used in the
documentation webpages.
"""

import os,sys

example_dirs = ["example-acoustics-1d", \
                "example-acoustics-2d/1drad", \
                "example-acoustics-2d", \
                "example-acoustics-2d-amr"]

rootdir = os.getcwd()

for dir in example_dirs:
    os.chdir(dir)
    os.system("rm -f .output .rst .htmls")
    os.system("make .output")
    os.system("make .rst")
    os.system("make .htmls")
    os.chdir(rootdir)

os.system("chmod -R og+rX *")
os.system("touch plotexamples*.rst")
