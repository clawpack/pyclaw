"""
Create gallery files for the documentation on the Clawpack website.

These tools assume that the examples have already been run and the plots
produced using "make .plots".  

You should use the script python/run_examples.py to do this first.
"""

import gallery as G
import os

G.claw_html_root='http://kingkong.amath.washington.edu/clawpack/trunk'
G.gallery_dir_default = os.path.join(G.clawdir_default,'doc','gallery')  
G.remake = True

G.make_all()

