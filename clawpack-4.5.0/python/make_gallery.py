"""
Create gallery files for local viewing.

These tools assume that the examples have already been run and the plots
produced using "make .plots".  

You should use the script python/run_examples.py to do this first.
"""

import gallery

gallery.make_all()

# Instead you can just make smaller galleries, e.g.
# gallery.make_1d()
