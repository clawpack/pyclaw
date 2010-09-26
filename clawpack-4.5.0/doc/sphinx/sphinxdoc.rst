
.. _sphinxdoc:

******************************************
Compiling the Sphinx documentation locally
******************************************

For most users, the best way to view the documentation is 
`online <http://kingkong.amath.washington.edu/clawpack/www/users>`_.

This should be up to date with the current revision in the 
`download directory
<http://kingkong.amath.washington.edu/clawpack/www/clawdownload>`_.

The source files that create this documentation are included in the
distribution in `$CLAW/doc/sphinx <claw/doc/sphinx>`_.  

For those who want
to build it locally in order to view the documentation when offline,
or for those writing better documentation, try the following::

  $ python run_doc_examples.py         # to create plots for some examples
  $ python $CLAW/python/make_htmls.py  # to make html files

  $ cd $CLAW/doc/sphinx
  $ ln -s $CLAW user/claw    # symbolic link so docs can find local files
  $ make sphinx

and then, with the Python web server running (see :ref:`startserver`)
navigate to http://localhost:50005/doc/sphinx/users/.

The documentation points to some of the $CLAW/apps and $CLAW/books directories. 
To have all the plots available in those directories, you would need to also do::

  $ python $CLAW/python/make_plots.py  

jsMath should  work properly with local documentation if you view it via the
a webserver started by following the instruction in  :ref:`startserver`.
If  you plan to post documentation elsewhere for wider viewing, you need to make
sure the parameters `root` in `$CLAW/doc/load.js` and `jsMathscript` in
`$CLAW/doc/clawcode2html.py` are set properly before running 
`$CLAW/python/make_htmls.py`.
