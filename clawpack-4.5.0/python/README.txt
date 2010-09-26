
begin_html [use: doc/doc.css]

<h2>Python scripts and modules for Clawpack</h2>

<b>Directory claw/python</b>

These need to be cleaned up, reorganized, and better documented...


<ul>
<li> [code: startserver.py] Use this to start a web server at port 50005.
<p>
<li> [code: make_clean.py] Execute "make clean" in all subdirectories of the
specified directory.  Useful to clean up executables, .o and .html files.
<p>
<li> [code: make_clobber.py] Execute "make clobber" in all subdirectories of the
specified directory.  This does "make clean" but also removes all output and
plot directories.
<p>
<li> [code: make_htmls.py] Creates .html files in all directories where
"make .htmls" can be used.
<p>
<li> [code: make_plots.py] Runs all examples and creates plots files by
doing "make .plots" whereever possible.  If a file make_all.sh exists then
this is executed rather than "make .plots", useful when other things need to
be done as well.
<p>
<li> [code: make_gallery.py] Create galleries showing thumbnails for
applications.
<p>
<li> [code: gallery.py] Module used by make_gallery.py.
<p>
<li> [code: convert43to44.py] Convert an application directory that works
with Clawpack 4.3 to a form that works with 4.4 or 4.5.  May not work
perfectly, but does much of the work.  Does not yet work in amr directories.
<p>

</ul>

end_html

