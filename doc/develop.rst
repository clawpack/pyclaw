.. _develop:

============================
Information for developers
============================

The PyClaw repository is hosted on Github at 
http://github.com/clawpack/pyclaw.  


Branching development model
============================
PyClaw development follows the model outlined at 
http://nvie.com/posts/a-successful-git-branching-model/, with
one important difference.  The *master* branch is used for development.
New branches are created for releases.
A nice cheat-sheet is available at
http://www.google.com/url?sa=D&q=http://www.globallinuxsecurity.pro/static/git/cheetsheet.pdf.

Contributions in the form of pull requests are also welcome; this approach
is probably the most convenient for occasional contributors.

Developer communication
============================

At the moment, developer communication takes place on the following
google groups:

  * http://groups.google.com/group/petclaw-dev/

  * http://groups.google.com/group/claw-dev/

Dependencies
============================

In general, additional dependencies on Python packages or other
software should be avoided.  If you wish to make a change that
will introduce a new dependency (including depending on a more
recent version of a particular package), it should be discussed
on the Google group first.

New versions of existing dependencies will typically be adopted 
either when new functionality provides an important benefit for
PyClaw or when the currently supported version is deemed to be
substantially outdated.

Committing
============================
Always run the tests before committing.

Be verbose in your commit messages.

It's helpful to always work in a named branch when
developing a new feature, and to merge with the --no-ff
option so that the history shows distinctly the development
of the feature.
