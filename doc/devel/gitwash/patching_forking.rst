.. _contribution:

==========================================
Contributing to *PyClaw* code development 
==========================================

If you are interested in contributing to *PyClaw* development, you
can add your contribution through making patches; see :ref:`making-patches`.
However, for substantial contributions, you may want to consider forking
`PyClaw github`_ repository, see :ref:`forking`, and making pull requests
for your changes. Here we introduce those two methods for contributing 
to PyClaw_ source code. For more information about the forking model,
please see :ref:`development-workflow`.


.. _following-latest:

Getting *PyClaw* from the development repository
================================================

After installing git, you can get the latest *PyClaw* source from the
`PyClaw github`_ repository using the following command: ::

   $ git clone git://github.com/clawpack/pyclaw.git

You now have a copy of the code tree in the new ``pyclaw`` directory. 
From time to time you may want to pull down the latest code.  Do this with::

   cd pyclaw
   git pull

The tree in ``pyclaw`` will now have the latest changes from ``pyclaw``
repository.

.. _making-patches:

Making a patch, for quick fixes
===============================

You've discovered a bug or something else you want to change
in PyClaw_ .. |emdash| excellent!

You've worked out a way to fix it |emdash| even better!

You want to tell us about it |emdash| best of all!

The easiest way is to make a *patch* or set of patches.  Here
we explain how.  Making a patch is the simplest and quickest,
but if you're going to be doing anything more than simple
quick things, please consider :ref:`forking` instead.

Following are the steps for making patches:

#. Tell git who you are so it can label the commits you've
   made::

      git config --global user.email you@yourdomain.example.com
      git config --global user.name "Your Name Comes Here"

#. If you don't already have one, clone a copy of the
   PyClaw_ repository::

      git clone git://github.com/clawpack/pyclaw.git
      cd pyclaw

#. Make a 'feature branch'.  This will be where you work on
   your bug fix.  It's nice and safe and leaves you with
   access to an unmodified copy of the code in the main
   branch::

      git branch the-fix-im-thinking-of
      git checkout the-fix-im-thinking-of

#. Do some edits, and commit them as you go::

      # hack, hack, hack
      # Tell git about any new files you've made
      git add somewhere/tests/test_my_bug.py
      # commit work in progress as you go
      git commit -am 'BF - added tests for Funny bug'
      # hack hack, hack
      git commit -am 'BF - added fix for Funny bug'

   Note the ``-am`` options to ``commit``. The ``m`` flag just
   signals that you're going to type a message on the command
   line.  The ``a`` flag |emdash| you can just take on faith |emdash|
   or see `why the -a flag?`_.

#. When you have finished, check you have committed all your
   changes::

      git status

#. Finally, make your commits into patches.  You want all the
   commits since you branched from the ``master`` branch::

      git format-patch -M -C master

   You will now have several files named for the commits::

      0001-BF-added-tests-for-Funny-bug.patch
      0002-BF-added-fix-for-Funny-bug.patch

   Send these files to the `PyClaw mailing list`_ |emdash| where
   we will thank you warmly.


When you are done, to switch back to the main copy of the
code, just return to the ``master`` branch::

   git checkout master

.. _forking:

Making a fork, for substantial contributions
============================================

If you find you have done some patches, and you have one or
more feature branches, you will probably want to start working
on a separate fork and make pull requests.

You need to make a fork only once. The instructions here are very similar
to the instructions at http://help.github.com/forking/ |emdash| please see
that page for more detail.  We're repeating some of it here just to give the
specifics for the PyClaw_ project, and to suggest some default names.

Set up and configure a github account
-------------------------------------

If you don't have a github account, go to the github page, and make one.

You then need to configure your account to allow write access |emdash| see
the `Generating SSH keys` help on `github help`_.

Create your own forked copy of PyClaw_
--------------------------------------

#. Log into your github account.
#. Go to the PyClaw_ github home at `PyClaw github`_.
#. Click on the *fork* button:

   .. image:: forking_button.png

   Now, after a short pause and some 'Hardcore forking action', you
   should find yourself at the home page for your own forked copy of PyClaw_.

.. _set-up-fork:

.. _linking-to-upstream:

Set your local clone of PyClaw_ to point to your fork and to the main repository
----------------------------------------------------------------------------------

If you already have a local clone of the main PyClaw_ repository, you can
set the clone origin to be your fork while keeping the main repository
as an additional remote that you can pull from. 
The following are the required steps: ::
   
   # checkout and refresh master branch from main repo
   git checkout master
   git pull origin master
   # rename pointer to main repository to 'upstream'
   git remote rename origin upstream
   # point your repo to default read / write to your fork on github
   git remote add origin git@github.com:your-user-name/pyclaw.git
   # push up any branches you've made and want to keep
   git push origin the-fix-im-thinking-of

Then you can, if you want, follow the
:ref:`development-workflow`.


You also can directly clone your fork and set it to additionally 
point to the main repository. Here are the required steps:

#. Clone your fork to your local machine with :: 

    $ git clone git@github.com:your-user-name/pyclaw.git

#. To view all the branches you have ::

    $ cd pyclaw
    $ git branch -a

   You'll get something like::

    * master
    remotes/origin/master


   This tells you that you are currently on the ``master`` branch, and
   that you also have a ``remote`` connection to ``origin/master``.

#. To view the URLs of the remote repositories you have, try ::

       $ git remote -v

   which will point to your github fork in this case.

  
#.  Now you want to connect to the upstream `PyClaw github`_ repository, so
    you can merge in changes from trunk. 
    To point your clone to the main PyClaw repository, do: ::

       $ git remote add upstream git://github.com/clawpack/pyclaw.git

   ``upstream`` here is just the arbitrary name we're using to refer to the
   main PyClaw_ repository at `PyClaw github`_.

   Note that we've used ``git://`` for the URL rather than ``git@``.  The
   ``git://`` URL is read only.  This means that we can't accidentally
   (or deliberately) write to the upstream repository, and we are 
   only going to use it to merge into our own code.

#. Just for your own satisfaction, show yourself that you now have a new
   'remote', with ``git remote -v show``, giving you something like::

       upstream	git://github.com/clawpack/pyclaw.git (fetch)
       upstream	git://github.com/clawpack/pyclaw.git (push)
       origin	git@github.com:your-user-name/pyclaw.git (fetch)
       origin	git@github.com:your-user-name/pyclaw.git (push)

.. include:: links.inc
