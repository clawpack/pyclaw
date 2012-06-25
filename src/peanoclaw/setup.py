#!/usr/bin/env python
"""PeanoClaw"""


DOCLINES = __doc__.split("\n")

import os
import sys
import warnings
import subprocess
import shutil
import re

if sys.version_info[0] < 3:
    import __builtin__ as builtins
else:
    import builtins

CLASSIFIERS = """\
Development Status :: 4 - Beta
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved
Programming Language :: C
Programming Language :: Python
Topic :: Software Development
Topic :: Scientific/Engineering
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS

"""

MAJOR               = 0
MINOR               = 1
MICRO               = 0
ISRELEASED          = False
VERSION             = '%d.%d.%d' % (MAJOR, MINOR, MICRO)

package_path='./'

version_file_path  =os.path.join(package_path,'version.py')

# Return the git revision as a string
def git_version():
    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        env = {}
        for k in ['SYSTEMROOT', 'PATH']:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v
        # LANGUAGE is used on win32
        env['LANGUAGE'] = 'C'
        env['LANG'] = 'C'
        env['LC_ALL'] = 'C'
        out = subprocess.Popen(cmd, stdout = subprocess.PIPE, env=env).communicate()[0]
        return out

    try:
        out = _minimal_ext_cmd(['git', 'rev-parse', 'HEAD'])
        GIT_REVISION = out.strip().decode('ascii')
    except OSError:
        GIT_REVISION = "Unknown"

    return GIT_REVISION

# BEFORE importing distutils, remove MANIFEST. distutils doesn't properly
# update it when the contents of directories change.
if os.path.exists('MANIFEST'):
    os.remove('MANIFEST')

def write_version_py(filename=version_file_path):
    cnt = """
# THIS FILE IS GENERATED FROM PEANOCLAW SETUP.PY
short_version = '%(version)s'
version = '%(version)s'
full_version = '%(full_version)s'
git_revision = '%(git_revision)s'
release = %(isrelease)s

if not release:
    version = full_version
"""
    # Adding the git rev number needs to be done inside
    # write_version_py(), otherwise the import of scipy.version messes
    # up the build under Python 3.
    FULLVERSION = VERSION
    if os.path.exists('.git'):
        GIT_REVISION = git_version()
    elif os.path.exists(version_file_path):
        # must be a source distribution, use existing version file
        from version import git_revision as GIT_REVISION
    else:
        GIT_REVISION = "Unknown"

    if not ISRELEASED:
        FULLVERSION += '.dev-' + GIT_REVISION[:7]

    a = open(filename, 'w')
    try:
        a.write(cnt % {'version': VERSION,
                       'full_version' : FULLVERSION,
                       'git_revision' : GIT_REVISION,
                       'isrelease': str(ISRELEASED)})
    finally:
        a.close()

def configuration(parent_package='clawpack',top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration('peanoclaw', parent_package, top_path)

    # ensure that peano library has been built
#    import subprocess
#    subprocess.check_call("make")

    # Python doesn't build the Peano library, so we treat it as a data file

    import os
    libname='libpeano-claw-2d'
    if os.path.exists(libname+'.so'):
        config.add_data_files(['',libname+'.so'])
    else:
        config.add_data_files(['',libname+'.dylib'])

    config.get_version('version.py')
    return config


def setup_package():    
    from numpy.distutils.core import setup

    old_path = os.getcwd()
    local_path = os.path.dirname(os.path.abspath(sys.argv[0]))
    src_path = local_path

    os.chdir(local_path)
    sys.path.insert(0, local_path)
    sys.path.insert(0, package_path)  # to retrieve version

    # Run build
    old_path = os.getcwd()
    os.chdir(src_path)
    sys.path.insert(0, src_path)

    # Rewrite the version file everytime
    write_version_py()

    try:
        setup(
            name = '',
            maintainer = "Peano and Clawpack Developers",
            maintainer_email = "unterweg@in.tum.de",
            description = DOCLINES[0],
            long_description = "\n".join(DOCLINES[2:]),
            url = "http://www5.in.tum.de/peano/releases/releases/downloads/PeanoClaw0_page_description.html",
            download_url = "http://www5.in.tum.de/persons/unterweg/peanoclaw/peanoclaw_src_20120618.tar.gz",
            license = 'BSD',
            classifiers=[_f for _f in CLASSIFIERS.split('\n') if _f],
            platforms = ["Linux", "Solaris", "Mac OS-X", "Unix"],
            configuration=configuration)
    except Exception as err:
        print err
    finally:
        del sys.path[0]
        os.chdir(old_path)
    return



if __name__ == '__main__':
    setup_package()
