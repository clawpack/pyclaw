#!/bin/bash

# Script to run tests on Travis.  This is separate from the file
# .travis.yml so that we can have if statements but make sure
# that failures at intermediate steps are caught (by setting -ev).
# See http://steven.casagrande.io/articles/travis-ci-and-if-statements/

set -ev

cd src/pyclaw
if [ "${TEST_PACKAGE}" == "pyclaw" ]; then
    # pyclaw doc-tests
    nosetests --first-pkg-wins --with-doctest --exclude=limiters --exclude=sharpclaw --exclude=fileio --exclude=example --with-coverage --cover-package=clawpack.pyclaw;
    mv .coverage .coverage.doctest;
    # pyclaw examples and I/O tests
    nosetests -v --first-pkg-wins --exclude=limiters --exclude=sharpclaw --with-coverage --cover-package=clawpack.pyclaw --include=IOTest;
    cd ../../examples;
    nosetests -v --with-coverage --cover-package=clawpack.pyclaw;
    mv .coverage ../src/pyclaw/.coverage.examples;
    cd ../src/pyclaw;
    coverage combine;
fi

if [[ "${TEST_PACKAGE}" == "petclaw" ]]; then
    # petclaw I/O tests
    cd ../petclaw/tests;
    mpirun -n 4 nosetests -v --first-pkg-wins --include=IOTest;
    # petclaw examples
    cd ../../../examples;
    mpirun -n 4 nosetests -v --first-pkg-wins;
fi

if [[ "${TEST_PACKAGE}" == "forestclaw" ]]; then
    # forestclaw tests (I/O only)
    cd ../forestclaw;
    nosetests --first-pkg-wins;
fi


