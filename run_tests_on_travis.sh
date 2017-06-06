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
    mv .coverage temp;
    # pyclaw examples and I/O tests
    nosetests -v --first-pkg-wins --exclude=limiters --exclude=sharpclaw --with-coverage --cover-package=clawpack.pyclaw --include=IOTest;
    mv temp .coverage.doctest;
    coverage combine;
fi

if [[ "${TEST_PACKAGE}" == "petclaw" ]]; then
    # petclaw examples
    mpirun -n 4 nosetests -v --first-pkg-wins --exclude=limiters --exclude=sharpclaw --exclude=fileio;
    cd ../petclaw/tests;
    # petclaw I/O tests
    mpirun -n 4 nosetests -v --first-pkg-wins --include=IOTest;
fi

if [[ "${TEST_PACKAGE}" == "forestclaw" ]]; then
    # forestclaw tests (I/O only)
    cd ../forestclaw;
    nosetests --first-pkg-wins;
fi


