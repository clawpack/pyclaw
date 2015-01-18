#!/usr/bin/env bash

# Get <#>
git clone https://github.com/hashdist/hashdist
( 
    cd hashdist
    git checkout v0.3
)
export PATH=`pwd`/hashdist/bin:$PATH
hit init-home

# a reasonable profile for building petsc4py on Travis
mv ../.hashdist_petsc4py_stack.yaml petsc4py_stack.yaml

# 75 MB, shouldn't take too long to download
wget https://dl.dropboxusercontent.com/u/65439/pyclaw_hashdist_bld.tar.gz
tar -zxf pyclaw_hashdist_bld.tar.gz -C ~/.hashdist

# This should now be very fast
set -o pipefail
hit develop -v petsc4py_stack.yaml 2>&1 | pv -i 15 -n | cat > log.txt
tail -n 60 log.txt
