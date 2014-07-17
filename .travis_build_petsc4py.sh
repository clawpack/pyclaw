#!/usr/bin/env bash

# Get <#>
git clone https://github.com/hashdist/hashdist
( 
    cd hashdist
    git checkout v0.3
)
export PATH=`pwd`/hashdist/bin:$PATH

# a reasonable profile for building petsc4py on Travis
cp ../.hashdist_petsc4py_stack.yaml petsc4py_stack.yaml

# This monitors the slow source-based Python/PETSc/petsc4py install
# Will be much faster with relocatable binary builds when available from <#>
set -o pipefail
hit develop -v petsc4py_stack.yaml 2>&1 | pv -i 15 -n | cat > log.txt
tail -n 60 log.txt
