#!/usr/bin/env bash

# Get <#>
git clone https://github.com/hashdist/hashdist
( 
    cd hashdist
    git checkout 231d3ade086c9f176e0bf19a9877d31703a50295
)
export PATH=`pwd`/hashdist/bin:$PATH
hit init-home

# a reasonable profile for building petsc4py on Travis
mv ../.stack.yaml stack.yaml

# 150 MB, shouldn't take too long to download
wget https://www.dropbox.com/s/h0d27lpft4ijpei/hashdist_profile_cache_5kbtllo3bawe.tar.gz
tar -zxf hashdist_profile_cache_5kbtllo3bawe.tar.gz -C ~/.hashdist

# This should now be very fast
set -o pipefail
hit develop -v stack.yaml 2>&1 | pv -i 15 -n | cat > log.txt
tail -n 60 log.txt
