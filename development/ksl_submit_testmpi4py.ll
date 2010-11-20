#!/usr/bin/env bash
#
# @ job_name            = test_mpi4py 
# @ job_type            = bluegene
# @ output              = ./$(job_name)_$(jobid).out
# @ error               = ./$(job_name)_$(jobid).err
# @ environment         = COPY_ALL; 
# @ wall_clock_limit    = 00:10:00,00:10:00
# @ notification        = always
# @ bg_size             = 1 
# @ cluster_list        = bgp
# @ class               = default
# @ account_no          = k47

# @ queue
/bgsys/drivers/ppcfloor/bin/mpirun -exp_env LD_LIBRARY_PATH -env PYTHONPATH="${PYTHONPATH_PPC450D}" -env BG_MAPPING=TXYZ -np 4 -mode VN /bgsys/drivers/ppcfloor/gnu-linux/bin/python ./test_mpi4py.py

