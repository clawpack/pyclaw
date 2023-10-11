#!/usr/bin/env bash
#
# @ job_name            = stegoton_shaheen
# @ job_type            = bluegene
# @ output              = ./$(job_name)_$(jobid).out
# @ error               = ./$(job_name)_$(jobid).err
# @ environment         = COPY_ALL; 
# @ wall_clock_limit    = 48:00:00,48:00:00
# @ notification        = always
# @ bg_size             = 1024
# @ cluster_list        = bgp
# @ class               = default
# @ account_no          = k47

# @ queue
/bgsys/drivers/ppcfloor/bin/mpirun -exp_env LD_LIBRARY_PATH -env PYTHONPATH="${PYTHONPATH_PPC450D}" -env BG_MAPPING=TXYZ -np 4096 -mode VN /bgsys/drivers/ppcfloor/gnu-linux/bin/python ./stegoton.py

