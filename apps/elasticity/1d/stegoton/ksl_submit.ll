#!/usr/bin/env bash
#
# @ job_name            = stegoton_shaheen
# @ job_type            = bluegene
# @ output              = ./$(job_name)_$(jobid).out
# @ error               = ./$(job_name)_$(jobid).err
# @ environment         = COPY_ALL; 
# @ wall_clock_limit    = 1:00:00,1:00:00
# @ notification        = always
# @ bg_size             = 512
# @ cluster_list        = bgp
# @ class               = default
# @ account_no          = k47

# @ queue
obgsys/drivers/ppcfloor/bin/mpirun -exp_env LD_LIBRARY_PATH -env PYTHONPATH="${PYTHONPATH_PPC450D}" -env BG_MAPPING=TXYZ -np 2048 -mode VN /bgsys/drivers/ppcfloor/gnu-linux/bin/python ./stegoton_shaheen.py

