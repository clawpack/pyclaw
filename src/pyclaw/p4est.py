
import ctypes

# Dynamically link in the pyclaw p4est interface
libp4est = ctypes.CDLL ("pyclaw_p4est.so")

# Initialize MPI. TODO: Do this in a generic routine
libp4est.pyclaw_MPI_Init ()

# Both 2D and 3D p4est objects default to the unitcube
# Create a 2D p4est internal state
p4est = libp4est.pyclaw_p4est_new ()
# Create a 3D p4est internal state
# p8est = libp4est.pyclaw_p8est_new ()

libp4est.pyclaw_p4est_destroy (p4est)

# Finalize MPI. TODO: Do this in a generic routine
libp4est.pyclaw_MPI_Finalize ()
