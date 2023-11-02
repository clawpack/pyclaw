#!/usr/bin/env python

"""Module contains class definitions related to dealing with gauge data"""

import os
import sys
import numpy

# See if pandas is available
try:
    import pandas
    pandas_available = True
except ImportError as e:
    pandas_available = False


class GaugeSolution(object):
    r"""Object representing a numerical gauge station

    The `GaugeSolution` class provides a generic means of reading and
    interacting with gauge data output from a run.  A header describing the
    contents of the data is read in and

    :Initialization:

        There are a few ways to create a `GaugeSolution`:
          1. Provide a `gauge_id` and a `path` to find that gauge.  This will
             initialize a gauge with the data at `path` for the gauge with gauge
             id `gauge_id`.  Note that `path` will default to the current
             working directory provided by `os.getcwd()` if `path` is not
             provided.
          2. Create an empty gauge by providing neither `gauge_id` or `path`.

        Additionally the parameter `use_pandas` can be provided which determines
        how the gauge's data is stored.  By default this is `False` and the
        `GaugeSolution` will use a `numpy` ndarray.
    """


    def __init__(self, gauge_id=None, path=None, use_pandas=False):
        r"""GaugeSolution Initiatlization Routine

        See :class:`GaugeSolution` for more info.
        """

        # Gauge descriptors
        self.id = None
        r"""(int) - Gauge id of this gauge.  Also known as gauge number"""
        self.location = None
        r"""(tuple) - Location of the gauge recorded in the header"""

        # Gauge data
        self.level = None
        r"""(ndarray(:) - int) - The level the given observation used"""
        self.t = None
        r"""(ndarray(:) - float) - The time the given observation was made at"""
        self.q = None
        r"""(ndarray(*, :) - float) - Observed data"""
        self.aux = None
        r"""(ndarray(*, :) - float) - Auxiliary data"""
        self.gtype = 'stationary'
        r"""(str) - 'stationary' or 'lagrangian'"""
        self.particle_path = None
        r"""(ndarray(:,3) - for lagrangian gauge columns are t,x,y,
            for stationary gauge, single row with [t0,x,y]
            Only implemented in 2d so far."""

        # Read in gauge data from file
        if gauge_id is not None:
            if path is None:
                path = os.getcwd()
            self.read(gauge_id, path, use_pandas)


    def read(self, gauge_id, path=None, use_pandas=False):
        r"""Read the gauge file at path into this object

        Read in the gauge with gauge id `gauge_id` located at `path`.  

        If `use_pandas` is `True` then `q` will be a Pandas frame.
            (not yet implemented)

        New in v5.9.0: If gauge00000N.txt file contains only a header,
        read data from corresponding .bin file.

        :Input:
         - *gauge_id* - (int) Gauge id to be read
         - *path* - (path) Path to directory containing the gauge data.
           Defaults to `path = os.getcwd()`.
         - *use_pandas* - (bool) If True then `q` will be a Pandas frame,
           otherwise it will be a `numpy.ndarary`.  Default to `False`.

        :Output:
         None
        """

        # Construct path to gauge
        if path is None:
            path = os.getcwd()

        # First read header from .txt file:
        gauge_file_name = "gauge%s.txt" % str(gauge_id).zfill(5)
        gauge_path = os.path.join(path, gauge_file_name)
        if not os.path.isfile(gauge_path):
            print('Did not find %s' % gauge_path)
            file_format = 'binary'

        with open(gauge_path, 'r') as gauge_file:
            # First line
            data = gauge_file.readline().split()
            self.id = int(data[2])
            if len(data) == 8:
                # 1d
                self.location = (float(data[4]),)
                num_eqn = int(data[7])
            if len(data) == 9:
                # 2d
                self.location = (float(data[4]), float(data[5]))
                num_eqn = int(data[8])
            elif len(data) == 10:
                # 3d
                self.location = (float(data[4]), float(data[5]), float(data[6]))
                num_eqn = int(data[9])

            line = gauge_file.readline()
            if 'lagrangian' in line.lower():
                # Lagrangian gauge (particle)
                self.gtype = 'lagrangian'
                line = gauge_file.readline() # third line of header
            elif 'stationary' in line.lower():
                # Standard stationary gauge
                self.gtype = 'stationary'
                line = gauge_file.readline() # third line of header
            else:
                # backward compatibility
                self.gtype = 'stationary'

            #data = line.split()
            #nvals = 1 + len(data)-6  # should agree with num_eqn
            
            # Check to see if there is also data in the .txt file,
            # otherwise perhaps it's in a binary .bin file.
            
            line = gauge_file.readline()
            if len(line) > 0:
                if 'binary32' in line:
                    file_format = 'binary32'
                elif 'binary' in line:
                    file_format = 'binary'  # also allows 'binary64'
                else:
                    # if 'ascii' in line, or no file_format line (data follows)
                    # (for backward compatibility)
                    file_format = 'ascii'
            else:
                # if file_format line missing and no data lines, try binary:
                file_format = 'binary'

        # Check to see if the gauge file name ID and that inside of the gauge
        # file are the same
        gauge_file_name = os.path.basename(path)
        if self.id != gauge_id:
            raise ValueError("Gauge ID requested does not match ID inside ",
                             "file!")

        # Read gauge data

        if use_pandas:
            raise NotImplementedError("Pandas data backend not implemented yet.")
            if not pandas_available:
                raise ImportError("Pandas not available.")


        if file_format == 'ascii':
            # data follows header in .txt file:
            data = numpy.loadtxt(gauge_path, comments="#")
            if data.ndim == 1:
                # only one line in gauge file, expand to 2d array
                data = data.reshape((1,len(data)))

        if file_format[:6] == 'binary':
            # data is in separate .bin file:
            gauge_file_name = "gauge%s.bin" % str(gauge_id).zfill(5)
            gauge_path = os.path.join(path, gauge_file_name)
            if not os.path.isfile(gauge_path):
                msg = 'No data in .txt file and did not find ' \
                         + '\n   binary file %s' %  gauge_path
                import warnings
                warnings.warn(msg)
                return

            if file_format in ['binary','binary64']:
                data = numpy.fromfile(gauge_path, dtype=numpy.float64)
            elif file_format == 'binary32':
                data = numpy.fromfile(gauge_path, dtype=numpy.float32)

            # assume rows are: level, t, q[0:num_eqn]
            nrows = 2 + num_eqn
            assert numpy.mod(len(data),nrows) == 0, \
                  '*** unexpected number of values in gauge file' \
                  + '\n*** expected nrows = %i rows'  % nrows
            ncols = int(len(data)/nrows)
            data = data.reshape((nrows,ncols), order='F').T

        self.level = data[:, 0].astype(numpy.int64)
        self.t = data[:, 1]
        self.q = data[:, 2:].transpose()

    
        if num_eqn != self.q.shape[0]:
            raise ValueError("Number of fields in gauge file does not match",
                             "recorded number in header.")

        if len(self.location) == 2:
            # lagrangian gauges only implemented in 2d so far
            if self.gtype == 'lagrangian':
                # for lagrangian gauges x,y paths are stored in q in place of u,v
                # note: only implemented in 2d geoclaw for now!
                self.particle_path = numpy.vstack((self.t, self.q[1,:], 
                                                  self.q[2,:])).T
            else:
                # set the lagrangian path to a single fixed location at t=t0:
                self.particle_path = numpy.array([[self.t[0], self.location[0], 
                                             self.location[1]]])



    def write(self, path=None, format="%+.15e"):
        r"""Write the data from this gauge to a file in `path`

        :Input:
         - *path* (path) Path to write the gauge file to.  Defaults to
           `path = os.getcwd()`.
         - *format* (str) Format string used for the field values.

        :Output:
         None
        """

        if path is None:
            path = os.getcwd()

        if not self.is_valid():
            raise ValueError("Gauge is not initialized properly.")

        gauge_file_name = "gauge%s.txt" % str(self.id).zfill(5)
        with open(os.path.join(path, gauge_file_name), "w") as gauge_file:

            gauge_file.write("# gauge_id= %s location=( %s %s ) num_eqn= %s\n" %
                 (self.id, self.location[0], self.location[1], self.q.shape[0]))
            gauge_file.write("# Columns: level time q(1 ... num_eqn)\n")

            # print(self.q.shape)
            for i in range(self.q.shape[1]):
                gauge_file.write("%02i %+.15e " % (self.level[i], self.t[i]))
                gauge_file.write(" ".join([format % value for value in self.q[:, i]]))
                gauge_file.write("\n")


    def is_valid(self):
        r"""Check to see if the gauge data has all been set.

        Simply checks to see if all attributes are not None.

        :Output:
         - (bool) Whether the gauge data has been set.
         """

        if ((self.id is not None) and (self.location is not None) and
            (self.level is not None) and (self.t is not None) and
            (self.q is not None)):

            return True

        return False


    def __repr__(self):
        if self.is_valid():
            output = "%4i" % self.id
            for j in range(len(self.location)):
                output = " ".join((output,"%17.10e" % self.location[j]))
            output = " ".join((output,"%13.6e" % self.t[0]))
            output = " ".join((output,"%13.6e\n" % self.t[-1]))
        else:
            output = None

        return output


    def __str__(self):
        return ("Gauge %s: location = %s, t = [%s, %s], %s" %
                                    (self.id, self.location,
                                     self.t[0], self.t[-1], self.gtype))


# ==============================
#  Utility Functions for Gauges
# ==============================
def compare_gauges(paths, gauge_id, fields='all'):
    r"""Make plots comparing gauges in different directories

    :Input:
     - *paths* (list) List of paths of length 2 pointing to the directories that
       the gauge files are to be taken from.
     - *gauge_id* (int) Gauge id to compare.
     - *fields* (int or list) Fields to be plotted.  If fields == 'all' then all
       available fields will be plotted.  Default is 'all'.

    :Output:
     - (matplotlib.figure.Figure) Figure object created by comparison
    """

    import matplotlib.pyplot as plt

    if len(paths) != 2:
        raise ValueError("Provide two paths to gauge files for comparison.")

    gauges = []
    for path in paths:
        gauges.append(GaugeSolution(path=path, gauge_id=gauge_id))

    if isinstance(fields, str):
        if fields.lower() == 'all':
            fields = range(gauges[0].q.shape[0])

    fig = plt.figure()
    fig.suptitle("Gauge %s" % gauges[0].id)
    for (i, n) in enumerate(fields):
        axes = fig.add_subplot(len(fields), 2, 2 * i + 1, )
        axes.plot(gauges[0].t, gauges[0].q[n, :], 'ko', label="%s" % paths[0])
        axes.plot(gauges[1].t, gauges[1].q[n, :], 'rx', label="%s" % paths[1])
        axes.set_xlabel("t")
        axes.set_ylabel("q[%s, :]" % n)
        axes.legend()

        axes = fig.add_subplot(len(fields), 2, 2 * i + 2)
        axes.plot(gauges[0].t, 
                  numpy.abs(gauges[0].q[n, :] - gauges[1].q[n, :]), 'r')
        axes.set_xlabel("t")
        axes.set_ylabel("$|q_{old}[%s, :] - q_{new}[%s, :]|$" % (n, n))

    return fig



# =============================================================
#  Utility Functions to Help Transition to the New Gauge Files
# =============================================================
def convert_gauges(path, output_path=None):
    r"""Convert gauge output data from fort.gauge file to new format

    This is provided for data files that were in the old format and will split
    up the data in the old `fort.gauge` files into separate files for each
    gauge.

    :Input:
     - *path* - (path) Path to old file.
     - *output_path* - (path) Path to directory the new gauge files will be
       placed.  Defaults to `output_path = os.getcwd()`.
    """

    if output_path is None:
        output_path = os.getcwd()

    # Load old data
    data = numpy.loadtxt(path)
    old_ids = numpy.asarray(data[:, 0], dtype=int)
    unique_ids = numpy.asarray(list(set(old_ids)))

    # Create new gauges and compare
    for gauge_id in unique_ids:
        gauge_indices = numpy.nonzero(old_ids == gauge_id)[0]
        new_gauge = GaugeSolution()
        new_gauge.id = gauge_id
        new_gauge.location = (numpy.nan, numpy.nan)
        new_gauge.level = numpy.asarray(data[gauge_indices, 1], dtype=int)
        new_gauge.t = data[gauge_indices, 2]
        new_gauge.q = data[gauge_indices, 3:].transpose()
        new_gauge.write(output_path)

        # Test gauge data
        if not compare_gauges(path, output_path, gauge_id):
            print("WARNING:  New gauge data does not match old gauge data.")
            print("          gauge id = %s" % gauge_id)



def compare_old_gauges(old_path, new_path, gauge_id, plot=False, abs_tol=1e-14,
                                                     rel_tol=0.0, 
                                                     verbose=False):
    r"""Compare old gauge data at `path` to new gauge data at same path

    Provided as a quick check to see if the function `convert_gauges` has
    correctly converted the data located in `fort.gauge`.  Really meant as a
    simple check but can be used to compare new test data to old.

    :Input:
     - *old_path* - (path) Path to old file.
     - *new_path* - (path) Path to new gauge files (directory).
     - *gauge_id* - (int) Gauge id to compare.
     - *plot* - (bool) Whether to plot the differences between the gauges.
     - *old_file_name* - (string) Name of the old gauge file, defaults to
       `old_file_name = "fort.gauge"`
     - *abs_tol* - (float) Absolute tolerance used to make the comparison.  See
       `numpy.allclose` for more info.
     - *rel_tol* - (float) Relative tolerance used to make the comparison.  See
       `numpy.allclose` for more info.

    :Output:
     - (bool) Whether the gauges agreed to double precision.  Uses
       `numpy.testing.assert_allequal` to check with the `abs_tol` and `rel_tol`
       specified above.
    """

    # Load old gauge data
    data = numpy.loadtxt(old_path)
    old_ids = numpy.asarray(data[:, 0], dtype=int)
    gauge_indices = numpy.nonzero(old_ids == gauge_id)[0]
    q = data[gauge_indices, 3:]

    # Load new data
    gauge = GaugeSolution(gauge_id, new_path)

    # Turn these into assertions or logicals
    if verbose:
        print("Comparison of gauge %s:" % gauge_id)
        print(r"           ||\Delta q||_2 = ",
                              numpy.linalg.norm(q - gauge.q.transpose(), ord=2))
        print(r"  arg(||\Delta q||_\infty = ",
                              numpy.argmax(q - gauge.q.transpose()))

    if plot:
        import matplotlib.pyplot as plt
        fig = plt.figure()
        for i in range(gauge.q.shape[0]):
            axes = fig.add_subplot(1, 3, i + 1)
            axes.plot(q[:, i] - gauge.q[i, :])
            axes.set_title("q[%s, :] comparison" % i)
            axes.set_xlabel("t (s)")
            axes.set_ylabel("q(%s, :)" % i)

    return numpy.allclose(q, gauge.q.transpose(), rtol=rel_tol, atol=abs_tol)


def check_old_gauge_data(path, gauge_id, new_gauge_path="./regression_data"):
    """Compare old gauge data to new gauge format

    Function is meant to check discrepancies between versions of the gauge
    files.  Note that this function also directly prints out some info.

    :Input:
     - *path* (string) - Path to old gauge data file
     - *gauge_id* (int) - Gauge id to compare
     - *new_gauge_path* (path) - Path to directory containing new gauge files, 
       defaults to './regression_data'.

    :Output:
     - (figure) Returns a matplotlib figure object plotting the differences in 
       time.
    """

    import matplotlib.pyplot as plt

    # Load old gauge data
    data = numpy.loadtxt(path)
    old_ids = numpy.asarray(data[:, 0], dtype=int)
    gauge_indices = numpy.nonzero(old_ids == gauge_id)[0]
    q = data[gauge_indices, 3:]

    # Load new data
    gauge = GaugeSolution(gauge_id, new_gauge_path)

    print(numpy.linalg.norm(q - gauge.q.transpose(), ord=2))
    print(numpy.argmax(q - gauge.q.transpose()))

    fig = plt.figure()
    for i in range(gauge.q.shape[0]):
        axes = fig.add_subplot(1, gauge.q.shape[0], i + 1)
        axes.plot(q[:, i] - gauge.q[i, :])
        axes.set_title("q[%s, :] comparison" % i)

    return fig


if __name__ == "__main__":

    import matplotlib.pyplot as plt

    help_msg = \
"""gauges.py path1 path2 [gauge_id] [fields...]

Plots a comparison between the gauges at path1 and path2 with gauge_id and the 
fields specified.  Only one gauge_id can be specified at a time but a number of 
fields can be specified including 'all'.
"""

    fields = 'all'
    gauge_id = 1
    if len(sys.argv) < 3:
        print(help_msg)
        sys.exit(0)
    elif len(sys.argv) >= 3:
        paths = [str(sys.argv[1]), str(sys.argv[2])]
        if len(sys.argv) > 3:
            gauge_id = int(sys.argv[3])
            if len(sys.argv) > 4:
                if sys.argv[4].lower() == 'all':
                    fields = 'all'
                else:
                    fields = [int(field for field in sys.argv[4:])]


    fig = compare_gauges(paths, gauge_id, fields)
    plt.show()
