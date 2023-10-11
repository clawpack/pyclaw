r"""
Module defining Pyclaw geometry objects.
"""

from __future__ import absolute_import
from __future__ import print_function
import numpy as np

import warnings
import six
from six.moves import range
from six.moves import zip
deprec_message = "'edges' has been deprecated; please use 'nodes' instead."
# ============================================================================
#  Default function definitions
# ============================================================================

# Default mapc2p functions
def identity_map_1d(x):
    return x,

def identity_map_2d(x,y):
    return x,y

def identity_map_3d(x,y,z):
    return x,y,z

identity_map={'1': identity_map_1d,
              '2': identity_map_2d,
              '3': identity_map_3d}

class Grid(object):
    r"""
    Representation of a single grid.

    :Dimension information:

        Each dimension has an associated name with it that can be accessed via
        that name such as ``grid.x.num_cells`` which would access the x dimension's
        number of cells.

    :Properties:

        If the requested property has multiple values, a list will be returned
        with the corresponding property belonging to the dimensions in order.

    :Initialization:

        Input:
         - *dimensions* - (list of :class:`Dimension`) Dimensions that are to
           be associated with this grid

        Output:
         - (:class:`grid`) Initialized grid object

    A PyClaw grid is usually constructed from a tuple of PyClaw Dimension objects:

        >>> from clawpack.pyclaw.geometry import Dimension, Grid
        >>> x = Dimension(0.,1.,10,name='x')
        >>> y = Dimension(-1.,1.,25,name='y')
        >>> grid = Grid((x,y))
        >>> print(grid)
        2-dimensional domain (x,y)
        No mapping
        Extent:  [0.0, 1.0] x [-1.0, 1.0]
        Cells:  10 x 25

    We can query various properties of the grid:

        >>> grid.num_dim
        2
        >>> grid.num_cells
        [10, 25]
        >>> grid.lower
        [0.0, -1.0]
        >>> grid.delta # Returns [dx, dy]
        [0.1, 0.08]

    A grid can be extended to higher dimensions using the add_dimension() method:

        >>> z=Dimension(-2.0,2.0,21,name='z')
        >>> grid.add_dimension(z)
        >>> grid.num_dim
        3
        >>> grid.num_cells
        [10, 25, 21]

    **Coordinates:**

    We can get the x, y, and z-coordinate arrays of cell nodes and centers from the grid.
    Properties beginning with 'c' refer to the computational (unmapped) domain, while
    properties beginning with 'p' refer to the physical (mapped) domain.  For grids with
    no mapping, the two are identical.  Also note the difference between 'center' and
    'centers'.

        >>> import numpy as np
        >>> np.set_printoptions(precision=2)  # avoid doctest issues with roundoff
        >>> grid.c_center([1,2,3])
        array([ 0.15, -0.8 , -1.33])
        >>> grid.p_nodes[0][0,0,0]
        0.0
        >>> grid.p_nodes[1][0,0,0]
        -1.0
        >>> grid.p_nodes[2][0,0,0]
        -2.0

    It's also possible to get coordinates for ghost cell arrays:

        >>> x = Dimension(0.,1.,5,name='x')
        >>> grid1d = Grid([x])
        >>> grid1d.c_centers
        [array([0.1, 0.3, 0.5, 0.7, 0.9])]
        >>> grid1d.c_centers_with_ghost(2)
        [array([-0.3, -0.1,  0.1,  0.3,  0.5,  0.7,  0.9,  1.1,  1.3])]

    **Mappings:**

    A grid mapping can be used to solve in a domain that is not rectangular,
    or to adjust the local spacing of grid cells.  For instance, we can
    use smaller cells on the left and larger cells on the right by doing:

        >>> double = lambda xarr : np.array([x**2 for x in xarr])
        >>> grid1d.mapc2p = double
        >>> grid1d.p_centers
        array([0.01, 0.09, 0.25, 0.49, 0.81])

    Note that the 'nodes' (or nodes) of the mapped grid are the mapped values
    of the computational nodes.  In general, they are not the midpoints between
    mapped centers:

        >>> grid1d.p_nodes
        array([0.  , 0.04, 0.16, 0.36, 0.64, 1.  ])
    """

    def __getattr__(self,key):
        # Provide dimension attribute lists when requested from Grid object.
        # Note that this only gets called when one requests an attribute
        # that the grid itself doesn't possess.
        if key in ['num_cells','lower','upper','delta','units','centers','nodes',
                    'on_lower_boundary','on_upper_boundary']:
            return self.get_dim_attribute(key)
        else:
            raise AttributeError("'Grid' object has no attribute '"+key+"'")

    # ========== Property Definitions ========================================
    @property
    def num_dim(self):
        r"""(int) - Number of dimensions"""
        return len(self._dimensions)
    @property
    def dimensions(self):
        r"""(list) - List of :class:`Dimension` objects defining the
                grid's extent and resolution"""
        return [getattr(self,name) for name in self._dimensions]
    @property
    def c_centers(self):
        r"""(list of ndarray(...)) - List containing the arrays locating
                  the computational locations of cell centers, see
                  :meth:`_compute_c_centers` for more info."""
        self._compute_c_centers()
        return self._c_centers
    @property
    def c_nodes(self):
        r"""(list of ndarray(...)) - List containing the arrays locating
                  the computational locations of cell nodes, see
                  :meth:`_compute_c_nodes` for more info."""
        self._compute_c_nodes()
        return self._c_nodes
    @property
    def p_centers(self):
        r"""(list of ndarray(...)) - List containing the arrays locating
                  the physical locations of cell centers, see
                  :meth:`_compute_p_centers` for more info."""
        self._compute_p_centers()
        return self._p_centers
    @property
    def p_nodes(self):
        r"""(list of ndarray(...)) - List containing the arrays locating
                  the physical locations of cell nodes, see
                  :meth:`_compute_p_nodes` for more info."""
        self._compute_p_nodes()
        return self._p_nodes
    @property
    def mapc2p(self):
        return self._mapc2p
    @mapc2p.setter
    def mapc2p(self,mapc2p):
        self._mapc2p = mapc2p
        self._clear_cached_values()


    # ========== Class Methods ===============================================
    def __init__(self,dimensions):
        r"""
        Instantiate a Grid object

        See :class:`Grid` for more info.
        """

        # ========== Attribute Definitions ===================================
        r"""(func) - Coordinate mapping function"""
        self.gauges = []
        r"""(list) - List of gauges' indices to be filled by add_gauges
        method.
        """
        self.gauge_file_names  = []
        r"""(list) - List of file names to write gauge values to"""
        self.gauge_files = []
        r"""(list) - List of file objects to write gauge values to"""
        self.gauge_dir_name = '_gauges'
        r"""(string) - Name of the output directory for gauges. If the
        `Controller` class is used to run the application, this directory by
        default will be created under the `Controller` `outdir` directory.
        """
        self._p_centers = None
        self._p_nodes = None
        self._c_centers = None
        self._c_nodes = None

        # Dimension parsing
        if isinstance(dimensions,Dimension):
            dimensions = [dimensions]
        self._dimensions = []
        for dim in dimensions:
            self.add_dimension(dim)

        super(Grid,self).__init__()

    def _clear_cached_values(self):
        self._p_centers = None
        self._p_nodes = None
        self._c_centers = None
        self._c_nodes = None

    # ========== Dimension Manipulation ======================================
    def add_dimension(self,dimension):
        r"""
        Add the specified dimension to this patch

        :Input:
         - *dimension* - (:class:`Dimension`) Dimension to be added
        """

        # Add dimension to name list and as an attribute
        if dimension.name in self._dimensions:
            raise Exception('Unable to add dimension. A dimension'+
                            ' of the same name: {name}, already exists.'
                            .format(name=dimension.name))

        self._dimensions.append(dimension.name)
        setattr(self,dimension.name,dimension)
        self._clear_cached_values()
        # Reset mapping as it presumably makes no sense now
        self.mapc2p = identity_map[str(self.num_dim)]

    def get_dim_attribute(self,attr):
        r"""
        Returns a tuple of all dimensions' attribute attr
        """
        return [getattr(dim,attr) for dim in self.dimensions]

    def __copy__(self):
        return self.__class__(self)

    def __str__(self):
        output = "%s-dimensional domain " % str(self.num_dim)
        output += "("+",".join([dim.name for dim in self.dimensions])+")\n"
        if self.mapc2p in list(identity_map.values()):
            output += "No mapping\n"
            output += "Extent:  "
        else:
            output += "Mapping function: "+self.mapc2p.__name__+"\n"
            output += "Computational domain: "
        output += " x ".join(["[{:.2}, {:.2}]".format(dim.lower, dim.upper)
                             for dim in self.dimensions])
        output += "\n"
        output += "Cells:  "
        output += " x ".join(["{}".format(dim.num_cells) for dim in self.dimensions])
        return output

    # ========== Coordinates =============================================
    def _compute_c_centers(self, recompute=False):
        r"""Calculate the coordinates of the centers in the computational domain.

        :Input:
         - *recompute* - (bool) Whether to force a recompute of the arrays
        """
        if recompute or (self._c_centers is None) or \
           any([c is None for c in self.get_dim_attribute('_centers')]):
            index = np.indices(self.num_cells)
            self._c_centers = []
            for i,center_array in enumerate(self.get_dim_attribute('centers')):
                self._c_centers.append(center_array[index[i,...]])

    def _compute_c_nodes(self, recompute=False):
        r"""Calculate the coordinates of the nodes in the computational domain.

        :Input:
         - *recompute* - (bool) Whether to force a recompute of the arrays
        """
        if recompute or (self._c_nodes is None) or \
           any([c is None for c in self.get_dim_attribute('_nodes')]):
            index = np.indices(n+1 for n in self.num_cells)
            self._c_nodes = []
            for i,edge_array in enumerate(self.get_dim_attribute('nodes')):
                self._c_nodes.append(edge_array[index[i,...]])

    def _compute_p_centers(self, recompute=False):
        r"""Calculate the coordinates of the centers in the physical domain.

        :Input:
         - *recompute* - (bool) Whether to force a recompute of the arrays
        """
        if recompute or (self._p_centers is None) or \
           any([c is None for c in self.get_dim_attribute('_centers')]):
            self._compute_c_centers(recompute=recompute)
            self._p_centers = self.mapc2p(*self._c_centers)

    def _compute_p_nodes(self, recompute=False):
        r"""Calculate the coordinates of the nodes (corners) in the physical domain.

        :Input:
         - *recompute* - (bool) Whether to force a recompute of the arrays
        """
        if recompute or (self._p_nodes is None) or \
           any([c is None for c in self.get_dim_attribute('_nodes')]):
            self._compute_c_nodes(recompute=recompute)
            self._p_nodes = self.mapc2p(*self._c_nodes)

    def c_center(self,ind):
        r"""Compute center of computational cell with index ind."""

        index = [np.array(i) for i in ind]
        return np.array([self.c_centers[i][index] for i in range(self.num_dim)])

    def p_center(self,ind):
        r"""Compute center of physical cell with index ind."""
        return self.mapc2p(*self.c_center(ind))

    def c_centers_with_ghost(self, num_ghost):
        r"""
        Calculate the coordinates of the cell centers, including
        ghost cells, in the computational domain.

        :Input:
         - *num_ghost* - (int) Number of ghost cell layers
        """
        index = np.indices(n+2*num_ghost for n in self.num_cells)
        centers = []
        for i,dim in enumerate(self.dimensions):
            center_array = dim.centers_with_ghost(num_ghost)
            centers.append(center_array[index[i,...]])
        return centers

    def c_nodes_with_ghost(self, num_ghost):
        r"""
        Calculate the coordinates of the cell nodes (corners), including
        ghost cells, in the computational domain.

        :Input:
         - *num_ghost* - (int) Number of ghost cell layers
        """
        index = np.indices(n+2*num_ghost+1 for n in self.num_cells)
        nodes = []
        for i,dim in enumerate(self.dimensions):
            edge_array = dim.nodes_with_ghost(num_ghost)
            nodes.append(edge_array[index[i,...]])
        return nodes

    def p_centers_with_ghost(self,num_ghost):
        return self.mapc2p(*self.c_centers_with_ghost(num_ghost))

    def p_nodes_with_ghost(self,num_ghost):
        return self.mapc2p(*self.c_nodes_with_ghost(num_ghost))

    # ========================================================================
    # Edges: deprecated; will be removed in 6.0
    @property
    def c_edges(self):
        warnings.warn(deprec_message)
        return self.c_nodes
    @property
    def p_edges(self):
        warnings.warn(deprec_message)
        return self.p_nodes
    def p_edges_with_ghost(self,num_ghost):
        warnings.warn(deprec_message)
        return self.p_nodes_with_ghost(num_ghost)
    def c_edges_with_ghost(self, num_ghost):
        warnings.warn(deprec_message)
        return self.c_nodes_with_ghost(num_ghost)
    # ========================================================================

    # ========================================================================
    #  Gauges
    # ========================================================================
    def add_gauges(self,gauge_coords):
        r"""
        Determine the cell indices of each gauge and make a list of all gauges
        with their cell indices.
        """
        for gauge in gauge_coords:
            # Check if gauge belongs to this grid:
            if all(self.lower[n]<=gauge[n]<self.upper[n] for n in range(self.num_dim)):
                # Set indices relative to this grid
                gauge_index = [int(round((gauge[n]-self.lower[n])/self.delta[n]))
                               for n in range(self.num_dim)]
                gauge_file_name = 'gauge'+'_'.join(str(coord) for coord in gauge)+'.txt'
                self.gauge_file_names.append(gauge_file_name)
                self.gauges.append(gauge_index)

    def setup_gauge_files(self,outdir):
        r"""
        Creates and opens file objects for gauges.
        """
        import os
        gauge_path = os.path.join(outdir,self.gauge_dir_name)
        if not os.path.exists(gauge_path):
            try:
                os.makedirs(gauge_path)
            except OSError:
                print("gauge directory already exists, ignoring")

        for gauge in self.gauge_file_names:
            gauge_file = os.path.join(gauge_path,gauge)
            if os.path.isfile(gauge_file):
                os.remove(gauge_file)
            self.gauge_files.append(open(gauge_file,'a'))

    def plot(self,num_ghost=0,mapped=True,mark_nodes=False,mark_centers=False):
        r"""Make a plot of the grid.

        By default the plot uses the mapping
        grid.mapc2p and does not show any ghost cells.  This can be modified
        via the arguments `mapped` and `num_ghost`.

        Returns a handle to the plot axis object.
        """
        import matplotlib.pyplot as plt
        if self.num_dim == 2:
            fig, ax = plt.subplots(1,1)
            if num_ghost>0:
                if mapped:
                    xe, ye = self.p_nodes_with_ghost(num_ghost)
                else:
                    xe, ye = self.c_nodes_with_ghost(num_ghost)
                p = ax.pcolormesh(xe,ye,0*xe,edgecolors='k',cmap='bwr',alpha=0.2)
                p.set_clim(-1,1)
            if mapped:
                xe, ye = self.p_nodes
                xc, yc = self.p_centers
            else:
                xe, ye = self.c_nodes
                xc, yc = self.c_centers
            p = ax.pcolormesh(xe,ye,0*xe,edgecolors='k',cmap='bwr')
            p.set_clim(-1,1)
            if mark_nodes:
                ax.plot(xe,ye,'or')
            if mark_centers:
                ax.plot(xc,yc,'ob')
            ax.axis('equal')
            ax.set_xlabel(self.dimensions[0].name)
            ax.set_ylabel(self.dimensions[1].name)
            return ax
        else:
            raise Exception('Grid plotting implemented for 2D grids only.')

    def _check_validity(self):
        for dim in self.dimensions:
            dim._check_validity()
        assert type(self.num_cells) is int, 'Dimension.num_cells must be an integer'
        assert type(self.lower) is float, 'Dimension.lower must be a float'
        assert type(self.upper) is float, 'Dimension.upper must be a float'
        assert self.num_cells>0, 'Dimension.num_cells must be positive'
        assert self.upper > self.lower, 'Dimension.upper must be greater than lower'


# ============================================================================
#  Dimension Object
# ============================================================================
class Dimension(object):
    r"""
    Basic class representing a dimension of a Patch object

    :Initialization:

    Required arguments, in order:
     - *lower* - (float) Lower extent of dimension
     - *upper* - (float) Upper extent of dimension
     - *num_cells* - (int) Number of cells

    Optional (keyword) arguments:
     - *name* - (string) string Name of dimension
     - *units* - (string) Type of units, used for informational purposes only

    Output:
     - (:class:`Dimension`) - Initialized Dimension object

    Example:

    >>> from clawpack.pyclaw.geometry import Dimension
    >>> x = Dimension(0.,1.,100,name='x')
    >>> print(x)
    Dimension x:  (num_cells,delta,[lower,upper]) = (100,0.01,[0.0,1.0])
    >>> x.name
    'x'
    >>> x.num_cells
    100
    >>> x.delta
    0.01
    >>> x.nodes[0]
    0.0
    >>> x.nodes[1]
    0.01
    >>> x.nodes[-1]
    1.0
    >>> x.centers[-1]
    0.995
    >>> len(x.centers)
    100
    >>> len(x.nodes)
    101
    """

    @property
    def delta(self):
        r"""(float) - Size of an individual, computational cell"""
        return (self.upper-self.lower) / float(self.num_cells)

    # ========== Edges: deprecated; will be removed in 6.0 =======
    @property
    def edges(self):
        warnings.warn(deprec_message)
        return self.nodes

    def edges_with_ghost(self,num_ghost):
        warnings.warn(deprec_message)
        return self.nodes_with_ghost(num_ghost)
    # ========================================================================

    # ========== Centers and nodes ========================================
    @property
    def nodes(self):
        r"""(ndarrary(:)) - Location of all cell edge coordinates
        for this dimension"""
        if self._nodes is None:
            self._nodes = np.empty(self.num_cells+1)
            for i in range(0,self.num_cells+1):
                self._nodes[i] = self.lower + i*self.delta
        return self._nodes
    @property
    def centers(self):
        r"""(ndarrary(:)) - Location of all cell center coordinates
        for this dimension"""
        if self._centers is None:
            self._centers = np.empty(self.num_cells)
            for i in range(0,self.num_cells):
                self._centers[i] = self.lower + (i+0.5)*self.delta
        return self._centers
    @property
    def lower(self):
        return self._lower
    @lower.setter
    def lower(self,lower):
        self._lower = float(lower)
        self._centers = None  # Reset cached arrays
        self._nodes = None
        self._check_validity()
    @property
    def upper(self):
        return self._upper
    @upper.setter
    def upper(self,upper):
        self._upper = float(upper)
        self._centers = None  # Reset cached arrays
        self._nodes = None
        self._check_validity()
    @property
    def num_cells(self):
        return self._num_cells
    @num_cells.setter
    def num_cells(self,num_cells):
        self._num_cells = int(num_cells)
        self._centers = None  # Reset cached arrays
        self._nodes = None
        self._check_validity()

    def centers_with_ghost(self,num_ghost):
        r"""(ndarrary(:)) - Location of all cell center coordinates
        for this dimension, including centers of ghost cells."""
        centers = self.centers
        pre = self.lower+(np.arange(-num_ghost,0)+0.5)*self.delta
        post = self.upper + self.delta * (np.arange(num_ghost) + 0.5)
        return np.hstack((pre,centers,post))

    def nodes_with_ghost(self,num_ghost):
        r"""(ndarrary(:)) - Location of all edge coordinates
        for this dimension, including nodes of ghost cells."""
        nodes   = self.nodes
        pre  = np.linspace(self.lower-num_ghost*self.delta,self.lower-self.delta,num_ghost)
        post = np.linspace(self.upper+self.delta, self.upper+num_ghost*self.delta,num_ghost)
        return np.hstack((pre,nodes,post))

    def __init__(self, lower, upper, num_cells, name='x',
                 on_lower_boundary=None,on_upper_boundary=None, units=None):
        r"""
        Create a Dimension object.

        See :class:`Dimension` for full documentation
        """
        if isinstance(lower,six.string_types):
            raise Exception('Passing dimension name as first argument is deprecated. \
                             Pass it as a keyword argument instead.')

        self._nodes = None
        self._centers = None
        self._centers_with_ghost = None
        self._nodes_with_ghost = None

        self._lower = float(lower)
        self._upper = float(upper)
        self._num_cells = int(num_cells)
        self.name = name
        self.on_lower_boundary = on_lower_boundary
        self.on_upper_boundary = on_upper_boundary
        self.units = units

        self._check_validity()

    def _check_validity(self):
        assert isinstance(self.num_cells,int), 'Dimension.num_cells must be an integer; got %s' % type(self.num_cells)
        assert isinstance(self.lower,float), 'Dimension.lower must be a float'
        assert isinstance(self.upper,float), 'Dimension.upper must be a float'
        assert self.num_cells>0, 'Dimension.num_cells must be positive'
        assert self.upper > self.lower, 'Dimension.upper must be greater than lower'

    def __str__(self):
        output = "Dimension %s" % self.name
        if self.units:
            output += " (%s)" % self.units
        output += ":  (num_cells,delta,[lower,upper]) = (%s,%s,[%s,%s])" \
            % (self.num_cells,self.delta,self.lower,self.upper)
        return output

    def __len__(self):
        return self.num_cells


# ============================================================================
#  Pyclaw Patch object definition
# ============================================================================
class Patch(object):
    """
    :Global Patch information:

        Each patch has a value for :attr:`level` and :attr:`patch_index`.
    """
    # Global properties
    @property
    def num_cells_global(self):
        r"""(list) - List of the number of cells in each dimension"""
        return self.get_dim_attribute('num_cells')
    @property
    def lower_global(self):
        r"""(list) - Lower coordinate extents of each dimension"""
        return self.get_dim_attribute('lower')
    @property
    def upper_global(self):
        r"""(list) - Upper coordinate extends of each dimension"""
        return self.get_dim_attribute('upper')
    @property
    def num_dim(self):
        r"""(int) - Number of dimensions"""
        return len(self._dimensions)
    @property
    def dimensions(self):
        r"""(list) - List of :class:`Dimension` objects defining the
                grid's extent and resolution"""
        return [getattr(self,name) for name in self._dimensions]
    @property
    def delta(self):
        r"""(list) - List of computational cell widths"""
        return self.get_dim_attribute('delta')
    @property
    def name(self):
        r"""(list) - List of names of each dimension"""
        return self._dimensions

    def __init__(self,dimensions):
        self.level = 1
        r"""(int) - AMR level this patch belongs to, ``default = 1``"""
        self.patch_index = 1
        r"""(int) - Patch number of current patch, ``default = 0``"""

        if isinstance(dimensions,Dimension):
            dimensions = [dimensions]
        self._dimensions = []
        for dim in dimensions:
            dim.on_lower_boundary = True
            dim.on_upper_boundary = True
            self.add_dimension(dim)

        self.grid = Grid(dimensions)

        super(Patch,self).__init__()

    def add_dimension(self,dimension):
        r"""
        Add the specified dimension to this patch

        :Input:
         - *dimension* - (:class:`Dimension`) Dimension to be added
        """

        # Add dimension to name list and as an attribute
        if dimension.name in self._dimensions:
            raise Exception('Unable to add dimension. A dimension'+
                            ' of the same name: {name}, already exists.'
                            .format(name=dimension.name))

        self._dimensions.append(dimension.name)
        setattr(self,dimension.name,dimension)

    def get_dim_attribute(self,attr):
        r"""
        Returns a tuple of all dimensions' attribute attr
        """
        return [getattr(getattr(self,name),attr) for name in self._dimensions]

    def __deepcopy__(self,memo={}):
        import copy
        result = self.__class__(copy.deepcopy(self.dimensions))
        result.__init__(copy.deepcopy(self.dimensions))
        result.grid.mapc2p = self.grid.mapc2p

        for attr in ('level','patch_index'):
            setattr(result,attr,copy.deepcopy(getattr(self,attr)))

        return result

    def __str__(self):
        output = "Patch %s:\n" % self.patch_index
        output += '\n'.join((str(getattr(self,dim)) for dim in self._dimensions))
        return output

# ============================================================================
#  Pyclaw Domain object definition
# ============================================================================
class Domain(object):
    r"""
    A Domain is a list of Patches.

    A Domain may be initialized in the following ways:

        1. Using 3 arguments, which are in order
            - A list of the lower boundaries in each dimension
            - A list of the upper boundaries in each dimension
            - A list of the number of cells to be used in each dimension

        2. Using a single argument, which is
            - A list of dimensions; or
            - A list of patches.

    :Examples:

        >>> from clawpack import pyclaw
        >>> domain = pyclaw.Domain( (0.,0.), (1.,1.), (100,100))
        >>> print(domain.num_dim)
        2
        >>> print(domain.grid.num_cells)
        [100, 100]
    """
    @property
    def num_dim(self):
        r"""(int) - :attr:`Patch.num_dim` of base patch"""
        return self._get_base_patch_attribute('num_dim')
    @property
    def patch(self):
        r"""(:class:`Patch`) - First patch is returned"""
        return self.patches[0]
    @property
    def grid(self):
        r"""(list) - :attr:`Patch.grid` of base patch"""
        return self._get_base_patch_attribute('grid')

    def __init__(self,*arg):
        if len(arg)>1:
            lower = arg[0]
            upper = arg[1]
            n     = arg[2]
            dims = []
            names = ['x','y','z']
            names = names[:len(n)+1]
            for low,up,nn,name in zip(lower,upper,n,names):
                dims.append(Dimension(low,up,nn,name=name))
            self.patches = [Patch(dims)]
        else:
            geom = arg[0]
            if not isinstance(geom,list) and not isinstance(geom,tuple):
                geom = [geom]
            if isinstance(geom[0],Patch):
                self.patches = geom
            elif isinstance(geom[0],Dimension):
                self.patches = [Patch(geom)]

    def _get_base_patch_attribute(self, name):
        r"""
        Return base patch attribute name

        :Output:
         - (id) - Value of attribute from ``self.patches[0]``
        """
        return getattr(self.patches[0],name)

    def __deepcopy__(self,memo={}):
        import copy
        result = self.__class__(copy.deepcopy(self.patches))
        result.__init__(copy.deepcopy(self.patches))

        return result


if __name__ == "__main__":
    import doctest
    doctest.testmod()
