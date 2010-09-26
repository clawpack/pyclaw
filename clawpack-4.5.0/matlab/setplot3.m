%  SETPLOT3 sets user defined plotting parameters
%
%      User defined Matlab script for setting various Clawpack plotting
%      parameters.  This script is called by PLOTCLAW3.  A default
%      version of this script can be found in claw/matlab/setplot3.m and
%      copied to users working directory and modifed to set things up
%      differently.
%
%      Parameters that can be set with SETPLOT3
%
%        OutputFlag        - set to 'ascii' (default) to read ascii output
%                                from files fort.qXXXX where XXXX is Frame #
%                            or to 'hdf' to read hdf output files.
%                                from files fort.qXXXX.hdf
%                            set to 'aux' to read aux array values from
%                                fort.aXXXX instead of q from fort.qXXXX
%                                These files are created only if
%                                outaux = .true.  in out3.f
%
%        OutputDir         - set to '.' to read data from current working
%                                directory (default) or set to path
%                                to output directory.
%
%        PlotType          - type of plot to produce:
% 			     - 1 = pcolor on slices (with optional contours
% 			         and isosurfaces)
% 			     - 2 = contour lines in 3d on white slices
% 			     - 3 = Schlieren plot on slices
% 			     - 4 = scatter plot of q vs. r
%
%        mq                  - which component of q to plot
%        UserVariable        - Set to 1 to specify a user defined variable.
%        UserVariableFile    - name of m-file mapping data to q
%        MappedGrid          - set to 1 if mapc2p.m exists for nonuniform grid
%        MaxFrames           - max number of frames
%        MaxLevels           - max number of AMR levels
%        PlotData            - Data on refinement level k is plotted only if
%                              PlotData(k) == 1
%        PlotGrid            - PLot grid lines on level k is PlotGrid(k) /= 0
%        PlotGridEdges       - Plot 2d patch borders if PlotGridEdges(k) /= 0
%        PlotCubeEdges       - Plot 3d patch cubes if PlotCubeEdges(k) /= 0
%        ContourValues       - Set to desired contour values, or [] for no ...
% 	                     lines.
%        xSliceCoords        - vector of x slice constants
%        ySliceCoords        - vector of y slice constants
%        zSliceCoords        - vector of z slice constants
%        x0,y0,z0            - center for scatter plots.
%        ScatterStyle        - symbols for scatter plots.
%        LineStyle           - same as ScatterStyle.
%        IsosurfValues       - constants for isosurfaces
%        IsosurfColors       - colors for isosurfaces.
%        UserMap1d           - set to 1 if 'map1d' file exists.
%
%      All parameters can be modified by typing 'k' at the PLOTCLAW3 prompt.
%
%      See also PLOTCLAW3, setPlotGrid, setPlotGridEdges.


setviews;  % set viewpoints so that view(xSlice), for example, can be used.


PlotType = 1;                % type of plot to produce:
			     % 1 = pcolor on slices (with optional contours,
			     % and isosurfaces)
			     % 2 = contour lines in 3d on transparent slices
			     % 3 = Schlieren plot on slices
			     % 4 = scatter plot of q vs. r

mq = 1;                      % which component of q to plot
UserVariable = 0;            % set to 1 to specify a user-defined variable
UserVariableFile = ' ';      % name of m-file mapping data to q
MappedGrid = 0;              % set to 1 if mapc2p.m exists for nonuniform grid
MaxFrames = 1000;            % max number of frames to loop over
MaxLevels = 6;               % max number of AMR levels

PlotData =  [1 1 1 0 0 0];       % Data on refinement level k is plotted only
			         % if k'th component is nonzero
PlotGrid =  [0 0 0 0 0 0];       % Plot grid lines on each level?
PlotGridEdges =  [0 0 0 0 0 0];  % Plot edges of patches of each grid at
                                 % this level on slices?
PlotCubeEdges = [0 0 0 0 0 0];   % Plot edges of cube of refinement patch at
                                 % this level?

% ---------------------------------------------------------------------
% The next three parameters are vectors of x,y,z coordinates of 2d slices
% to be displayed for PlotType = 1,2,3.
% Empty ==> no slices in that direction.

xSliceCoords = 0.5;
ySliceCoords = 0.5;
zSliceCoords = 0.5;

% ---------------------------------------------------------------------
% ContourValues is a vector of values used to draw contour lines.
% The valid settings for this parameter are identical to those used by the
% Matlab contour plotting routine.  See also CONTOUR.  In particular:
%   If ContourValues is the empty matrix, no contour lines will be drawn.
%   If ContourValues is a vector, these values will be used for contours.
%   If ContourValues is an integer, this number of levels will be drawn,
%        with values chosen based on the data.  (May not work well with
%        AMR data since different levels may be chosen on different grids.

ContourValues = [];

% ---------------------------------------------------------------------
% Isosurfaces.  Plot surfaces at q = IsosurfValue(i)
% If empty, no isosurfaces will be drawn.

IsosurfValues    =  [];  

IsosurfColors    = 'b';  % Colors for each surface.
                         % Set to 'q' to get colors from the
                         % current colormaps.  Use STRVCAT to get
                         % multiple colors, i.e. strvcat('b','r','g','y');

% ---------------------------------------------------------------------
% plot q(r) vs. r = sqrt((x-x0)^2 + (y-y0)^2 + (z-z0)^2);
% for scatter plot (PlotType==4):
% The default is to plot q(i,j,k) vs. 
%             r(i,j,k) = (x(i,j,k)-x0)^2 + (y(i,j,k)-y0)^2 + (z(i,j,k)-z0)^2
x0 = 0.5;
y0 = 0.5;
z0 = 0.5;

UserMap1d = 0;  % set to 1 and provide map1d.m file to specify a different
                % mapping of (x,y,z,q) to (r,q) for scatter plots.

ScatterStyle = setplotstyle('ro','bx','k.','rs','bv','k^');
%               determines symbol and color for plotting scatter data on
%               each refinement level.
