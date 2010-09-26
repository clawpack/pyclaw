%  SETPLOT1, SETPLOT2, SETPLOT3  sets user defined plotting parameters for
%           1d, 2d or 3d plots.
%
% The three scripts are SETPLOT1, SETPLOT2, SETPLOT3.m
%
%      User defined Matlab script for setting various Clawpack plotting
%      parameters.  This script is called by PLOTCLAW<N>.  A default
%      version of this script can be found in claw/matlab/setplot<N>.m and
%      copied to users working directory and modifed to set things up
%      differently.
%
% ---------------------------------------------------------------------------
%      Parameters that can be set with SETPLOT1
%
%        OutputFlag        - set to 'ascii' (default) to read ascii files,
%                            'hdf' to read hdf files.
%        PlotStyle           - used in plot command for line color and type.
%        mq                  - which component of q to plot
%        UserVariable        - Set to 1 to specify a user defined variable.
%        UserVariableFile    - name of m-file mapping data to q
%        MappedGrid          - set to 1 if mapc2p.m exists for nonuniform
%                              grid
%        MaxFrames           - max number of frames
%
%      All parameters can be modified by typing 'k' at the PLOTCLAW3 prompt.
%
%      See also PLOTCLAW1, SETPLOT1, SETSYMBOLS.
% ---------------------------------------------------------------------------
%      Parameters that can be set with SETPLOT2
%
%        OutputFlag        - set to 'ascii' (default) to read ascii files,
%                            'hdf' to read hdf files.
%        PlotType          - type of plot to produce:
% 			     - 1 = pcolor on slices (with optional contours)
% 			     - 2 = contour lines in 2d on white slices
% 			     - 3 = Schlieren plot on slices
% 			     - 4 = scatter plot of q vs. r
%
%        mq                  - which component of q to plot
%        UserView            - camera viewing angle.  Can be set to anything
%                              that is valid for the Matlab VIEW command.
%        UserVariable        - Set to 1 to specify a user defined variable.
%        UserVariableFile    - name of m-file mapping data to q
%        MappedGrid          - set to 1 if mapc2p.m exists for nonuniform
%                              grid
%        Manifold            - set to 1 if mapc2m.m exists for manifold plot.
%        MaxFrames           - max number of frames
%        MaxLevels           - max number of AMR levels
%        PlotData            - Data on refinement level k is plotted only if
%                              PlotData(k) == 1
%        PlotGrid            - PLot grid lines on level k is PlotGrid(k) /= 0
%        PlotGridEdges       - Plot 2d patch borders if PlotGridEdges(k) /= 0
%        ContourValues       - Set to desired contour values, or [] for no ...
% 	                       lines.
%        ScatterStyle        - symbols to be used for scatter plots.
%        LineStyle           - same as ScatterStyle
%        x0,y0               - center for scatter plots.
%        UserMap1d           - set to 1 if file 'map1d' exists.
%
%      All parameters can be modified by typing 'k' at the PLOTCLAW3 prompt.
%
%      See also PLOTCLAW2, SETPLOT2, MANIFOLD, MAPPEDGRID, SETSYMBOLS.
% ----------------------------------------------------------------------------
%      Parameters that can be set with SETPLOT3
%
%        OutputFlag        - set to 'ascii' (default) to read ascii files,
%                            'hdf' to read hdf files.
%        PlotType          - type of plot to produce:
% 			     - 1 = pcolor on slices (with optional contours)
% 			     - 2 = contour lines in 3d on white slices
% 			     - 3 = Schlieren plot on slices
% 			     - 4 = scatter plot of q vs. r
%
%        mq                  - which component of q to plot
%        UserView            - camera viewing angle.  Can be set to anything
%                              that is valid for the Matlab VIEW command.
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
%        ContourValues       - Set to desired contour values, or [] for no
% 	                       contour lines.
%        xSliceCoords        - vector of x slice constants
%        ySliceCoords        - vector of y slice constants
%        zSliceCoords        - vector of z slice constants
%        ScatterStyle        - symbols to be used with scatter plots.
%        LineStyle           - same as ScatterStyle.
%        x0,y0,z0            - center for scatter plots.
%        UserMap1d           - set to 1 if file 'map1d' exists.
x%        IsosurfValues       - constants for isosurfaces
%        IsosurfColors       - colors for isosurfaces.
%
%      All parameters can be modified by typing 'k' at the PLOTCLAW3 prompt.
%
%      See also PLOTCLAW3, SETPLOT3, MANIFOLD, MAPPEDGRID, SETPLOTSTYLE.


error(['setplot : This is routine cannot be called on its own. Please ',...
       'use SETPLOT1, SETPLOT2, or SETPLOT3 to set parameters');
