% Routines for plotting data output from Clawpack package
%
% Basic Clawpack plotting routines and options.
%
%      plotclaw1              - plot 1d results
%      plotclaw2              - plot 2d results
%      plotclaw3              - plot 3d results
%      setplot1               - script for setting parameters for 1d plots
%      setplot2               - script for setting parameters for 2d plots
%      setplot3               - script for setting parameters for 3d plots
%      setplot                - help on using SETPLOT<N>
%      readamrdata            - reads output files produced by Clawpack.
%      plotframe1ez           - functional form of plotframe1.
%                               useful for adding a 1d plot to existing plot
%      setplotstyle           - sets symbols and colors for line/scatter plots
%      queryframe             - queries user for next frame information
%      NoQuery                - setting NoQuery=1 avoids query between plots
%      getlegendinfo          - returns legend information on line plots.
%      printgif               - prints a gif file from current figure
%      printjpg               - prints a jpg file from current figure
%      makeframegif           - print current figure to file frame0000N.gif
%      makeframejpg           - print current figure to file frame0000N.jpg
%                                 where N = current value of Frame counter.
%
% Data analysis routines (for use with UserVariable==1)
%      pressure               - returns pressure given input data (gamma law)
%      xvelocity              - returns x velocity
%      yvelocity              - returns y velocity
%      zvelocity              - returns z velocity
%      mach                   - return mach number.
%
% General graph properties for 2d and 3d plots.
%
%      showcontourlines       - shows contour lines
%      hidecontourlines       - hides contour lines
%      showgridlines          - shows computational grid
%      hidegridlines          - hides computational grid
%      showpatchborders       - shows patch borders
%      hidepatchborders       - hides patch borders
%      showmesh               - shows a coarsened mesh on specified levels
%      hidemesh               - hides a coarsened mesh on specified levels
%      mappedgrid             - parameter for plotting mapped grids.
%      setcolors              - gives user control over how colors are set
%      setopengl              - Sets OpenGL renderer
%      showslices             - shows slices/manifold
%      hideslices             - hides slices/manifold
%      setslicecolor          - sets color of slice/manifold
%      setslicealpha          - set transparency value of slice/manifold
%      yrbcolormap            - yellow/red/blue colormap
%      redwhite               - red/white colormap
%      rybcolormap            - red/yellow/blue colormap
%      getblocknumber         - get block number for plotting results of
%                               multi-block calculations.
%
% 2d specific graphics routines
%
%      showlevels             - shows specified levels (2d only)
%      hidelevels             - hides specified levels (2d only)
%      manifold               - parameter for plotting manifold.
%      projectcontours        - projects contour lines to user-specified plane.
%
% 3d specific graphics routines
%
%      sliceloop              - loop over slices on 3d plots.
%      surfloop               - loop over isosurfaces
%      setviews               - sets pre-defined viewing angles
%      showcubes              - shows 3d amr patch cube borders
%      hidecubes              - hides 3d amr patch cube borders
%      showsurfs              - shows isosurfaces created with ISOSURFVALUES
%      hidesurfs              - hides isosurfaces created with ISOSURFVALUES
%      showsurflevels         - show isosurfaces at specified AMR levels
%      hidesurflevels         - hide isosurfaces at specified AMR levels
%      reducesurf             - reduces number of faces on isosurface.
%      showsurfmesh           - shows isosurface mesh
%      hidesurfmesh           - hides isosurface mesh
%      setsurfalpha           - sets isosurface transparency value
%      setsurfcolor           - sets isosurface color
%
% Type 'help' on any one of the individual topics above for more help.
%
