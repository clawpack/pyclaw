function setcolors(p,x,y,z,q)

% SETCOLORS determines mapping between q data and color map
%
%      This subroutine, called by the Clawpack graphics routines, sets the
%      way in which q data is mapped to the graphics colormap.  This routine
%      is in the same location as the rest of the Clawpack graphics files;
%      if you wish to modify this file, copy claw/matlab/setcolors.m to your
%      working directory, and make any desired changes.
%
%      The syntax for this routine is
%
%              setcolors(p,x,y,z,q);
%
%      where p is the handle to the patch whose colors are being set,
%      (x,y,z) are the Cartesian locations of cell centers with corresponding
%      data values in q.  If a mapped grid or manifold is being plotted,
%      (x,y,z) are given in the Cartesian locations, NOT the physical
%      locations of cell centers.
%
%      By default, the q values are mapped linearly into current colormap,
%      with min and max values clamped to limits set by caxis.
%
%      This default behavior is accomplished with the following Matlab
%      commands :
%
%            set(p,'CData',q);                % Data to use for coloring.
%            set(p,'CDataMapping','scaled');  % Scale into current color map.
%            set(p,'FaceColor','flat');       % Single color per cell
%
%      Other color mapping schemes can be provided, for example to mask out
%      embedded boundary regions, or to flag certain values that lie out
%      side of a given data range.
%
%      For example, to highlight all values that lie outside of a given
%      range [a,b] (to see where overshoots and undershoots occur, for example).
%      Assume for example that the current colormap has length 'n', and that
%      in location 1 you have assigned the color black ([0 0 0]) and in
%      location n you have assigned the color white ([1 1 1]).  You want to
%      color all values q < a the color black, and all values q > b the
%      color white.  Everything in the range [a,b] should be mapped into the
%      colormap 'default'.  The following code will do this:
%
%      Example :
%
%           colormap('default');
%           cmap = [[0 0 0]; colormap; [1 1 1]];  % set black, white values
%           colormap(cmap);
%           n = length(cmap);
%
%           % Assign a color into map (2:n-1) based on value of q:
%           a = 0;  % value data range.
%           b = 1;
%
%           % Map a -> index value 2
%           % Map b -> index value n-1
%           idx = ((n-1)-2)/(b-a)*(q - a) + 2;
%
%           % Map all values that fall outside of the range
%           % [a,b] the color black or white
%           qcolors(q < a) = 1;   % assign values < a the color black.
%           qcolors(q > b) = n;   % assign values > b the color white.
%           qcolors(q >= a & q <= b) = round(idx);
%
%           % Set CData property for patch
%           set(p,'CData',qcolors);
%
%           % Interpret values in CData as indices which can be used
%           % directly into colormap.
%           set(p,'CDataMapping','direct');
%
%      See also PATCH, COLORMAP.

set(p,'CData',q);                % Data to use for coloring.
set(p,'CDataMapping','scaled');  % Scale into current color map.
set(p,'FaceColor','flat');       % Single color per cell
