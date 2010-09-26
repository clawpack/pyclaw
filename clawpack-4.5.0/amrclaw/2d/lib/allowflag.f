
c     =========================================
      logical function allowflag(x,y,t,level)
c     =========================================

c     # Indicate whether the grid point at (x,y,t) at this refinement level
c     # is allowed to be flagged for further refinement.
c
c     # This is useful if you wish to zoom in on some structure in a 
c     # known location but don't want the same level of refinement elsewhere.  
c     # Points are flagged only if one of the errors is greater than the 
c     # corresponding tolerance.
c
c     # For example, to allow refinement of Level 1 grids everywhere but
c     # of finer grids only for  y >= 0.4:
c     # allowed(x,y,t,level) = (level.le.1 .or. y.ge.0.4d0) 
c
c     # This routine is called from routine flag2refine.
c     # If Richardson error estimates are used (if tol>0) then this routine
c     # is also called from errf1.

      implicit double precision (a-h,o-z)

c     # default version allows refinement anywhere:
      allowflag = .true.

      return
      end
