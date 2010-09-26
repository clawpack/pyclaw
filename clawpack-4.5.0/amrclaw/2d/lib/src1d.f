c
c
c =========================================================
      subroutine src1d(meqn,mbc,mx1d,q1d,maux,aux1d,t,dt)
c =========================================================
      implicit double precision (a-h,o-z)
      dimension   q1d(mx1d, meqn)
      dimension aux1d(mx1d, maux)
c
c     # dummy source routine... does nothing
c
c     # This routine should be a simplified version of src2 
c     # which applies source terms for a 1-d slice of data along the
c     # edge of a grid.  This is called only from qad where the conservative
c     # fix-up is applied and is used to apply source terms over partial
c     # time steps to the coarse grid cell values used in solving Riemann 
c     # problems at the interface between coarse and fine grids.
c
c     # If the source terms depend only on q, it should be easy to 
c     # adapt src2 to create this routine, just loop over 1:mx1d.
c     # If the source terms are more complicated, it may not be easy.
c
c     # The code may work fine without applying source terms in this
c     # context, so using this dummy routine might be successful even when
c     # source terms are present. 
c
      return
      end
