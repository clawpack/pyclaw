c
c     =================================================
      function stream(x,y,z)
c     =================================================
c
c     # stream function defining the velocity.
c     # (x,y,z) is a point on the sphere.
c     # rotation about an axis through 0 and (xaxis,yaxis,zaxis)

      implicit double precision (a-h,o-z)
      common /comaxis/ xaxis,yaxis,zaxis

      
      stream = x*xaxis + y*yaxis + z*zaxis
c
      return
      end
