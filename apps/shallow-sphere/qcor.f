c
c     =====================================================
      subroutine qcor(ixy,i,m,aux,q,maxm,meqn,mbc,qc)
c     =====================================================
c
c     # maintain conservation on sphere
c
      implicit double precision (a-h,o-z)

      dimension aux(16, 1-mbc:maxm+mbc)
      dimension q(meqn, 1-mbc:maxm+mbc)
      dimension qc(4)

      common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom
      common /sw/  g


      if(ixy.eq.1) then
        in = 2
        dy = dycom
      else
        in = 8
        dy = dxcom 
      endif
c     # left interface 
      etxl = aux(in+3,i)
      etyl = aux(in+4,i)
      etzl = aux(in+5,i)
      gammal = dsqrt(etxl**2.d0 + etyl**2.d0 + etzl**2.d0) / dy
      enxl = aux(in,i) * gammal
      enyl = aux(in+1,i) * gammal
      enzl = aux(in+2,i) * gammal
    

c     # right interface
      etxr = aux(in+3,i+1)
      etyr = aux(in+4,i+1)
      etzr = aux(in+5,i+1)        
      gammar = dsqrt(etxr**2.d0 + etyr**2.d0 + etzr**2.d0) / dy
      enxr = aux(in,i+1) * gammar
      enyr = aux(in+1,i+1) * gammar
      enzr = aux(in+2,i+1)  * gammar
c
      qc(1) = (enxr-enxl)*q(2,i)
     &       +(enyr-enyl)*q(3,i)
     &       +(enzr-enzl)*q(4,i)

      qc(2) = (enxr-enxl)*(q(2,i)**2.d0/q(1,i)+0.5d0*g*q(1,i)**2.d0)
     &       +(enyr-enyl)*(q(2,i)*q(3,i)/q(1,i))
     &       +(enzr-enzl)*(q(2,i)*q(4,i)/q(1,i))
    
      qc(3) = (enxr-enxl)*(q(2,i)*q(3,i)/q(1,i))
     &       +(enyr-enyl)*(q(3,i)**2.d0/q(1,i)+0.5d0*g*q(1,i)**2.d0)
     &       +(enzr-enzl)*(q(3,i)*q(4,i)/q(1,i))
   
      qc(4) = (enxr-enxl)*(q(2,i)*q(4,i)/q(1,i))
     &       +(enyr-enyl)*(q(3,i)*q(4,i)/q(1,i))
     &       +(enzr-enzl)*(q(4,i)**2.d0/q(1,i)+0.5d0*g*q(1,i)**2.d0)

c     # project qc to tangent plane:
      erx = aux(14,i)
      ery = aux(15,i)
      erz = aux(16,i)

      qcn = erx*qc(2) + ery*qc(3) + erz*qc(4)
      
      qc(2) = qc(2) - qcn*erx
      qc(3) = qc(3) - qcn*ery
      qc(4) = qc(4) - qcn*erz

      return
      end
