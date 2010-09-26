c
c --------------------------------------------------------------------
c
       subroutine signs(badpts,npts,iscr,jscr,idim,jdim,ist,iend,
     &                  ilo,ihi,jlo,jhi)
c
       implicit double precision (a-h,o-z)
       dimension badpts(2,npts)
       dimension iscr(idim), jscr(jdim)
c
c :::::::::::::::::::::::::::: SIGNS ::::::::::::::::::::::::::::::
c  compute signatures = number of flagged cells in each row/column.
c  also return first and last nonzero row/column, so don't have
c  to waste time over entire region.
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
       ilo= 1
       ihi= idim
       jlo= 1
       jhi= jdim
       do 10 i = 1, idim
 10       iscr(i) = 0
       do 15 j = 1, jdim
 15       jscr(j) = 0
c
c count all flagged points in a given row/column in one pass through
c the points, i.e. a bin count
c
       do 20 ipt = ist, iend
          iloc  =  badpts(1,ipt)+1.1
          jloc  =  badpts(2,ipt)+1.1
          iscr(iloc) = iscr(iloc)+1
          jscr(jloc) = jscr(jloc)+1
 20    continue
c
       do 30 ipt = 1, idim
         if (iscr(ipt) .ne. 0) then
            ilo = ipt
            go to 40
         endif
 30    continue
 40    do 50 ipt = 1, idim
          if (iscr(idim+1-ipt) .ne. 0) then
             ihi = idim+1-ipt
             go to 60
          endif
 50    continue

 60    do 70 ipt = 1, jdim
          if (jscr(ipt) .ne. 0) then
             jlo = ipt
             go to 80
          endif
 70    continue
 80    do 90 ipt = 1, jdim
          if (jscr(jdim+1-ipt) .ne. 0) then
             jhi = jdim+1-ipt
             go to 99
          endif
 90    continue

 99    return
       end
