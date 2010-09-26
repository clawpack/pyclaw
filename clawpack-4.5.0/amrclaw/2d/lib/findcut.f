c
c -----------------------------------------------------------
c
      subroutine findcut(icl,iscr,jscr,idim,jdim,index,iside,
     1                   ilo,ihi,jlo,jhi)
c
c ::::::::::::::::::::: FINDCUT ::::::::::::::::::::::::::::;
c   find best place to split the 2D array of flagged points
c   either split at a hole, or use signatures to find
c   zero crossing of laplacian.
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;
c
      implicit double precision (a-h,o-z)

      include "call.i"

      dimension iscr(idim), jscr(jdim)

c Modified 6/02:
c Include call.i to get def's of horizontal/vertical.
c      integer   horizontal, vertical
c      parameter(horizontal = 1)
c      parameter(vertical = 2)

      parameter(ithres = 2)
      parameter(minoff = 2)
c
c  look for holes first in horizontal then vertical direction
c
       do 10 i = ilo, ihi
          if (iscr(i) .eq. 0) then
             index = i
             iside = horizontal
             return
          endif
 10    continue

       do 20 j = jlo, jhi
          if (jscr(j) .eq. 0) then
              index = j
              iside = vertical
              return
          endif
 20    continue

c
c  no holes - find 2nd derivative of signatures for best cut.
c  overwrite signature arrays. don't make cuts less than minoff
c  from boundary
c
      ipre = iscr(ilo)
      do 50 i = ilo+1, ihi-1
         icur = iscr(i)
         iscr(i) = iscr(i+1)-2*icur+ipre
         ipre = icur
 50   continue

      locmaxi = 0
      indexi  = 0
      imid    = (ilo + ihi) / 2
      do 60 i = ilo+minoff, ihi-minoff+1
           itemp1 = iscr(i-1)
           itemp2 = iscr(i)
           locdif = iabs(itemp1-itemp2)
           if (itemp1*itemp2.lt.0) then
                if (locdif .gt. locmaxi) then
                 locmaxi = locdif
                 indexi = i
                else if (locdif .eq. locmaxi) then
                    if (iabs(i-imid).lt.iabs(indexi-imid)) indexi = i
                endif
           endif
 60   continue


       jpre   = jscr(jlo)
       do 130 j = jlo+1, jhi-1
           jcur = jscr(j)
           jscr(j) = jscr(j+1)-2*jcur+jpre
           jpre = jcur
 130   continue

       locmaxj = 0
       indexj  = 0
       jmid    = (jlo + jhi) / 2
       do 160 j = jlo+minoff, jhi-minoff+1
           jtemp1 = jscr(j-1)
           jtemp2 = jscr(j)
           locdif = iabs(jtemp1-jtemp2)
           if (jtemp1*jtemp2.lt.0) then
               if (locdif .gt. locmaxj) then
                  locmaxj = locdif
                  indexj = j
                else if (locdif .eq. locmaxj) then
                       if (iabs(j-jmid).lt.iabs(indexj-jmid)) indexj = j
               endif
           endif
 160   continue
c
c      ::::: choose max dif for splitting
c
 220   if (locmaxi .gt. locmaxj) then
            index  = indexi
            iside  = horizontal
            locmax = locmaxi
       else if (locmaxi .lt. locmaxj) then
            index  = indexj
            iside  = vertical
            locmax = locmaxj
        else if (iabs(indexi-imid).lt.iabs(indexj-jmid)) then
                 index  = indexi
                 iside  = horizontal
                 locmax = locmaxi
             else
               index  = indexj
               iside  = vertical
               locmax = locmaxj
       endif

c      ::::: if inflection pt. not over the threshold, signal
c      ::::: with index 0 (out of range)
       if (locmax .lt. ithres) index = 0

       return
       end
