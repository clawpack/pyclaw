c
c ---------------------------------------------------------------
c
        subroutine trimbd(used,nrow,ncol,set,il,ir,jb,jt)
c
c :::::::::::::::::::::::: TRIMBD :::::::::::::::::::::::::::;
c  if used array is completely set (=1.) then return set=true, 
c  otherwise return false, alogn with the dimensions of the smallest 
c  rectangle containing all unset points in il,ir,jb,jt.
c ::::::::::::::::::::::::::::::::::::;::::::::::::::::::::::;
c
        implicit double precision (a-h,o-z)
        dimension  used(nrow,ncol)
        logical   set

        utot = 0.
        do 100 j = 1,ncol
        do 100 i = 1,nrow
100        utot = utot + used(i,j)
        if (utot .ge. dble(nrow*ncol)) then
                set = .true.
                return
        endif
 
        set = .false.
 
        uleft = 1.
        do 200 i = 1,nrow
           do 220 j = 1,ncol
              uleft = dmin1(uleft,used(i,j))
220        continue
           il = i
           if (uleft .eq. 0.) go to 230
200     continue

230     uright = 1.
        do 300 i = 1,nrow
           do 320 j = 1,ncol
              uright = dmin1(uright,used(nrow - i + 1,j))
320        continue
           ir = nrow - i + 1
           if (uright .eq. 0.) go to 330
300     continue

330     ubot = 1.
        do 400 j = 1,ncol
           do 420 i = 1,nrow
              ubot = dmin1(ubot,used(i,j))
420        continue
           jb = j
           if (ubot .eq. 0.) go to 430
400        continue
 
430     utop = 1.
        do 500 j = 1,ncol
           do 520 i = 1,nrow
              utop = dmin1(utop,used(i,ncol - j + 1))
520        continue
           jt = ncol - j + 1
           if (utop .eq. 0.) go to 530
500     continue
 
530     return
        end
