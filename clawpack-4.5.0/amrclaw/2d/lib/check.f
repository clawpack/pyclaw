c
c ---------------------------------------------------------
c
      subroutine check(nsteps,time,nvar,naux)
c
c :::::::::::::::::::::: CHECK ::::::::::::::::::::::::::::::::;
c   check point routine - can only call at end of coarse grid cycle
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;

      implicit double precision (a-h,o-z)
      integer tchkunit
      parameter (tchkunit = 12)
      character  chkname*12
      character  tchkname*12
      include  "call.i"
c
c     ###  make the file name showing the time step
c
      chkname = 'fort.chkxxxx'
      tchkname = 'fort.tckxxxx'
      nstp = nsteps
      do 20 ipos = 12, 9, -1
         idigit = mod(nstp,10)
         chkname(ipos:ipos) = char(ichar('0') + idigit)
         tchkname(ipos:ipos) = char(ichar('0') + idigit)
         nstp = nstp / 10
 20   continue

      open(unit=tchkunit,file=tchkname,status='unknown',
     .     form='formatted')
      open(unit=chkunit,file=chkname,status='unknown',
     .     form='unformatted')
c
c     ###  dump the data
c
      write(tchkunit,*) 'Checkpoint file at time t = ',time
      write(tchkunit,*) 'alloc size memsize = ',memsize
      write(tchkunit,*) 'Number of steps taken = ',nsteps
      write(chkunit) lenmax,lendim,memsize
      write(chkunit) (alloc(i),i=1,lendim)
      write(chkunit) hxposs,hyposs,possk,icheck
      write(chkunit) lfree,lenf
      write(chkunit) rnode,node,lstart,newstl,listsp,tol,
     1          ibuff,mstart,ndfree,lfine,iorder,mxnest,
     2          intratx,intraty,kratio,iregsz,jregsz,
     2          iregst,jregst,iregend,jregend, 
     3          kcheck,nsteps,
     3          time,matlabu
      write(chkunit) evol, rvol, rvoll, lentot
c
      close(chkunit)
      close(tchkunit)
c
      return
      end
