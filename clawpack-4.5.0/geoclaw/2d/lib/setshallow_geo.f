c=========================================================================
      subroutine settsunami
c=========================================================================

      implicit double precision (a-h,o-z)
      character*25 fname
      logical foundFile

      include "call.i"
      include "geo.i"

      write(parmunit,*) ' '
      write(parmunit,*) '--------------------------------------------'
      write(parmunit,*) 'SETSHALLOW:'
      write(parmunit,*) '-----------'

c       # read user parameters from setshallow.data

      fname  = 'setshallow.data'
      inquire(file=fname,exist=foundFile)
      if (.not. foundFile) then
        write(*,*) 'You must provide a file ', fname
        stop
      endif

      iunit = 7
      call opendatafile(iunit, fname)

      read(7,*) sealevel
      read(7,*) drytolerance
      read(7,*) depthdeep
      read(7,*) maxleveldeep
      read(7,*) ifriction
      read(7,*) coeffmanning
      read(7,*) frictiondepth
      close(7)

      write(parmunit,*) '   drytolerance:',drytolerance
      write(parmunit,*) '   maxleveldeep:', maxleveldeep
      write(parmunit,*) '   depthdeep:', depthdeep
      write(parmunit,*) '   ifriction:', ifriction
      write(parmunit,*) '   Manning coefficient:',coeffmanning
      write(parmunit,*) '   frictiondepth:',frictiondepth

      return
      end
