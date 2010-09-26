
      subroutine opendatafile(iunit, fname)
c     #
c     # Open the file fname and determine how many leading lines are
c     # comments.  Then rewind and skip over the comment lines so that the
c     # file is ready for reading data from.
c     #
c     # All comment lines must start with # in the first column.

      integer iunit, line, commentlines
      character*12 fname
      character*1 firstchar

c     write(6,*) 'Check correct filename: XXX',fname,'XXX'
      open(unit=iunit,file=fname,status='old',form='formatted')
      firstchar = '#'
      commentlines = -1
      do while (firstchar == '#')
         read(iunit,*) firstchar
         commentlines = commentlines + 1
         enddo

      write(6,601) commentlines, fname
  601 format('Reading data file, first',i2,' lines are comments: ',a12)

      rewind(iunit)
      do line=1,commentlines
         read(iunit,*) firstchar
         enddo

      return
      end

