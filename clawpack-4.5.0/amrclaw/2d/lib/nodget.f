c
c ------------------------------------------------------------
c
      integer function nodget(dummy)
c
      implicit double precision (a-h,o-z)

      include  "call.i"
c
c ::::::::::::::::: NODGET ::::::::::::::::::::::::::::::::::::;
c nodget =  get first free node of the linked list kept in node
c            array. adjust pointers accordingly.
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;
c
      if (ndfree .ne. null) go to 10
          write(outunit,100) maxgr
          write(*,100)       maxgr
100       format(' out of nodal space - allowed ',i5,' grids')
          stop
c
 10     nodget         = ndfree
        ndfree         = node(nextfree,ndfree)
c
c  initialize nodal block
c
        do 20 i        = 1, nsize
           node(i,nodget) = 0
 20     continue
c
        do 30 i         = 1, rsize
           rnode(i,nodget) = 0.0d0
 30     continue
c
      return
      end
