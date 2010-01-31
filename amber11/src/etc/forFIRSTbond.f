C=====================================================================
C
      subroutine corbondl(MAXATOM, NATOM, IGRAPH, NRES, IPRES, LBRES, 
     +                    nbond, ib, jb, c)
C
C Corrects bond list:
C   - removes bonds between H's of water
C   - adds bonds between metal ion and surroundings
C
C Holger Gohlke
C   06.12.2001
C
      implicit none
C
      integer MAXATOM, natom, nres, ipres, nbond, ib, jb
      integer i, i1, i2, j, k, k1, k2, l
      double precision c, dist, dist2 
      double precision x1, y1, z1, x2, y2, z2, xd, yd, zd, tmp
      character*4 igraph, lbres, lbr
      logical found
C
      dimension igraph(*), ipres(*), lbres(*), ib(*), jb(*), c(*)
C
C     --- Distance and square of distance of metal - surroundings ---
      data dist, dist2/2.3, 5.29/
C
C     --- Remove bonds between H's in water ---
C
      do i=1,nres
        lbr = lbres(i)
        if(lbr.eq.'WAT' .or. 
     +     lbr.eq.'HOH' .or. lbr.eq.'H2O' .or. 
     +     lbr.eq.'OH2' .or.
     +     lbr.eq.'DOD' .or. lbr.eq.'D2O' .or. 
     +     lbr.eq.'OD2') then
          i1 = ipres(i)
          i2 = ipres(i+1) - 1
          l = 0
          do j=i1,i2
            if(igraph(j)(1:1).eq.'H' .or. igraph(j)(1:1).eq.'D') then
              do k=1, nbond
                k1 = ib(k)/3+1
                k2 = jb(k)/3+1
                if((k1.eq.j .and. (igraph(k2)(1:1).eq.'H' .or. 
     +                             igraph(k2)(1:1).eq.'D')) .or.
     +             (k2.eq.j .and. (igraph(k1)(1:1).eq.'H' .or. 
     +                             igraph(k1)(1:1).eq.'D'))) then
                  l = k
                  goto 10
                endif
              enddo
            endif
          enddo
   10     continue
          if(l.gt.0) then
ccc            write(6,*) 'Removing bond between ', ib(l)/3+1, jb(l)/3+1
            ib(l) = ib(nbond)
            jb(l) = jb(nbond)
            ib(nbond) = 0
            jb(nbond) = 0
            nbond = nbond - 1
          endif
        endif
      enddo
C
C     --- Add bonds between metal ions and surroundings ---
C
      do i=1,natom
        found = .false.
        do j=1,nbond
          i1 = ib(j)/3+1
          i2 = jb(j)/3+1
          if(i.eq.i1 .or. i.eq.i2) then
            found = .true.
            goto 20
          endif
        enddo
   20   continue
        if(.not. found) then
C         --- Atom w/o bonds found => assuming metal ion ---
ccc          write(6,*) 'Found metal ion ', i
          x1 = c(3*i-2)
          y1 = c(3*i-1)
          z1 = c(3*i  )
          do k=1,natom
            if(i.ne.k) then
              x2 = c(3*k-2)
              y2 = c(3*k-1)
              z2 = c(3*k  )
              xd = x2 - x1
              yd = y2 - y1
              zd = z2 - z1
              if(abs(xd).lt.dist .or. 
     +           abs(yd).lt.dist .or. 
     +           abs(zd).lt.dist) then
                tmp = xd*xd + yd*yd + zd*zd
                if(tmp.lt.dist2) then
ccc                  write(6,*) '  Found surrounding ', k
                  nbond = nbond + 1
                  if(nbond.gt.MAXATOM) then
                    write(6,*) 'Too many bonds: ', nbond
                    stop
                  endif
                  ib(nbond) = (i-1)*3
                  jb(nbond) = (k-1)*3
                endif
              endif
            endif
          enddo
        endif
      enddo
C
      return
      end
C        
C=====================================================================
C
      subroutine findtf(MAXATOM, nbond, ib, jb, fhybrid, igraph,
     +                  ntf, itf, jtf)
C
C Finds bonds which have to be tethered
C
C Holger Gohlke
C   06.12.2001
C
      implicit none
C
      integer MAXATOM, MAXRSIZE
      integer nbond, ib, jb, fhybrid, fh1, fh2
      integer ntf, itf, jtf
      integer i, i1, i2, j, j1, j2, n1, n2
      character*4 igraph
      logical tether, isinring, inring1, inring2, webstyle
C
      parameter(MAXRSIZE = 6)
C
      dimension ib(*), jb(*), fhybrid(*), itf(*), jtf(*), igraph(*)
C
C     --- TRUE: CZ-NH{1,2} in ARG are tethered,
C               but not C(O)-NH2 bonds in ASN, GLN ---
      data webstyle/.true./
C
      ntf = 0
      do i=1,nbond
        i1 = ib(i)/3 + 1
        i2 = jb(i)/3 + 1
        if(fhybrid(i1).eq.2 .and. fhybrid(i2).eq.2) then
          tether = .true.
C         --- Don't tether sp2-sp2 bonds with one terminal atom ---
          n1 = 0
          n2 = 0
          do j=1,nbond
            if(i.ne.j) then
              j1 = ib(j)/3 + 1
              j2 = jb(j)/3 + 1
              fh1 = fhybrid(j1)
              fh2 = fhybrid(j2)
              if(webstyle) then
                if((fh2.gt.0 .or. igraph(j1)(1:2).eq.'NH') .and. 
     +             i1.eq.j1 .or.
     +             (fh1.gt.0 .or. igraph(j2)(1:2).eq.'NH') .and.
     +             i1.eq.j2) n1 = n1 + 1
                if((fh2.gt.0 .or. igraph(j1)(1:2).eq.'NH') .and.
     +             i2.eq.j1 .or.
     +             (fh1.gt.0 .or. igraph(j2)(1:2).eq.'NH') .and.
     +             i2.eq.j2) n2 = n2 + 1
              else
                if(fh2.gt.0 .and. i1.eq.j1 .or.
     +             fh1.gt.0 .and. i1.eq.j2) n1 = n1 + 1
                if(fh2.gt.0 .and. i2.eq.j1 .or.
     +             fh1.gt.0 .and. i2.eq.j2) n2 = n2 + 1
              endif
              if(n1.gt.0 .and. n2.gt.0) goto 10
            endif
          enddo
   10     continue
          if(n1.eq.0 .or. n2.eq.0) tether = .false.
C         --- Don't tether bonds in rings ---
          if(tether) then
            inring1 = isinring(MAXRSIZE, 0, nbond, ib, jb, i1, 0, i1)
            inring2 = isinring(MAXRSIZE, 0, nbond, ib, jb, i2, 0, i2)
            if(inring1 .and. inring2) tether = .false.
          endif
C         --- Tether bond ---
          if(tether) then
            ntf = ntf + 1
            if(ntf.gt.MAXATOM) then
              write(6,*) 'Too many tether bonds: ', ntf
              stop
            endif
            itf(ntf) = (i1-1)*3
            jtf(ntf) = (i2-1)*3
          endif
        endif
      enddo
C
      return
      end
C
C=====================================================================
C
      recursive
     +  function isinring(MAXRSIZE, rsize, 
     +                    nbond, ib, jb, at, from, curr)
     +  result(inring)
C
C Determines if atom at is in ring of maximum size MAXRSIZE
C
C Holger Gohlke
C   06.12.2001
C
      implicit none
C
      integer MAXRSIZE, rsize, nbond, ib, jb
      integer at, from, curr
      integer i, i1, i2
      logical inring
C
      dimension ib(*), jb(*)
C
      inring = .false.
      if(rsize.gt.0 .and. curr.eq.at) then
ccc        write(6,*) 'Atom in ring ', at
        inring = .true.
      else if((rsize+1).le.MAXRSIZE) then
        do i=1,nbond
          i1 = ib(i)/3+1
          i2 = jb(i)/3+1
          if(i1.eq.curr .and. i2.ne.from) then
            inring = inring .or.
     +               isinring(MAXRSIZE, rsize+1, 
     +                        nbond, ib, jb, at, curr, i2)
          else if(i2.eq.curr .and. i1.ne.from) then
            inring = inring .or. 
     +               isinring(MAXRSIZE, rsize+1, 
     +                        nbond, ib, jb, at, curr, i1)
          endif
          if(inring) goto 10
        enddo
   10   continue
      endif
C
      return
      end
C
C=====================================================================
C
