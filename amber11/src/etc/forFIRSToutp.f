C=====================================================================
C
      subroutine writebond(NF, nbond, ib, jb)
C
C Writes bond list for FIRST input file
C
C Holger Gohlke
C   06.12.2001
C
      implicit none
C
      integer NF, nbond, ib, jb
      integer i, i1, i2
C
      dimension ib(*), jb(*)
C
      write(NF, 100)
      write(NF, 110)
C
      do i=1,nbond
        i1 = ib(i)/3+1
        i2 = jb(i)/3+1
        if(i1.lt.i2) then
          write(NF, 120) i1, i2
        else
          write(NF, 120) i2, i1
        endif
      enddo
C
 100  format('REMARK:L:----------------------------------------------',
     +       '----------------------------------') 
 110  format('REMARK:cf','    so','    sf','     (atom-label pairs)') 
 120  format('REMARK:CF',2i6)
C
      return
      end
C
C=====================================================================
C
      subroutine writetf(NF, ntf, itf, jtf)
C
C Writes tf list for FIRST input file
C
C Holger Gohlke
C   06.12.2001
C
      implicit none
C
      integer NF, ntf, itf, jtf
      integer i, i1, i2
C
      dimension itf(*), jtf(*)
C
      write(NF, 100)
      write(NF, 110)
C
      do i=1,ntf
        i1 = itf(i)/3+1
        i2 = jtf(i)/3+1
        if(i1.lt.i2) then
          write(NF, 120) i1, i2
        else
          write(NF, 120) i2, i1
        endif
      enddo
C
 100  format('REMARK:L:----------------------------------------------',
     +       '----------------------------------') 
 110  format('REMARK:tf','    so','    sf','     (atom-label pairs)') 
 120  format('REMARK:TF',2i6)
C
      return
      end
C
C=====================================================================
C
      subroutine writehb(NF, nhb, ftype, fhybrid, fhbdon, fhbh, 
     +                   fhbacc, fhbene, igraph, hbene, fhbunit)
C
C Writes HB for FIRST input file
C
C Holger Gohlke
C   06.12.2001
C
      implicit none
C
      integer NF, nhb, fhybrid, fhbdon, fhbh, fhbacc
      integer i, iold, icurr, iacc, idon, j, cnt
      integer fhbunit
      double precision fhbene, ene, eneold, tmp, hbene
      character*4 igraph
      character*1 ftype
C
      dimension ftype(*), fhybrid(*), fhbdon(*), fhbh(*),
     +          fhbacc(*), fhbene(*), igraph(*), fhbunit(*)
C
      write(NF, 100)
      write(NF, 110)
C
      ene = fhbene(1)
      icurr = 1
      do i=2,nhb
        if(fhbene(i).gt.ene) then
          ene = fhbene(i)
          icurr = i
        endif
      enddo
C
      cnt = 0
      do i=1,nhb
        iacc = fhbacc(icurr)
        idon = fhbdon(icurr)
        if(fhbene(icurr) .le. hbene) then
          cnt = cnt + 1
          if(ftype(iacc).eq.'E' .and. ftype(idon).eq.'C' .or.
     +       ftype(iacc).eq.'C' .and. ftype(idon).eq.'E') then
            write(NF, 130) cnt, fhbene(icurr), fhbdon(icurr),
     +        fhbh(icurr), fhbacc(icurr)
          else if(igraph(iacc)(1:1).eq.'X' .or. 
     +            igraph(idon)(1:1).eq.'X') then
            write(NF, 140) cnt, fhbene(icurr), fhbdon(icurr),
     +        fhbh(icurr), fhbacc(icurr)
          else
            write(NF, 120) cnt, fhbene(icurr), fhbdon(icurr),
     +        fhbh(icurr), fhbacc(icurr), fhybrid(idon), fhybrid(iacc), 
     +        fhbunit(icurr)
          endif
        endif
C
        eneold = ene
        iold = icurr
        ene = -1.0d10
        do j=1,nhb
          tmp = fhbene(j)
          if(tmp.eq.eneold .and. j.gt.iold) then
            icurr = j
            ene = tmp
            goto 10
          else if(tmp.lt.eneold .and. tmp.gt.ene) then
            icurr = j
            ene = tmp
          endif          
        enddo
   10   continue
      enddo   
C
 100  format('REMARK:L:----------------------------------------------',
     +       '----------------------------------') 
 110  format('REMARK:hb','   ID','   E(Kcal/Mol)','  Donor',
     +       '  Hydrogen',' Acceptor','  Description') 
 120  format('REMARK:HB',i5,1x,f12.5,3(1x,i7),4x,
     +       'HB',1x,'Dsp',i1,1x,'Asp',i1,2x,'hbtype:',1x,i1)
 130  format('REMARK:HB',i5,1x,f12.5,3(1x,i7),4x,
     +       'SB',1x,'no energy')
 140  format('REMARK:HB',i5,1x,f12.5,3(1x,i7),4x,
     +       'PH',1x,'hydr',1x,'phob')
C
      return
      end
C
C=====================================================================
