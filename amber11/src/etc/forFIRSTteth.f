C=====================================================================
C
      subroutine addtether(MAXFHB, MAXATOM, MAXRES, natom, nres, 
     .                     igraph, ipres, lbres, ib, jb, nbond, c, 
     .                     ftype, fhybrid, 
     .                     nhb, fhbdon, fhbh, fhbacc, fhbene)
C
C Adds tether atoms to simulate hydrophobic interactions.
C
C Holger Gohlke - 06.12.2001
C
C
C Base stacking interactions between bases are limited to 1. 
C   This avoids over-rigidification of nucleic acid structures.
C The threshold for hydrophobic interactions in proteins is set 
C   to 0.25, and in nucleic acids to 0.15. 
C Note that the parametrization for nucleic acids has been
C   tested so far only on RNA structures!
C
C Simone Fulle - 26.02.2008
C
C
      implicit none
C
      integer MAXFHB, MAXATOM, MAXRES, 
     +        MAXNEIGH, MAXDEPTH, MAXTETH, MAXDOF
      integer natom, nres, ipres, ib, jb, nbond
      integer fhybrid, nhb, fhbdon, fhbh, fhbacc
      integer d, i, i2, ith, j, k, l, t, n, nat, locali, localj
      integer neigh, neighnum, neighnumfst, tethnum, residjth
      integer residi, residj, tethercut, count, residith, nrs
      double precision c, fhbene, vecfac, cutoff
      double precision thresh, threshprot, threshrna
      double precision vdw, vdw1, vdw2, sum, sum2, dist2, ene
      double precision x1, y1, z1, x2, y2, z2, xd, yd, zd
      double precision atinfo, array
      character*4 igraph, lbres, aname
      character*1 ftype
      character*3 modnuc
      logical go
C
C     --- Either MAXTETH or MAXDOF must be >= 0 ---
C     --- MAXTETH >= 0: fixed number of tethering atoms applied ---
C     --- MAXDOF  >= 0: number of degrees of freedom REMOVED per tether from the system ---
      parameter(MAXNEIGH = 14, MAXDEPTH = 3, MAXTETH = 3, MAXDOF = -1)
C
      dimension igraph(*), ipres(*), lbres(*), ib(*), jb(*)
      dimension ftype(*), fhybrid(*), fhbdon(*), fhbh(*), 
     +          fhbacc(*), fhbene(*)
      dimension c(*), vdw(2), modnuc(7) 
      dimension neigh(MAXNEIGH, MAXATOM), neighnum(MAXATOM), 
     +          neighnumfst(MAXATOM)
      dimension atinfo(natom,3)
      dimension array((nres*2),5)
C          
C     --- Bondi radii for C, S ---
      data vdw/1.7, 1.8/
C     --- Threshold: dist < vdw1 + vdw2 + thresh <=> hydrophobic contact ---
      data threshprot/0.25/
      data threshrna/0.15/
C     --- Default energy value for tether ---
      data ene/-9.999999/
C     --- number of tethers considered for base stacking = {1,2} ---
      data tethercut/1/
C     --- explicitely considered modified RNA nucleosides ---
      data modnuc/'PSU', 'OMA', 'OMG', 'OMC', 'OMU', '1MA', 'UR3'/
c
      do i=1,MAXATOM                                                  
         neighnum(i) = 0
      enddo
C
C     --- Check consistency of MAXTETH, MAXDOF
C
      if((MAXTETH.ge.0 .and. MAXDOF.ge.0) .or.
     +   (MAXTETH.lt.0 .and. MAXDOF.lt.0)) then
        write(6,*) 'Either MAXTETH or MAXDOF must be >= 0', 
     +             MAXTETH, MAXDOF
        stop
      endif
C
C     --- Precalc some values for fixed MAXTETH
C
      if(MAXTETH.ge.0) then
        tethnum = MAXTETH
C       --- Factor to calc coords of tether atoms ---
        vecfac = 1.0 / dble(tethnum + 1)
      endif
C
C     --- Build list of bonded neighb up to MAXDEPTH ---
C
      do i=1,natom
        if(igraph(i)(1:1).eq.'C' .or. igraph(i)(1:1).eq.'S') then
          neighnum(i) = 0
          neighnumfst(i) = 0
          call buildneigh(MAXNEIGH, MAXATOM, MAXDEPTH, 
     +                    i, 0, i, nbond, ib, jb, igraph,
     +                    0, neigh, neighnum, neighnumfst)
        endif
      enddo
C
C     --- Find tether atoms. Tethers can not be between atoms 
C         of bond distance <= MAXDEPTH ---
C
      nat = natom
      nrs = nres
      residi = 1
      residj = 2
      residjth = 0
      residith = nres+1
      cutoff = 10000.0
      do i=1, (nrs*2)
        array(i,5) = cutoff
      enddo
C
C     -----------------------------
C     --- Loop over all atoms i ---
      do i=1,natom-1
C
        if(i.eq.ipres(residi+1)) then
          residi = residi+1
          residjth=0
        endif
        if(i.eq.(ipres(residi+1)-1)) then
          residj = residi+1
        else 
          residj=residi
        endif
C
C       --- set threshold for hc ---
        if(lbres(residi)(1:1).eq.'R' .or.
     +     lbres(residi)(1:1).eq.'D' .or.
     +     lbres(residi).eq.'PSU' .or.
     +     lbres(residi)(1:2).eq.'OM' .or.
     +     lbres(residi).eq.'1MA' .or.
     +     lbres(residi).eq.'UR3') then
           thresh = threshrna
        else
           thresh = threshprot 
        endif
C
        aname = igraph(i)
        if(aname(1:1).eq.'C' .or. aname(1:1).eq.'S') then
          if(aname(1:1).eq.'C') then
            vdw1 = vdw(1)
          else
	    vdw1 = vdw(2)
          endif
          x1=c(3*i-2)
          y1=c(3*i-1)
          z1=c(3*i  )
C
C         ---------------------------------
C         --- Loop over all atoms j > i ---
          do j=i+1,natom
            if(j.eq.ipres(residj+1)) then
              residj = residj+1
            endif
C
C 	   ---  set threshold for hc between nucleic acid and protein units ---
            if(lbres(residj)(1:1).ne.'R' .and.
     +         lbres(residj)(1:1).ne.'D' .and.
     +         lbres(residi).ne.'PSU' .and.
     +         lbres(residi)(1:2).ne.'OM' .and.
     +         lbres(residi).ne.'1MA' .and.
     +         lbres(residi).ne.'UR3') then
               thresh = threshprot
            endif
C-----------------------------------------------------------------
C           --- for all carbon atoms in nucleic acid bases ---
            if ((igraph(i).eq.'C2' .or. igraph(i).eq.'C4' .or.
     +          igraph(i).eq.'C5' .or. igraph(i).eq.'C6' .or.
     +          igraph(i).eq.'C8') .and.
     +          (igraph(j).eq.'C2' .or. igraph(j).eq.'C4' .or.
     +          igraph(j).eq.'C5' .or. igraph(j).eq.'C6' .or.
     +          igraph(j).eq.'C8')) then
C
C               --- which do not belong to the same nucleotide---
                if (residi.ne.residj) then
                  aname = igraph(j)
                  go = .true.
                  n = neighnum(j)
C                 --- Check that no bond distance <= MAXDEPTH ---
                  do k=1,n
                    if(neigh(k,j).eq.i) then
                      go = .false.
                      goto 10
                    endif
                  enddo
   10             continue
                  if(go) then
                    vdw2 = vdw(1)
C                   --- Check distance ---
                    sum = vdw1 + vdw2 + threshrna
                    sum2= sum * sum
                    x2=c(3*j-2)
                    y2=c(3*j-1)
                    z2=c(3*j  )
                    xd = x2 - x1
                    yd = y2 - y1
                    zd = z2 - z1

	            if(abs(xd).lt.sum .or. 
     +                abs(yd).lt.sum .or.
     +                abs(zd).lt.sum) then
                      dist2 = xd*xd + yd*yd + zd*zd                                 
C	           
                      if(dist2.lt.sum2) then
C                     --- candidate for hc ---
                        do ith=1, residjth
                          if(residj.eq.array(ith*2-1,2)) then
                            residith = ith
                            goto 20
                          endif
C                         --- new residue
                          residith = residjth+1 
                        enddo
   20                   continue
C		--- save per base stacking 2 hc with smallest distance ---
                        if(residith.le.residjth) then
C                          --- old residue j ---
                           if(dist2.lt.array(residith*2-1,1)) then
                              array(residith*2,1)=array(residith*2-1,1)
                              array(residith*2,2)=array(residith*2-1,2)
                              array(residith*2,3)=array(residith*2-1,3)
                              array(residith*2,4)=array(residith*2-1,4)
                              array(residith*2-1,1)=dist2
                              array(residith*2-1,2)=residj
                              array(residith*2-1,3)=i
                              array(residith*2-1,4)=j
                           else if(dist2.lt.array(residith*2,1)) then
                                array(residith*2,1)=dist2
                                array(residith*2,2)=residj
                                array(residith*2,3)=i
                                array(residith*2,4)=j
                           endif
                        else
C                          --- new residue j ---
                           residjth=residjth+1
                           array(residjth*2-1,1)=dist2
                           array(residjth*2-1,2)=residj
                           array(residjth*2-1,3)=i
                           array(residjth*2-1,4)=j
                        endif  
                      endif
                    endif  
                  endif
                endif  
CC              --- end different residues ---
C                
C-----------------------------------------------------------------
C           --- within the sugar rings ---
            else if((igraph(i).eq.'C1''' .and. igraph(j).eq.'C5''' .or.
     +              igraph(i).eq.'C5''' .and. igraph(j).eq.'C1''') .and.
     +              residi.eq.residj) then 
C           forbit tether between C1' and C5' of the same residue    
C           required if 2 pseudoatoms are inserted into one sugar ring
C
C-----------------------------------------------------------------
C           --- insert tether for all other cases ---
            else
		aname = igraph(j)
                go = .true.
                n = neighnum(j)
C               --- Check that no bond distance <= MAXDEPTH ---
                do k=1,n
                  if(neigh(k,j).eq.i) then
                    go = .false.
                    goto 30
                  endif
                enddo
   30           continue
                if(go .and.
     +	           (aname(1:1).eq.'C' .or. aname(1:1).eq.'S')) then
                  if(aname(1:1).eq.'C') then
                    vdw2 = vdw(1)
                  else
	            vdw2 = vdw(2)
                  endif
C                 --- Check distance ---
                  sum = vdw1 + vdw2 + thresh
                  sum2= sum * sum
                  x2=c(3*j-2)
                  y2=c(3*j-1)
                  z2=c(3*j  )
                  xd = x2 - x1
                  yd = y2 - y1
                  zd = z2 - z1
	          if(abs(xd).lt.sum .or. 
     +               abs(yd).lt.sum .or.
     +               abs(zd).lt.sum) then
                    dist2 = xd*xd + yd*yd + zd*zd
	            if(dist2.lt.sum2) then
C
C                     --- Add tether atoms ---
C  
C                     --- Update residue information ---
                      nres = nres + 1
                      if(nres.gt.MAXRES) then
                        write(6,*) 'Too many residues: ', nres
                        stop
                      endif
                      lbres(nres) = 'BMH'
                      ipres(nres) = nat + 1
C                     --- Calc some values if MAXDOF >= 0 --- 
                      if(MAXDOF.ge.0) then
                        tethnum = -MAXDOF + 1 + neighnumfst(i) + 
     +                            neighnumfst(j)
                        vecfac = 1.0 / dble(tethnum + 1)
                      endif
C                     --- Update atom and bond information ---                  
                      xd = xd * vecfac
                      yd = yd * vecfac
                      zd = zd * vecfac
                      do k=1,tethnum
                        nat = nat + 1
                        if(nat.gt.MAXATOM) then
                          write(6,*) 'Too many atoms: ', nat
                          stop
                        endif
                        c(3*nat-2) = x1 + k * xd
                        c(3*nat-1) = y1 + k * yd
                        c(3*nat  ) = z1 + k * zd
                        igraph(nat) = 'X'
                        ftype(nat) = 'N'
                        fhybrid(nat) = 0
                        nbond = nbond + 1
                        if(nbond.gt.MAXATOM) then
                          write(6,*) 'Too many bonds: ', nbond
                          stop
                        endif
                        if(k.eq.1) then
                          ib(nbond) = (i-1)*3
                        else
	                  ib(nbond) = (nat-2)*3
                        endif
                        jb(nbond) = (nat-1)*3
                      enddo
C                     --- Update H-bond information ---
                      nhb = nhb+1
                      if(nhb.gt.MAXFHB) then
                        write(6,*) 'Too many H-bonds: ', nhb
                        stop
                      endif
                      fhbdon(nhb) = nat - 1
                      fhbh(nhb) = nat
                      fhbacc(nhb) = j
                      fhbene(nhb) = ene
                    endif
                  endif
                endif 
            endif  
          enddo  
C-----------------------------------------------------------------
C   
C         --- add hc between bases --
          if(((lbres(residi)(1:2).eq.'RC' .or. 
     +         lbres(residi)(1:2).eq.'RU' .or.
     +         lbres(residi)(1:2).eq.'DC' .or. 
     +         lbres(residi)(1:2).eq.'DT' .or. 
     +         lbres(residi).eq.'PSU' .or.
     +         lbres(residi).eq.'OMC' .or.
     +         lbres(residi).eq.'OMU' .or.
     +         lbres(residi).eq.'UR3') .and.
     +         igraph(i).eq.'C2') .or.
     +       ((lbres(residi)(1:2).eq.'RA' .or. 
     +         lbres(residi)(1:2).eq.'RG' .or.
     +         lbres(residi)(1:2).eq.'DA' .or. 
     +         lbres(residi)(1:2).eq.'DG' .or. 
     +         lbres(residi).eq.'OMA' .or.
     +         lbres(residi).eq.'OMG' .or.
     +         lbres(residi).eq.'1MA') .and.
     +         igraph(i).eq.'C4')) then
C     
              do t=1,tethercut
                if(t.eq.1) then
                  d=1
                else
                  d=0  
                endif
C
                do count=1,residjth
                   locali=array(count*2-d,3)
                   localj=array(count*2-d,4)
                   vdw1=vdw(1)
                   vdw2=vdw(1)
                   sum = vdw1 + vdw2 + threshrna
                   sum2= sum * sum
                   x1 = c(3*locali-2)
                   y1 = c(3*locali-1)
                   z1 = c(3*locali  )
                   x2 = c(3*localj-2)
                   y2 = c(3*localj-1)
                   z2 = c(3*localj  )
                   xd = x2 - x1
                   yd = y2 - y1
                   zd = z2 - z1
                   dist2 = array(count*2-d,1)
C
C                  --- Add tether atoms ---
                   if(dist2.lt.cutoff .and. dist2.lt.sum2) then
C                    --- Update residue information ---
                     nres = nres + 1
                     if(nres.gt.MAXRES) then
                       write(6,*) 'Too many residues: ', nres
                       stop
                     endif
                     lbres(nres) = 'BMH'
                     ipres(nres) = nat + 1
C                    --- Calc some values if MAXDOF >= 0 --- 
                     if(MAXDOF.ge.0) then
                       tethnum = -MAXDOF + 1 + neighnumfst(locali) + 
     +                 neighnumfst(localj)
                       vecfac = 1.0 / dble(tethnum + 1)
                     endif
C                    --- Update atom and bond information ---                  
                     xd = xd * vecfac
                     yd = yd * vecfac
                     zd = zd * vecfac
                     do k=1,tethnum
                       nat = nat + 1
                       if(nat.gt.MAXATOM) then
                         write(6,*) 'Too many atoms: ', nat
                         stop
                       endif
                       c(3*nat-2) = x1 + k * xd
                       c(3*nat-1) = y1 + k * yd
                       c(3*nat  ) = z1 + k * zd
                       igraph(nat) = 'X'
                       ftype(nat) = 'N'
                       fhybrid(nat) = 0
                       nbond = nbond + 1
                       if(nbond.gt.MAXATOM) then
                         write(6,*) 'Too many bonds: ', nbond
                         stop
                       endif
                       if(k.eq.1) then
                         ib(nbond) = (locali-1)*3
                       else
	                 ib(nbond) = (nat-2)*3
                       endif
                       jb(nbond) = (nat-1)*3
                     enddo
C                    --- Update H-bond information ---
                     nhb = nhb+1
                     if(nhb.gt.MAXFHB) then
                       write(6,*) 'Too many H-bonds: ', nhb
                       stop
                     endif
                     fhbdon(nhb) = nat - 1
                     fhbh(nhb) = nat
                     fhbacc(nhb) = localj
                     fhbene(nhb) = ene
                   endif  
                enddo
              enddo
              do i2=1, (nrs*2)
                 array(i2,1) = cutoff
                 array(i2,2) = 0
              enddo
          endif         
        endif  
      enddo
C     
      natom = nat
      if((nres+1).gt.MAXRES) then
        write(6,*) 'Too many residues: ', nres+1
        stop
      endif
      ipres(nres+1) = natom+1
C
      return
      end
C
C=====================================================================
C
      recursive 
     +  subroutine buildneigh(MAXNEIGH, MAXATOM, MAXDEPTH, 
     +                        at, from, curr, nbond, ib, jb, igraph,
     +                        ndepth, neigh, neighnum, neighnumfst)
C
C Recursively find neighbors of iat
C
C Holger Gohlke
C   06.12.2001
C
      implicit none
C
      integer MAXNEIGH, MAXATOM, MAXDEPTH
      integer at, from, curr, to, nbond, ib, jb
      integer ndepth, neigh, neighnum, neighnumfst, nn
      integer currdepth
      integer i, i1, i2, j
      character*4 igraph
      logical go
C
      dimension neigh(MAXNEIGH, MAXATOM)
      dimension neighnum(MAXATOM), neighnumfst(MAXATOM)
      dimension ib(*), jb(*), igraph(*)
C
      currdepth = ndepth + 1
      if(currdepth.le.MAXDEPTH) then
        do i=1,nbond
          i1 = ib(i)/3+1
          i2 = jb(i)/3+1
          go = .false.
          if(i1.eq.curr .and. i2.ne.from) then
            to = i2
            go = .true.
          else if(i2.eq.curr .and. i1.ne.from) then
            to = i1
            go = .true.
          endif

          if(go) then
C           --- Store number of neighbors in first shell ---
            if(currdepth.eq.1) then
              neighnumfst(at) = neighnumfst(at) + 1
            endif
C           --- Store neighbor information ---
            if(igraph(to)(1:1).eq.'C' .or. 
     +         igraph(to)(1:1).eq.'S') then
              nn = neighnum(at)
              nn = nn + 1
              if(nn.gt.MAXNEIGH) then
                write(6,*) 'Too many neighbors: ', nn
                stop
              endif
              neighnum(at) = nn
              neigh(nn, at) = to
ccc              write(6,*) ndepth+1, at, ' (', nn,') ', to
            endif
            call buildneigh(MAXNEIGH, MAXATOM, MAXDEPTH,
     +                      at, curr, to, nbond, ib, jb, igraph,
     +                      currdepth, neigh, neighnum, neighnumfst)
          endif
        enddo
      endif
C
      return
      end
C
C=====================================================================
