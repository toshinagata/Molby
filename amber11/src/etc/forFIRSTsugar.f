C=====================================================================
C
      subroutine addsugarpsatoms(MAXFHB, MAXATOM, MAXRES, natom, nres, 
     .                     igraph, ipres, lbres, ib, jb, nbond, c, 
     .                     ftype, fhybrid, 
     .                     nhb, fhbdon, fhbh, fhbacc, fhbene)
C
C Two pseudoatoms are inserted into the sugar rings of nucleic acids
C   resulting into a 7 membered ring with 1 dof.
C   The pseudoatoms are inserted between C1'-O4' and C4'-O4'. 
C   The former corresponding bonds are removed, respectively. 
C Due to the increased ring size hc between C1' and C5' of the same 
C   sugar ring are explicitely forbidden in forFIRSTteth.f
C Note that only standard nucleosides are considered.
C
C Simone Fulle - 26.02.2008
C
C
      implicit none
C
      integer MAXFHB, MAXATOM, MAXRES, 
     +        MAXDEPTH, MAXTETH, MAXDOF, pseudonum
      integer natom, nres, ipres, ib, jb, nbond, nb
      integer nrs, nat, m
      integer fhybrid, nhb, fhbdon, fhbh, fhbacc
      integer d, i, r, j, k, l, t, n
      integer neigh, neighnum, neighnumfst
      integer resid, atom
      double precision c, fhbene, vecfac
      double precision ene
      double precision x1, y1, z1, x2, y2, z2, xd, yd, zd
      double precision atinfo, array
      character*4 igraph, lbres, aname
      character*1 ftype
      logical go
C
C     --- fixed number of pseudo atoms applied in sugar rings *2---
      parameter(pseudonum = 1)
C
      dimension igraph(*), ipres(*), lbres(*), ib(*), jb(*)
      dimension ftype(*), fhybrid(*), fhbdon(*), fhbh(*), 
     +          fhbacc(*), fhbene(*)
      dimension c(*)
      dimension neigh(MAXATOM), neighnum(MAXATOM), 
     +          neighnumfst(MAXATOM)
      dimension atinfo(natom,3)
C          
C     --- Default energy value for tether ---
      data ene/-9.999999/
C
C       --- Factor to calc coords of tether atoms ---
      vecfac = 1.0 / dble(pseudonum + 1)
C
      nrs = nres
      nat = natom
      nb = nbond
C      
C       --- Loop over all residues ---
      do resid=1,nrs
C       -- in the case of a nucleic acid
        if(lbres(resid)(1:1).eq.'R' .or. 
     +     lbres(resid)(1:1).eq.'D') then
C         --- Loop over all atoms i in residue r ---
           do i=ipres(resid),(ipres(resid+1)-1)
              if(igraph(i).eq.'O4''') then
	         atom = i
                 x1 = c(3*i-2)
                 y1 = c(3*i-1)
                 z1 = c(3*i  )
              endif
           enddo
C
C          -- neighbor of C4' has index "atom-2" and of C1' has index "atom+1" --
           do n=1,2
             if(n.eq.1) then
               j = atom-2
             else
               j = atom+1
             endif
           
             x2 = c(3*j-2)
             y2 = c(3*j-1)
             z2 = c(3*j  )
C                                                    
             xd = x2 - x1
             yd = y2 - y1
             zd = z2 - z1
C         
             nres = nres + 1
             lbres(nres) = 'BMH'
             ipres(nres) = nat + 1
C
C            --- Update atom and bond information ---                  
             xd = xd * vecfac
             yd = yd * vecfac
             zd = zd * vecfac
C
             nat = nat + 1
             if(nat.gt.MAXATOM) then 
               write(6,*) 'Too many atoms: ', nat
               stop
             endif
             k = pseudonum
C
             c(3*nat-2) = x1 + k * xd
             c(3*nat-1) = y1 + k * yd
             c(3*nat  ) = z1 + k * zd
C                      
             igraph(nat) = 'X'
             ftype(nat) = 'N'
             fhybrid(nat) = 0
             nbond = nbond + 1
             if(nbond.gt.MAXATOM) then 
               write(6,*) 'Too many bonds: ', nbond
               stop
             endif
C            -- insert bond between pseudoatom = nat and O4' = atom --         
             ib(nbond) = (atom-1)*3 
             jb(nbond) = (nat-1)*3
C         
C            -- instead of bond between C4' and O4 --
C            -- insert either bond between pseudoatom = nat and C4' = atom-2 --
             if(n.eq.1) then 
               do m=1,nb
                 if(ib(m).eq.((atom-3)*3) .and.
     +              jb(m).eq.((atom-1)*3)) then
C                     ib(m) = (atom-3)*3   
                      jb(m) = (nat-1)*3
                 endif
               enddo
C            -- ... or bond between pseudoatom = nat and C1' = atom+1 --
             else
               do m=1,nb
                 if(ib(m).eq.((atom-1)*3) .and.
     +             jb(m).eq.((atom)*3)) then
                     ib(m) = (nat-1)*3   
                 endif
               enddo           
             endif
C         
             natom = nat
             if((nres+1).gt.MAXRES) then
               write(6,*) 'Too many residues: ', nres+1
               stop
             endif
             ipres(nres+1) = natom+1
           enddo
        endif
      enddo
      return
      end
C
C=====================================================================

