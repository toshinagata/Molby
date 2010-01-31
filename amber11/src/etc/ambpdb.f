c
c   Simple filter program to convert AMBER coordinate files into PDB files
c   (and to do other similar conversions).
c
      implicit none
      character(len=80) arg,prmtop,title
      character(len=8)  hbenec
      character(len=6)  arg1,arg2
      character(len=4), dimension(:), allocatable ::  igraph,lbres
      character(len=1), dimension(:), allocatable ::  ftype
      integer access, ntf, ioffset, iarg, indx, iargc, natom, nres
      integer natomlen, nbondlen, nreslen
      integer nbond, ier, kcform, idbl, j, nhb
      character(len=4) ititl(20)
      logical alttit,center,aatm,bres,ter,bin,hasradii,remediate
      
C     --- required for -first option ---
      integer MAXATOM, MAXRES, MAXFHB
c
      double precision, dimension(:), allocatable :: C,fhbene
      real, dimension(:), allocatable :: chg
      real, dimension(:), allocatable :: radius
      double precision hbene
c
      integer, allocatable, dimension(:) :: residue_number
      character(len=4), allocatable, dimension(:) :: residue_chainid,
     &   atom_altloc, atom_element, residue_icode
      logical :: ext_pdb_data
c 
      integer, dimension(:), allocatable :: ipres,lastat,ib,jb,fhybrid,
     .     fhbdon,fhbh,fhbacc,fhbunit,itf,jtf
c
c      ------ check argument options
c
      arg1 = '-pdb'
      prmtop = 'prmtop'
      alttit = .false.
      hbene = 1.0E10
      center = .false.
      aatm = .false.
      bres = .false.
      ter = .true.
      bin = .false.
      remediate = .true.
      ext_pdb_data = .false.
      ioffset = 0
      iarg = 0
      indx = iargc()
      if (indx.eq.0) go to 20
   10 continue
           iarg = iarg + 1
           call getarg(iarg,arg)
           if (arg .eq. '-p') then
                iarg = iarg + 1
                call getarg(iarg,prmtop)
           else if (arg .eq. '-tit') then
                iarg = iarg + 1
                call getarg(iarg,title)
                alttit = .true.
           else if (arg .eq. '-ene') then
                iarg = iarg + 1
                call getarg(iarg,hbenec)
                read(hbenec,'(f8.0)' ) hbene                
           else if (arg .eq. '-sas') then
                arg1 = '-sas'
           else if (arg .eq. '-pqr') then
                arg1 = '-pqr'
           else if (arg .eq. '-bnd') then
                arg1 = '-bnd'
           else if (arg .eq. '-atm') then
                arg1 = '-atm'
           else if (arg .eq. '-first') then
                 arg1 = '-first'
           else if (arg .eq. '-ctr') then
                center = .true.
           else if (arg .eq. '-aatm') then
                aatm = .true.
           else if (arg .eq. '-bres') then
                bres = .true.
           else if (arg .eq. '-noter') then
                ter = .false.
           else if (arg .eq. '-noremediate') then
                remediate = .false.
           else if (arg .eq. '-ext') then
                ext_pdb_data = .true.
           else if (arg .eq. '-bin') then
                bin = .true.
           else if (arg .eq. '-offset') then
                iarg = iarg + 1
                call getarg(iarg,arg2)
                read(arg2,'(i5)' ) ioffset
           else if (arg .eq. '-h' .or. arg .eq. '-help') then
                call usage()
           else
                write(6,'(/,5x,a,a)') 'unknown flag: ',arg
                stop
           endif
      if (iarg .lt. indx) go to 10
c
   20 continue
C
C       ----- OPEN THE PARM FILE AND LOAD THE NECESSARY STUFF -----
C
      call amopen(10,prmtop,'O','F','R')
      CALL GETNAM0(NATOM,NRES,NBOND,10)
      REWIND(UNIT=10)
c
c      ------ Allocate memory: ------
c    
      natomlen = natom
      nbondlen = nbond
      nreslen = nres
      if (arg1.eq.'-first') then
        MAXATOM = 2*natom
        MAXRES = 2*natom
        MAXFHB = natom
        natomlen = MAXATOM
        nbondlen = MAXATOM
        nreslen = MAXRES
      endif

      allocate( c(3*natomlen), igraph(natomlen), ipres(nreslen+1), 
     .          lbres(nreslen), lastat(nreslen), 
     .          ib(nbondlen), jb(nbondlen),
     .          chg(natomlen), ftype(natomlen), fhybrid(natomlen),
     .          fhbdon(natomlen), fhbh(natomlen), fhbacc(natomlen),
     .          fhbene(natomlen), fhbunit(natomlen), 
     .          itf(natomlen), jtf(natomlen),
     .          radius(natomlen),
     .          residue_number(nreslen),
     .          residue_chainid(nreslen), atom_altloc(natomlen),
     .          atom_element(natomlen), residue_icode(nreslen),
     .          stat = ier)
      if (ier /= 0 ) then
         write(0,*) 'memory allocation failure'
         call mexit(6,1)
      end if
c
      CALL GETNAM(NATOM,NRES,IGRAPH,IPRES,LBRES,ITITL,10,1,
     1           C,C,C,ib,jb,nbond,chg,radius,lastat,ter,hasradii,
     2           residue_number, residue_chainid, atom_altloc,
     3           atom_element, residue_icode, ext_pdb_data)
      CLOSE(UNIT=10)
C
C       ----- READ THE COORDINATE FILE -----
C
      CALL GETCOR(bin,NATOM,C,C,KCFORM,IDBL,5)
c
      if(arg1.eq.'-first') then      
c
c       ----- Prepare for FIRST output -----
c
        call getftype(NATOM, NRES, IGRAPH, IPRES, LBRES, 
     .                ib, jb, nbond, ftype)
        call gethybrid(NATOM, NRES, IGRAPH, IPRES, ib, jb, 
     .                 nbond, c, fhybrid)
        call findtf(MAXATOM, nbond, ib, jb, fhybrid, igraph,
     .              ntf, itf, jtf)
        call corbondl(MAXATOM, NATOM, IGRAPH, NRES, IPRES, LBRES,
     .                nbond, ib, jb, c)
        call findhbond(MAXFHB, NATOM, NRES, IGRAPH, IPRES, LBRES, 
     .                 ib, jb, nbond, c, ftype, fhybrid, 
     .                 nhb, fhbdon, fhbh, fhbacc, fhbene, fhbunit)
C       call addsugarpsatoms(MAXFHB, MAXATOM, MAXRES, NATOM, NRES,
C    .                 IGRAPH, IPRES, LBRES, ib, jb, nbond, c,
C    .                 ftype, fhybrid,
C    .                 nhb, fhbdon, fhbh, fhbacc, fhbene)
        call addtether(MAXFHB, MAXATOM, MAXRES, NATOM, NRES,
     .                 IGRAPH, IPRES, LBRES, ib, jb, nbond, c, 
     .                 ftype, fhybrid, 
     .                 nhb, fhbdon, fhbh, fhbacc, fhbene)
      endif
C
C       ----- OUTPUT THE PDB, atm, pqr, first, bnd or, cpk  FILE -----
C
      if (arg1.eq.'-pdb' .or. arg1.eq.'-atm' 
     .      .or. arg1.eq.'-pqr' .or. arg1.eq.'-sas'
     .      .or. arg1.eq.'-first') then
        CALL GENPDB(NATOM,NRES,C,IGRAPH,IPRES,LBRES,ITITL,6,arg1,
     .     alttit,title,center,chg,aatm,bres,ioffset,lastat,
     .     ftype,radius,hasradii,remediate,
     .     residue_number, residue_chainid, atom_altloc,
     .     atom_element, residue_icode, ext_pdb_data)
        if(arg1.eq.'-first') then
          call writebond(6, nbond, ib, jb)
          call writetf(6, ntf, itf, jtf)
          call writehb(6, nhb, ftype, fhybrid, fhbdon, fhbh, 
     .                 fhbacc, fhbene, igraph, hbene, fhbunit)
        endif
C
      else if (arg1.eq.'-bnd') then
        do j=1,nbond
          write(6,'(2i5)') ib(j)/3+1,jb(j)/3+1
        end do
C
      else if (arg1.eq.'-cpk') then
        call cpkgen(c,natom,igraph)
C
      else if (arg1.eq.'-plu') then
        call pluto(c,natom,igraph,alttit,title)
C
      else
        write(6,*) 'Needs -pdb, -atm, -bnd, -cpk, or -first flag'
        stop
      end if
C
      deallocate( c, igraph, ipres, lbres, lastat, ib, jb, chg,
     .          ftype, fhybrid, fhbdon,fhbh, fhbacc, fhbene, itf, jtf,
     .          stat = ier)
      if (ier /= 0 ) then
         write(0,*) 'memory deallocation failure'
         call mexit(6,1)
      end if
      call mexit(6,0)
      END
C
C=====================================================================
C
      SUBROUTINE GETNAM0(NATOM,NRES,MBONA,NF)
c
      CHARACTER*80 FMT,FMTIN,IFMT,AFMT,RFMT,TYPE
      character(len=4) ITITL(20)
      IFMT = '(12I6)'
      AFMT = '(20A4)'
      RFMT = '(5E16.8)'
C
C     ----- READ THE MOLECULAR TOPOLOGY -----
C
      FMTIN = AFMT
      TYPE = 'TITLE'
      CALL NXTSEC(NF,  0,  0,FMTIN,  TYPE,  FMT,  IOK)
      READ(NF,FMT) (ITITL(I),I=1,20)
c
      FMTIN = IFMT
      TYPE = 'POINTERS'
      CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
      READ(NF,FMT) NATOM,NTYPES,NBONH,MBONA,NTHETH,MTHETA,NPHIH,MPHIA,
     1              NHPARM,NPARM,NNB,NRES
      !---- fix for allocating bond arrays ------
      mbona = mbona+nbonh

      RETURN
      END
C
C=====================================================================
C
      SUBROUTINE GETNAM(NATOM,NRES,IGRAPH,IPRES,LBRES,ITITL,NF,KPF,
     1           X,IX,IH,ib,jb,nbond,chg,radius,lastat,ter,hasradii,
     2           residue_number, residue_chainid, atom_altloc,
     3           atom_element, residue_icode, ext_pdb_data)
      character(len=4) igraph,lbres,ititl
      logical ter,hasradii
      DIMENSION IGRAPH(*),IPRES(*),LBRES(*),ITITL(20)
      DIMENSION X(*),IX(*),IH(*),ib(*),jb(*),chg(*),lastat(*),radius(*)
      integer :: residue_number(NRES)
      character(len=4) :: residue_chainid(NRES), atom_altloc(NATOM),
     &         atom_element(NATOM), residue_icode(NRES)
      logical :: ext_pdb_data
c
      CHARACTER*80 FMT
      CHARACTER*80 FMTIN,IFMT,AFMT,RFMT,TYPE
      IFMT = '(12I6)'
      AFMT = '(20A4)'
      RFMT = '(5E16.8)'
      hasradii = .false.
C
C     ----- READ THE MOLECULAR TOPOLOGY -----
C
      FMTIN = AFMT
      TYPE = 'TITLE'
      CALL NXTSEC(NF,  0,  0,FMTIN,  TYPE,  FMT,  IOK)
      READ(NF,FMT) (ITITL(I),I=1,20)
c
      FMTIN = IFMT
      TYPE = 'POINTERS'
      CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
      READ(NF,FMT) NATOM,NTYPES,NBONH,MBONA,NTHETH,MTHETA,NPHIH,MPHIA,
     1              NHPARM,NPARM,NNB,NRES,NBONA,NTHETA,NPHIA,
     1              NUMBND,NUMANG,NPTRA,NATYP,NPHB,IFPERT,NBPER,NGPER,
     1              NDPER,MBPER,MGPER,MDPER,IFBOX,NMXRS,IFCAP
C
      NTYPE = NTYPES*NTYPES
C
C     ----- READ THE SYMBOLS AND THE CHARGES AND THE MASSES -----
C
      FMTIN = AFMT
      TYPE = 'ATOM_NAME'
      CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
      READ(NF,FMT) (IGRAPH(I),I = 1,NATOM)
c
      FMTIN = RFMT
      TYPE = 'CHARGE'
      CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
      READ(NF,FMT) (CHG(I),I = 1,NATOM)
      do i=1,natom
        chg(i) = chg(i)/18.2223
      end do
c
      if( iok.eq.-1 ) then   ! this is an old-style prmtop file
c
        READ(NF,40) (X(I),I = 1,NATOM)
        READ(NF,30) (IX(I),I = 1,NATOM)
        READ(NF,30) (IX(I),I = 1,NATOM)
        READ(NF,30) (IX(I),I = 1,NTYPE)
        READ(NF,20) (LBRES(I),I=1,NRES)
        READ(NF,30) (IPRES(I),I=1,NRES)
        IPRES(NRES+1) = NATOM+1
   20   FORMAT(20A4)
   30   FORMAT(12I6)
   40   FORMAT(5E16.8)
C
C       ----- READ THE PARAMETERS -----
C
        READ(NF,40) (X(I),    I = 1,NUMBND)
        READ(NF,40) (X(I),   I = 1,NUMBND)
        READ(NF,40) (X(I),    I = 1,NUMANG)
        READ(NF,40) (X(I),   I = 1,NUMANG)
        READ(NF,40) (X(I),    I = 1,NPTRA)
        READ(NF,40) (X(I),    I = 1,NPTRA)
        READ(NF,40) (X(I), I = 1,NPTRA)
        READ(NF,40) (X(I), I = 1,NATYP)
C
        NTTYP = NTYPES*(NTYPES+1)/2
C
        READ(NF,40) (X(I),   I = 1,NTTYP)
        READ(NF,40) (X(I),   I = 1,NTTYP)
C
C     ----- READ THE BONDING INFORMATION -----
C
        READ(NF,30) (IB(I),JB(I),IX(I),I = 1,NBONH)
        READ(NF,30) (IB(I),JB(I),IX(I),I = NBONH+1,NBONH+NBONA)
        nbond = nbonh + nbona
c
      else    !   this is new-style prmtop file:
c
        FMTIN = AFMT
        TYPE = 'RESIDUE_LABEL'
        CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
        READ(NF,FMT) (LBRES(I),I=1,NRES)
c
        FMTIN = IFMT
        TYPE = 'RESIDUE_POINTER'
        CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
        READ(NF,FMT) (IPRES(I),I=1,NRES)
        IPRES(NRES+1) = NATOM+1
c
        FMTIN = IFMT
        TYPE = 'BONDS_INC_HYDROGEN'
        CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
        READ(NF,FMT) (IB(I),JB(I),IX(I),I = 1,NBONH)
c
        FMTIN = IFMT
        TYPE = 'BONDS_WITHOUT_HYDROGEN'
        CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
        READ(NF,FMT) (IB(I),JB(I),IX(I),I = NBONH+1,NBONH+NBONA)
        nbond = nbonh + nbona
c
        FMTIN = RFMT
        TYPE = 'RADII'
        CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
        READ(NF,FMT) (RADIUS(I), I=1,NATOM)
        hasradii = .true.
c
        if (ext_pdb_data) then
           call nxtsec(nf,6,1,'*','RESIDUE_NUMBER',fmt,iok)
           if (iok==0) then
              read(nf,fmt) residue_number
              call nxtsec(nf,6,0,'*','RESIDUE_CHAINID',fmt,iok)
              read(nf,fmt) residue_chainid
              call nxtsec(nf,6,0,'*','ATOM_ELEMENT',fmt,iok)
              read(nf,fmt) atom_element
              call nxtsec(nf,6,1,'*','RESIDUE_ICODE',fmt,iok)
              if (iok==0) then
                 read(nf,fmt) residue_icode
              else
                 residue_icode=' '
              end if
              call nxtsec(nf,6,1,'*','ATOM_ALTLOC',fmt,iok)
              if (iok==0) then
                 read(nf,fmt) atom_altloc
              else
                 atom_altloc=' '
              end if
           else
              write(*,'(A)')
     &             'PRMTOP file did not have extended PDB data.'
              ext_pdb_data = .false.
           end if
        end if
c
c
      endif
c
      if( ifbox.gt.0 ) then
c
c   --- prmtop file will have information needed to generate where the
c        TER cards in the output should be:
c
        if( iok.eq.-1 ) then     !  old-style prmtop:
          read(NF,30) (idummy,jdummy,kdummy,ldummy,i=1,ntheth)
          read(NF,30) (idummy,jdummy,kdummy,ldummy,i=1,ntheta)
          read(NF,30) (idummy,jdummy,kdummy,ldummy,mdummy,i=1,nphih)
          read(NF,30) (idummy,jdummy,kdummy,ldummy,mdummy,i=1,nphia)
          read(NF,30) (idummy, i=1,nnb)
          read(NF,40) (xdummy,i=1,nphb)
          read(NF,40) (xdummy,i=1,nphb)
          read(NF,40) (xdummy,i=1,nphb)
          read(NF,20) (idummy,i=1,natom)
          read(NF,20) (idummy,i=1,natom)
          read(NF,30) (idummy,i=1,natom)
          read(NF,30) (idummy,i=1,natom)
        endif
c
c     ----get molecule info to put TER cards in output:
c
        lastat(1) = natom
        FMTIN = IFMT
        TYPE = 'SOLVENT_POINTERS'
        CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
        read(NF,FMT) iptres,nspm,nspsol
c
        FMTIN = IFMT
        TYPE = 'ATOMS_PER_MOLECULE'
        CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
        read(NF,FMT) (lastat(i), i=1,nspm)
        do i=2,nspm
          lastat(i) = lastat(i-1) + lastat(i)
        end do
c
      else if (ter) then
c
c  --- check bonding file to find end of molecules: assume that if there
c       are no bonds between two adjacent residues, there should be a TER
c       card:
c
        imol = 1
        do 50 ires=2,nres
          i1 = ipres(ires-1)
          i2 = ipres(ires) - 1
          j1 = ipres(ires)
          j2 = ipres(ires+1) - 1
          i13 = 3*(i1-1)
          i23 = 3*(i2-1)
          j13 = 3*(j1-1)
          j23 = 3*(j2-1)
          do ibond=1,nbond
            if(    ib(ibond).ge.i13 .and. ib(ibond).le.i23
     +       .and. jb(ibond).ge.j13 .and. jb(ibond).le.j23 ) go to 50
            if(    jb(ibond).ge.i13 .and. jb(ibond).le.i23
     +       .and. ib(ibond).ge.j13 .and. ib(ibond).le.j23 ) go to 50
          end do
c
c  --- here if there is no bond between the two residues:
c
          lastat(imol) = i2
          imol = imol + 1
   50   continue
        lastat(imol) = natom
      else
c
c   ----  user has requested no TER cards:
c
        lastat(1) = 0
      end if
      RETURN
      END
C
C=====================================================================
C
      SUBROUTINE GENPDB(NATOM,NRES,C,IGRAPH,IPRES,LBRES,ITITL,NF,arg1,
     .    alttit,title,center,chg,aatm,bres,ioffset,lastat,
     .    ftype,radius,hasradii,remediate,
     .    residue_number, residue_chainid, atom_altloc,
     .    atom_element, residue_icode, ext_pdb_data)
      integer :: residue_number(NRES)
      character(len=4) :: residue_chainid(NRES), atom_altloc(NATOM),
     &         atom_element(NATOM), residue_icode(NRES)
      logical :: ext_pdb_data

      integer jpat
      character*6 arg1
      character*4 igraph,code,lbres,atnam,tmpnam
      character*3 resnam
      character*1 type, ftype
      character*40 title
      character*23 occb
      character*3  element
      character(len=4) ititl(20)
      logical alttit,center,aatm,bres,hasradii,remediate
      double precision c
      DIMENSION C(*),IGRAPH(*),IPRES(*),LBRES(*),chg(*),radius(*),
     .     lastat(*)
      DIMENSION ftype(*)
      dimension elrad(15)
c
c    from Sitkoff, Sharp & Honig, PARSE values, J. Phys. Chem. 98:1978(1994)
c
c     data elrad/2.0, 1.4, 1.0, 1.5, 1.0, 0.0, 1.7, 1.85, 0.0, 0.0, 0.0,
c    .     0.0, 0.0, 0.0, 1.0/
c
c    Bondi radii, from Huheey, Inorganic Chemistry: Principles of Structure
c        and Reactivity, 2nd ed. p. 232.
c
      data elrad/1.85, 1.5, 1.0, 1.55, 1.0, 0.0, 1.7, 1.8, 0.0, 0.0,
     .     0.0, 0.0, 0.0, 0.0, 1.2/
c
c     dummy occupation and b-factor string:
      data occb/ '  1.00  0.00           '/
C
c   ---isrn =2 makes this a "normal" atom for ms
      isrn = 2
      nn = 1
      if (alttit) then
        write(nf,100) title
      else
        WRITE(NF,90) (ITITL(i),i=1,15)
      end if
c
c  -- if -ctr was requested, need to center molecule at origin
c
      if (center) then
        xs = 0.0
        ys = 0.0
        zs = 0.0
        k = 0
        do iat=1,natom
          xs = xs + c(k+1)
          ys = ys + c(k+2)
          zs = zs + c(k+3)
          k = k + 3
        end do
        xs = xs /natom
        ys = ys /natom
        zs = zs /natom
        k = 0
        do iat=1,natom
          c(k+1) = c(k+1) - xs
          c(k+2) = c(k+2) - ys
          c(k+3) = c(k+3) - zs
          k = k + 3
        end do
      end if
      imol = 1
      DO J = 1,NRES
        J1 = IPRES(J)
        J2 = IPRES(J+1)-1
        if (bres) then
c
c       ---convert protein residue names back to more like Brookhaven format:
c
          if (lbres(j).eq.'HID ' .or. lbres(j).eq.'HIE ' .or.
     .     lbres(j).eq.'HIP '.or. lbres(j).eq.'HIC') lbres(j) = 'HIS '
          if (lbres(j).eq.'CYX ') lbres(j) = 'CYS '
          if (lbres(j).eq.'CYM ') lbres(j) = 'CYS '
          if (lbres(j).eq.'MEM ') lbres(j) = 'MET '
          if (lbres(j).eq.'ASH ') lbres(j) = 'ASP '
          if (lbres(j).eq.'GLH ') lbres(j) = 'GLU '
c
c       ---also for nucleic acid names:
c
          if( lbres(j)(1:2).eq.'RG' ) lbres(j) = '  G '
          if( lbres(j)(1:2).eq.'DG' ) lbres(j) = ' DG '
          if( lbres(j)(1:2).eq.'RC' ) lbres(j) = '  C '
          if( lbres(j)(1:2).eq.'DC' ) lbres(j) = ' DC '
          if( lbres(j)(1:2).eq.'RA' ) lbres(j) = '  A '
          if( lbres(j)(1:2).eq.'DA' ) lbres(j) = ' DA '
          if( lbres(j)(1:2).eq.'RU' ) lbres(j) = '  U '
          if( lbres(j)(1:2).eq.'DT' ) lbres(j) = ' DT '
c
        end if
        DO K = J1,J2
          if ( aatm ) then
            write(atnam,'(a4)') igraph(k)
          else if ( remediate ) then
c
c        ---convert atom names to closely resemble those used by the
c           wwPDB in its version 3 (aka "remdiated") files:
c
c        ---First, assume that there are no two-character element names
c           (like Fe or Ca or Na).  Then, according to Brookhaven rules,
c           column 13 will be blank, and the name will be left-justified
c           starting in column 14.  UNLESS, the name is four characters
c           long!  In that case, don't use the first blank.
c
            resnam = lbres(j)(1:3)
            write(tmpnam,'(a4)') igraph(k)
            element(1:1) = tmpnam(1:1)
            element(2:3) = '  '
            if (tmpnam(4:4) .eq. ' ') then 
               atnam(1:1) = ' '
               atnam(2:4) = tmpnam(1:3)
            else
               atnam(1:4) = tmpnam(1:4)
            endif
c
c       --- Special fixes where old Amber nucleic acid atom names differ from
c           version 3 pdb names:
c           (N.B.: this little section is no longer necessary if ff09 is used)
c
            if( atnam(1:4) .eq. 'H5''1' ) atnam(1:4) = ' H5'''
            if( atnam(1:4) .eq. 'H5''2' ) atnam(1:4) = 'H5'''''
            if( atnam(1:4) .eq. 'H2''1' ) atnam(1:4) = ' H2'''
            if( atnam(1:4) .eq. 'H2''2' ) atnam(1:4) = 'H2'''''
            if( atnam(1:4) .eq. ' O1P' ) atnam(1:4) = ' OP1'
            if( atnam(1:4) .eq. ' O2P' ) atnam(1:4) = ' OP2'
            if( atnam(1:4) .eq. ' H5T' ) atnam(1:4) = 'HO5'''
            if( atnam(1:4) .eq. ' H3T' ) atnam(1:4) = 'HO3'''
            if( atnam(1:4) .eq. 'HO''2' ) atnam(1:4) = 'HO2'''
c
c       --- Now, special case out the two-character element names:
c
            if( atnam(1:4).eq.' Na+' .or. atnam(1:4).eq.' NA+' .or.
     +          atnam(1:3).eq.' Fe' .or. atnam(1:3).eq.' FE' .or.
     +          atnam(1:3).eq.' Cl' .or. atnam(1:3).eq.' CL' .or.
     +          atnam(1:3).eq.' Zn' .or. atnam(1:3).eq.' ZN' .or.
     +          atnam(1:4).eq.' Li+' .or. atnam(1:4).eq.' LI+' .or.
     +          atnam(1:4).eq.' Ca+' .or. atnam(1:4).eq.' CA+' .or.
     +          atnam(1:4).eq.' Mg+' .or. atnam(1:4).eq.' MG+' .or.
     +          atnam(1:4).eq.' Br-' .or. atnam(1:4).eq.' BR-' ) then
              atnam(1:1) = atnam(2:2)
              atnam(2:2) = atnam(3:3)
              atnam(3:3) = atnam(4:4)
              atnam(4:4) = ' '
            end if
c
          else
c
c        ---convert atom names to closely resemble those used by Brookhaven
c           *before* the "remediation" that changed hydrogen names
c
c        ---First, assume that there are no two-character element names
c           (like Fe or Ca or Na).  Then, according to Brookhaven rules,
c           column 13 will be blank, and the name will be left-justified
c           starting in column 14.  UNLESS, the name is four characters
c           long!  In that case, wrap around the final character into
c           column 13.
c
            atnam = '    '
            resnam = lbres(j)(1:3)
            write(tmpnam,'(a4)') igraph(k)
            atnam(2:4) = tmpnam(1:3)
            if (tmpnam(4:4) .ne. ' ') atnam(1:1) = tmpnam(4:4)
c           write(6,*) 'converting ', tmpnam, '->', atnam
c
c       --- here are some more Brookhaven wraparounds:
c           This gives files that look very much like Brookhaven, EXCEPT
c           that Brookhaven uses "1" and "2" for beta-protons (for example)
c           whereas the standard Amber database (along with many in
c           the NMR field) use "2" and "3", i.e. we have 2HB and 3HB,
c           whereas Brookhaven files use 1HB and 2HB.
c
            if (atnam(1:2) .eq. ' H' .and.
     .          (atnam(4:4) .eq. '1' .or. atnam(4:4).eq.'2' .or.
     .           atnam(4:4) .eq. '3')) then
              if (atnam(2:3) .eq. 'HB' .or. atnam(2:3) .eq. 'H7' .or.
     .            atnam(2:3) .eq. 'H6' ) then
                atnam(1:1) = atnam(4:4)
                atnam(4:4) = ' '
              end if
              if (atnam(2:3) .eq. 'HG' .and. resnam.ne.'THR') then
                atnam(1:1) = atnam(4:4)
                atnam(4:4) = ' '
              end if
              if (atnam(2:3) .eq. 'HD' .and. (resnam .ne. 'PHE'
     .               .and. resnam .ne. 'TYR' .and. resnam .ne. 'TRP'
     .               .and. resnam(1:2) .ne. 'HI')) then
                atnam(1:1) = atnam(4:4)
                atnam(4:4) = ' '
              end if
              if (atnam(2:3) .eq. 'HE' .and. (resnam .ne. 'PHE'
     .               .and. resnam .ne. 'TYR' .and. resnam .ne. 'TRP'
     .               .and. resnam(1:2) .ne. 'HI')) then
                atnam(1:1) = atnam(4:4)
                atnam(4:4) = ' '
              end if
            end if
c
            if (atnam.eq.' H1 ') atnam = '1H  '
            if (atnam.eq.' H2 ') then
               if( resnam .ne. 'ADE' .and. resnam(1:2) .ne. 'DA' .and.
     .             resnam(1:2) .ne. 'RA' ) atnam = '2H  '
            end if
            if (atnam.eq.' H3 ') then
               if( resnam .ne. 'THY' .and. resnam(1:2) .ne. 'DT' .and.
     .             resnam(1:2) .ne. 'RU' .and. resnam .ne. 'URA' ) 
     .             atnam = '3H  '
            end if
c
c       --- Convert nucleic acid primed names into asterisk: these
c           are always in the fourth column of the atom name:
c
            if (atnam(4:4).eq.'''') atnam(4:4) = '*'
c
c       --- Now, special case out the two-character element names:
c
            if( atnam(1:4).eq.' Na+' .or. atnam(1:4).eq.' NA+' .or.
     +          atnam(1:3).eq.' Fe' .or. atnam(1:3).eq.' FE' .or.
     +          atnam(1:3).eq.' Cl' .or. atnam(1:3).eq.' CL' .or.
     +          atnam(1:3).eq.' Zn' .or. atnam(1:3).eq.' ZN' .or.
     +          atnam(1:4).eq.' Li+' .or. atnam(1:4).eq.' LI+' .or.
     +          atnam(1:4).eq.' Ca+' .or. atnam(1:4).eq.' CA+' .or.
     +          atnam(1:4).eq.' Mg+' .or. atnam(1:4).eq.' MG+' .or.
     +          atnam(1:4).eq.' Br-' .or. atnam(1:4).eq.' BR-' ) then
              atnam(1:1) = atnam(2:2)
              atnam(2:2) = atnam(3:3)
              atnam(3:3) = atnam(4:4)
              atnam(4:4) = ' '
              element(1:2) = atnam(1:2)
            end if
c
          end if
c
          JC = 3*K-3
          J_print = mod(J+ioffset,10000)
          if (arg1.eq.'-pdb') then
            if (ext_pdb_data) then
               write(nf,
     &         '(A4,I7,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,2A4)'),
     &             'ATOM',k,atnam,atom_altloc(j)(1:1),lbres(j),
     &             residue_chainid(j)(1:1),residue_number(j),
     &             residue_icode(j)(1:1),C(jc+1:jc+3),1.0,0.0,
     &             atom_element(k)
            else if( remediate ) then
               WRITE(NF,60) K,atnam,LBRES(J),J_print,(C(JC+M),M=1,3), 
     .                   occb,element
            else
               WRITE(NF,61) K,atnam,LBRES(J),J_print,(C(JC+M),M=1,3) 
            end if
          else if (arg1.eq.'-atm' .or. arg1.eq.'-pqr' 
     .                            .or. arg1.eq.'-sas') then
            code = igraph(k)
            type = code(1:1)
            if (type.eq.'H') then
              itype = 15
            else if (type.eq.'C') then
              itype = 7
            else if (type.eq.'N') then
              itype = 4
            else if (type.eq.'O') then
              itype = 2
            else if (type.eq.'P') then
              itype = 1
            else if (type.eq.'S') then
              itype = 8
            else if (type.eq.'L') then
              itype = 16
            else if (type.eq.'Z') then
              itype = 5
            else
              itype = 3
            end if
            if (arg1.eq.'-atm') then
              write(nf,50) (c(jc+m),m=1,3),
     .          itype,isrn,nn,lbres(j),j+ioffset,igraph(k)
            else if (arg1.eq.'-pqr') then
              if (hasradii) then
                WRITE(NF,81) K,atnam,LBRES(J),J_print,
     .             (C(JC+M),M=1,3),chg(k),radius(k)
              else
                WRITE(NF,81) K,atnam,LBRES(J),J_print,
     .             (C(JC+M),M=1,3),chg(k),elrad(itype)
              endif
            else if (arg1.eq.'-sas') then
              WRITE(NF,81) K,atnam,LBRES(J),J_print,
     .           (C(JC+M),M=1,3),chg(k),elrad(itype)+1.4
            end if
          else if(arg1.eq.'-first') then
            if(atnam.eq.' X' .and. LBRES(J).eq.'BMH') then
              WRITE(NF,67) K,atnam,LBRES(J),J_print - jpat,
     .                     (C(JC+M),M=1,3), 0.00, 0.00, 
     .                     J_print - jpat, ftype(k)
            else
              WRITE(NF,65) K,atnam,LBRES(J),J_print,(C(JC+M),M=1,3),
     .                     0.00, 0.00, J_print, ftype(k)
              jpat = J_print
            endif
          end if
          if( k.eq.lastat(imol)) then
            write(NF,120)
            imol = imol + 1
          end if
        end do
      end do
C
      write(NF,130)
C
      RETURN
   50 FORMAT(3f10.3,3i5,5x,a4,i5,1x,a4)
   60 FORMAT('ATOM',2X,I5,1X,A4,1X,A4,1X, I4,4X,3F8.3,A23,A3)
   61 FORMAT('ATOM',2X,I5,1X,A4,1X,A4,1X, I4,4X,3F8.3)
   65 FORMAT('ATOM',2X,I5,1X,A4,1X,A4,1X, I4,4X,3F8.3,2F6.2,I5,1X,A1)
   67 FORMAT('HETATM', I5,1X,A4,1X,A4,'G',I4,4X,3F8.3,2F6.2,I5,1X,A1)
   70 FORMAT('ATOM',2X,I5,1X,A4,1X,A4,1X, I4,4X,3F8.3,16x,f5.3)
   80 FORMAT('ATOM',2X,I5,1X,A4,1X,A4,1X, I4,4X,3F8.3,2f10.5)
   81 FORMAT('ATOM',2X,I5,1X,A4,1X,A4,1X, I4,4X,3F8.3,f8.4,f8.3)
   90 FORMAT('REMARK  ',15a4)
  100 format('REMARK  ',a40)
  110 FORMAT('REMARK  ',10I5)
  120 FORMAT('TER   ')
  130 FORMAT('END   ')
      END
C
C=====================================================================
C
      SUBROUTINE GETCOR(bin,NATOM,C,CD,KCF,IDBL,NF)
      character(len=4) ititl(20)
      double precision C(*),CD(*)
      integer fh
      character*80 line
      logical bin
C
      NAT3 = 3*NATOM
C
      if(.not. bin) then
        READ(NF,30) ITITL
        read(nf,'(a)') line
        if( line(6:6).eq.' ' ) then
          READ(line,'(i5)') MATOM
        else
          READ(line,'(i6)') MATOM
        end if
        IF(MATOM.NE.NATOM) then
          write(0,20) natom,matom
          stop
        end if
        READ(NF,50) (C(I),I=1,NAT3)
        RETURN
      else
        read(nf) ititl
        read(nf) matom
        IF(MATOM.NE.NATOM) then
          write(0,20) natom,matom
          stop
        end if
        read(nf, end=15, err=15)(c(i), i=1,nat3)
        return       
      endif
c
   15 continue
      write(0,'(a,a)') 'COULD NOT READ COORDINATES FROM RESTRT FILE'
      stop
c
   20 FORMAT(/2X,'ATOMS DO NOT MATCH BETWEEN PARM AND MIN FILES',2i8)
   30 FORMAT(20A4)
   50 FORMAT(6F12.7)
c
      END
C
C=====================================================================
C
      subroutine cpkgen(x,natom,igraph)
      character*4 igraph(*),atnam
      double precision x(*)
c
      character*6 color
      character*1 code
      character*2 code2
c
      k = 0
      do i=1,natom
c
c---codes for labelling atoms by atom type:
c
        atnam = igraph(i)
        code = atnam(1:1)
        if (code.eq.'H') then
          rad = 1.1
          color = 'white'
        else if (code.eq.'C') then
          rad = 1.7
          color = 'green'
        else if (code.eq.'N' .or. code.eq.'P') then
          rad = 1.6
          color = 'blue '
        else if (code.eq.'O') then
          rad = 1.6
          color = 'red  '
        else if (code.eq.'S') then
          rad = 2.0
          color = 'black'
        else if (code.eq.'L') then
          rad = 1.6
          color = 'white'
        else if (code.eq.'F') then
          rad = 1.8
          color = 'red  '
        else
          write (0,'(a,a1)') 'bad code:', code
          stop
        end if
c
c--- special section to give small, blue water molecules
c
        code2 = atnam(1:2)
        if (code2.eq.'OW') then
          rad = 0.8
          color = 'blue '
        else if (code2.eq.'HW') then
          rad = 0.4
          color = 'blue '
        end if
c
c---rotate by some amount about x-axis
c
c       xp =  x(k+1)
c       yp = ct*x(k+2) - st*x(k+3)
c       zp = st*x(k+2) + ct*x(k+3)
c
c---followed by rotation about new y-axis
c
c       xp2 = ct2*xp - st2*zp
c       yp2 = yp
c       zp2 = st2*xp + ct2*zp
c
        if (code.ne.'L')
     .    write(6,'(3f8.3,f5.2,1x,a6)') x(k+1),x(k+2),x(k+3),rad,color
        k = k + 3
      end do
      return
      end
C
C=====================================================================
C
      subroutine pluto(x,natom,igraph,alttit,title)
      character*1 blnk,c1
      character*4 igraph,atnam
      character*40 title
      double precision x
      dimension x(*),igraph(*)
      logical alttit
c
c     ------ program to write pluto input file -----
c
      write(6,120)
      blnk = ' '
      k = 0
      do j=1,natom
         atnam = igraph(j)
         c1 = atnam(1:1)
         if (j.le.9) then
           write(6,130) c1,j,x(k+1),x(k+2),x(k+3)
         else if (j.le.99) then
           write(6,140) c1,j,x(k+1),x(k+2),x(k+3)
         else
           write(6,150) c1,j,x(k+1),x(k+2),x(k+3)
         end if
         k = k + 3
      end do
      write(6,30)
      write(6,20) title
      write(6,50)
      write(6,90)
      write(6,100)
      write(6,110)
      write(6,60)
      write(6,70)
      write(6,80)
      return
   20 format('TITLE',1x,a40)
   30 format('*')
   40 format('LABEL',1X,A40)
   50 format('JOIN RADII C 0.8 N 0.8 H 0.4 O 0.74 P 1.1 S 1.1')
   60 format('STEREO')
   70 format('VIEW YORIGIN')
   80 format('PLOT')
   90 format('SOLID')
  100 format('RADII ATOMS C 0.3 N 0.4 O 0.4 H 0.1')
  110 format('RADII BONDS 0.05 8 TAPER 5')
  120 format('DATA',1x,a40)
  130 format(a1,i1,t9,3f8.3)
  140 format(a1,i2,t9,3f8.3)
  150 format(a1,i3,t9,3f8.3)
      end
C
C=====================================================================
C
      subroutine usage()
       write(*,'(A)')
     & 'Usage:',
     & 'ambpdb [OPTION]... < restrt > out.pdb',
     & '',
     & 'Options:',
     & ' -p PRMTOP     Define PRMTOP filename (default:"prmtop").',
     & ' -tit TITLE    Write a REMARK record containing TITLE.',
     & '                   (default: use prmtop title)',
     & ' -aatm         Left-justified Amber atom names.',
     & ' -bres         Brookhaven Residue names (HIE->HIS, etc.).',
     & ' -ctr          Center molecule on (0,0,0).',
     & ' -noter        Do not write TER records.',
     & ' -ext          Use PRMTOP extended PDB info, if present.',
     & ' -ene FLOAT    Define H-bond energy cutoff for FIRST.',
     & ' -bin          The coordinate file is in binary form.',
     & ' -offset INT   Add offset to residue numbers.',
     & '',
     & 'Options for alternate output format (give only one of these):',
     & ' -pqr          PQR (MEAD) format with charges and radii.',
     & ' -sas          PQR with 1.4 added to atom radii.',
     & ' -bnd          list bonds from the PRMTOP.',
     & ' -atm          Mike Connolly surface/volume format.',
     & ' -first        Add REMARKs for input to FIRST.'
c
c    & 'Disabled/obsolete output options:'
c    & ' -cpk          X,Y,X,R,color'
c    & ' -plu          Write a PLUTO plot file'
       call mexit(6,0)
      end subroutine usage

