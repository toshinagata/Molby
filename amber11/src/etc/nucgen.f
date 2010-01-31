C
C************************************************************************
C                              AMBER                                   **
C                                                                      **
C                  Copyright (c) 1986, 1991, 1995                      **
C             Regents of the University of California                  **
C                       All Rights Reserved.                           ** 
C                                                                      **
C  This software provided pursuant to a license agreement containing   **
C  restrictions on its disclosure, duplication, and use. This software **
C  contains confidential and proprietary information, and may not be   **
C  extracted or distributed, in whole or in part, for any purpose      **
C  whatsoever, without the express written permission of the authors.  **
C  This notice, and the associated author list, must be attached to    **
C  all copies, or extracts, of this software. Any additional           **
C  restrictions set forth in the license agreement also apply to this  **
C  software.                                                           **
C************************************************************************
C
      PROGRAM NUCGEN
c   
c    copyright @1985,1989 Regents of the University of California
c    Rev A revision: George Seibel, 1989.
c    AUTHORS: U.C.Singh, N.Pattabiraman, and S.N.Rao, 1985.
c    Director: P.A. Kollman
C
C     ----- PROGRAM TO GENERATE DOUBLE HELICAL DNA OF DIFFERENT
C           HELICAL CONFORMATION -----
C
      common /files/ ngin, ngout, ngdat, pdbout, owrite
      character*80   ngin, ngout, ngdat, pdbout
      character(len=4) ititl(20),kend,iblank,ilbmol
      character owrite
      COMMON/IOFILE/IN,IOUT,IOUTB
      character(len=4) LBRES(1000),LBDUM(20)
      DATA KEND,IBLANK/'END ','    '/
c
c     --- get file names ---
c
      call ngfil
C
C     --- open all files ---
C
      IN = 5
      IOUT = 6
      IOUTB = 10
c     subr amopen(lun,fname,fstat,fform,facc)
c     -- output
      call amopen(iout, ngout, owrite,'F','W')
c     -- input
      call amopen(in,   ngin,  'O','F','R')
c     -- dna data
      call amopen(7,    ngdat, 'O','F','R')
c     -- pdb output
      call amopen(ioutb,pdbout,owrite,'F','W')
C
C     ----- READ THE RESIDUE INFORMATION UNTIL AN END CARD
C           IS FOUND -----
C
      WRITE(IOUT,9108)
      WRITE(IOUT,9118)
C
      NMOL = 0
      NRES = 0
  100 CONTINUE
      NMOL = NMOL+1
      READ(IN,9008) ITITL
      IF(ITITL(1).EQ.KEND) GO TO 200
      WRITE(IOUT,9008) ITITL
      READ(IN,9008) ILBMOL
      IF(NMOL.GT.2) GO TO 900
      NRESM = 0
C
  120 CONTINUE
        READ(IN,9018) (LBDUM(K),K=1,16)
        IF(LBDUM(1).EQ.IBLANK) GO TO 180
        DO 140 K = 1,16
          IF(LBDUM(K).EQ.IBLANK) GO TO 160
          NRES = NRES+1
          NRESM = NRESM+1
          LBRES(NRES) = LBDUM(K)
  140   CONTINUE
  160   CONTINUE
      GO TO 120
  180 CONTINUE
C
C     ----- OUTPUT THE MOLECULE INFORMATION -----
C
      WRITE(IOUT,9128) NMOL,NRESM
      WRITE(IOUT,9158) (LBRES(K),K=nres-nresm+1,nres)
C
      GO TO 100
  200 CONTINUE
C
C     ----- NOW CALL THE GENERATION ROUTINE -----
C
      CALL GENNUC(NRES,LBRES,ILBMOL)
C
      call mexit(0, 0)
  900 CONTINUE
      WRITE(IOUT,9208)
      call mexit(iout, 1)
 9008 FORMAT(20A4)
 9018 FORMAT(16(A4,1X))
 9108 FORMAT(/10X,54(1H-),/10X,
     +       'Amber 5   NUCGEN                           UCSF 1997',
     +       /10X,54(1H-))
 9118 FORMAT(/15X,'INPUT MOLECULES INFORMATION',/)
 9128 FORMAT(/5X,'MOLECULE ',I2,' CONTAINS ',I4,' RESIDUES:',/)
 9158 FORMAT(5X,12(A4,1X))
 9208 FORMAT(/5X,'INPUT CONTAINS MORE THAN TWO MOLECULES AND ASSUMED',
     +       /5X,'THAT IT IS NOT A STANDARD DNA OR RNA',/)
      END
C
C************************************************************************
C                              AMBER                                   **
C                                                                      **
C                  Copyright (c) 1986, 1991, 1995                      **
C             Regents of the University of California                  **
C                       All Rights Reserved.                           ** 
C                                                                      **
C  This software provided pursuant to a license agreement containing   **
C  restrictions on its disclosure, duplication, and use. This software **
C  contains confidential and proprietary information, and may not be   **
C  extracted or distributed, in whole or in part, for any purpose      **
C  whatsoever, without the express written permission of the authors.  **
C  This notice, and the associated author list, must be attached to    **
C  all copies, or extracts, of this software. Any additional           **
C  restrictions set forth in the license agreement also apply to this  **
C  software.                                                           **
C************************************************************************
C
      SUBROUTINE GENNUC(NRES,LBRES,ILBMOL)
      character*8 NAME,NAMA,JTYPM,ISPEC
      character*4 tmpnam, atnam, lbres, idna, irna, jseq, katn, 
     .            natn, kseq, ilbmol, namat
C
C     ----- ROUTINE TO GENERATE STANDARD DNA AND RNA DOUBLE
C           HELICES. SIX TYPES ARE AVAILABLE -----
C
C           ARNA    ...   RIGHT HANDED A-RNA (ARNOTT)
C           APRNA   ...   RIGHT HANDED A PRIME RNA (ARNOTT)
C           BDNAL   ...   RIGHT HANDED BDNA (LANGRIDGE)
C           BDNAA   ...   RIGHT HANDED BDNA (ARNOTT)
C           SDNA    ...   LEFT HANDED BDNA (SASISEKHARAN)
C           ADNA    ...   RIGHT HANDED A-DNA (ARNOTT)
C
      COMMON/IOFILE/IN,IOUT,IOUTB
C
      DIMENSION R(60),PHI(60),ZZ(60),NAMAT(60)
      DIMENSION IR1(9),IR2(9),JSEQ(9),HXREP(6),HXHT(6)
      DIMENSION KATN(3),NATN(3)
      DIMENSION NAME(6),IRDS(9,2),IRDE(9,2)
      DIMENSION LBRES(*),IPMOL(500)
      DIMENSION KSEQ(15),INDEXD(5,15),INDEXU(5,15)
C     
      DATA NAME/'$ARNA'   ,'$APRNA'  ,'$LBDNA'  ,'$ABDNA'  ,
     +          '$SBDNA'  ,'$ADNA'   /
      DATA ISPEC/'$SPECIAL'/
      DATA IDNA,IRNA/'D   ','R   '/
      DATA JSEQ/'HB  ','HE  ','POM ','RIBO','ADE ','GUA ','THY ',
     +          'CYT ','URA '/
      DATA KATN/'OA  ','OB  ','O1  '/
      DATA NATN/'O1P ','O2P ','O4  '/
      DATA IRDS/   1,   2,   3,   6,  14,  26,  40,  50,   0,
     +             1,   2,   3,   6,  15,  27,   0,  41,  51/
      DATA IRDE/   1,   2,   5,  13,  25,  39,  49,  59,   0,
     +             1,   2,   5,  14,  26,  40,   0,  50,  59/
      DATA HXREP/ 32.7, 30.0, 36.0, 36.0, 36.0, 32.7/
      DATA HXHT / 2.81, 3.00, 3.38, 3.38,-3.38, 2.56/
      DATA KSEQ/'A5  ','A   ','A3  ','G5  ','G   ','G3  ','T5  ',
     +          'T   ','T3  ','C5  ','C   ','C3  ','U5  ','U   ',
     +          'U3  '/
C
C     ----- EVALUATE SOME CONSTANTS -----
C
      PI = 4.0E0*ATAN(1.0E0)
      CONV = PI/180.0E0
      IATOM = 0
      NBASP = 0
      IACID = 0
      NF = 7
C
C     ----- SET THE CORRESPONDENCE BETWEEN LBRES AND
C           IPMOL -----
C
      CALL SETSEQ(NRES,LBRES,IPMOL,JSEQ,KSEQ,inew)
      if (inew.eq.1)
     .  write(iout, '(4x, a)') 'New (1994) residue naming convention'
C
C     ----- READ THE TYPE OF NUCLEIC ACID NEEDED -----
C
      READ(IN,9008) JTYPM
C
C     ----- IACID IS 1 FOR DNA OR 2 FOR RNA -----
C
      IF(ILBMOL.EQ.IDNA) IACID = 1
      IF(ILBMOL.EQ.IRNA) IACID = 2
      IF(IACID.EQ.0) GO TO 900
C
C     ----- TRANSFERING THE NUMBER OF ATOMS FOR DEOXY OR RIBOSE SUGAR
C           AND BASES -----
C
      DO 120 I = 1,9
        IR1(I) = IRDS(I,IACID)
        IR2(I) = IRDE(I,IACID)
  120 CONTINUE
      IF (inew.eq.0) GO TO 139
C
C     --- SET UP NEW STYLE OF NUCLEIC ACIDS RESIDUES NOMENCLATURE 
C
      DO 125 I = 1,5
        DO 125 J = 1,12
          INDEXD(I,J) = 0
          INDEXU(I,J) = 0
  125 CONTINUE
C
C     ----- FIND NEW ATOM INDEXES FOR a) HB, b) HE, c) SUGARS, -----
C           d) POM and e) BASES 
C     HB
      DO 126 J = 1,15,3
        INDEXD(1,J) = IR1(1)
        INDEXU(1,J) = IR2(1)
  126 CONTINUE
C     HE
      DO 127 J = 3,15,3
        INDEXD(5,J) = IR1(2)
        INDEXU(5,J) = IR2(2)
  127 CONTINUE
C     POM
      DO 128 J = 2,15,3
        INDEXD(2,J) = IR1(3)
        INDEXU(2,J) = IR2(3)
  128 CONTINUE
      DO 129 J = 3,15,3
        INDEXD(2,J) = IR1(3)
        INDEXU(2,J) = IR2(3)
  129 CONTINUE
C     SUGAR
      DO 130 J = 1,15
        INDEXD(3,J) = IR1(4)
        INDEXU(3,J) = IR2(4)
  130 CONTINUE
C     BASES
      J = 1
      DO 132 K = 5,9
        DO 131 I = 1,3
          INDEXD(4,J) = IR1(K)
          INDEXU(4,J) = IR2(K)
          J = J + 1
  131   CONTINUE
  132 CONTINUE
c     do 135 i=1,5
c     write(6,1010) (indexd(i,j),j=1,15) 
c     write(6,1010) (indexu(i,j),j=1,15) 
c 135 continue
c1010 format(1x,15i3)
C 
C
C     ----- BRINGOUT THE UNIT TWIST AND HEIGHT FROM THE DATA -----
C
  139 CONTINUE
      DO 140 I = 1,6
        IF(JTYPM.EQ.NAME(I)) GO TO 180
  140 CONTINUE
      IF(JTYPM.EQ.ISPEC) GO TO 220
      write(iout,*) 'jtypm .ne. ispec: ', JTYPM, ISPEC
      GO TO 900
  180 CONTINUE
C
      IF(JTYPM.EQ.NAME(1)) WRITE(6,9218)
      IF(JTYPM.EQ.NAME(2)) WRITE(6,9228)
      IF(JTYPM.EQ.NAME(3)) WRITE(6,9238)
      IF(JTYPM.EQ.NAME(4)) WRITE(6,9248)
      IF(JTYPM.EQ.NAME(5)) WRITE(6,9258)
      IF(JTYPM.EQ.NAME(6)) WRITE(6,9268)
      UROT = HXREP(I)
      UHIGH = HXHT(I)
      GO TO 240
C
  220 CONTINUE
      READ(IN,9018) UROT,UHIGH
      WRITE(IOUT,9108) UROT,UHIGH
C
C     ----- NOW THE HELICAL PARAMETERS ARE AVAILABLE -----
C
  240 CONTINUE
C
C     ----- READING THE MONOMER UNITS' R ,PHI, ZZ COORDINATES -----
C
      write(iout,'(2x,a)') 'reading monomer parameters'
      IF (inew.EQ.0) THEN
  260   CONTINUE
          READ(NF,9008,END=900) NAMA
          IF(NAMA.EQ.JTYPM) GO TO 270
          GO TO 260
  270   CONTINUE
      ELSE
C
C       ----- NEW-AMBER ATOM NAMES TO BE USED -------
C
  272   CONTINUE
c         write(iout,*)' searching for ', JTYPM
          READ(NF,9008,END=900) NAMA
          IF(NAMA.EQ.JTYPM) GO TO 275
        GO TO 272
  275 CONTINUE
      ENDIF
C
      DO 280 I = 1,60
        READ(NF,9028)  NAMAT(I),R(I),PHI(I),ZZ(I)
  280 CONTINUE
C
      ROT = 0.0
      HIGH = 0.0
      IRES = 0
      NRESH = NRES/2
C
C     ----- GENERATE THE COORDINATES OF THE TWO CHAINS -----
C
      HXMUL = -1.0E0
C
      IF (inew.eq.0) THEN
C
C       ---- OLD DNA/RNA NOMENCLATURE ------------
C
        DO 400 ICHAIN = 1,2
          IF(ICHAIN.EQ.2) HXMUL = 1.0E0
C
          DO 420 I = 1,NRESH
            IF(ICHAIN.EQ.1) NNSEQ = IPMOL(I)
            IF(ICHAIN.EQ.2) NNSEQ = IPMOL(I+NRESH)
            IRES = IRES+1
C
C           ----- IF THE RESIDUE IS A POM,HB OR HE SKIP THE SUGAR -----
C
            IF(NNSEQ.LT.4) GO TO 320
C
C           ----- GENERATE THE SUGAR BACKBONE FIRST -----
C
            IDOI = IR1(4)
            IDOF = IR2(4)
C
            DO 440 IP = IDOI,IDOF
              XRAD = R(IP)
              YYR = (HXMUL*PHI(IP)+ROT)*CONV
              X = XRAD*COS(YYR)
              Y = XRAD*SIN(YYR)
              Z = HXMUL*ZZ(IP)+HIGH
              IATOM = IATOM+1
c
c             --- 3 char names start col. 14, 
c                 4 char names start col 13. ---
c
              atnam = ' '
              tmpnam = namat(ip)
              if (tmpnam(4:4) .eq. ' ') then
                atnam(2:4) = tmpnam
              else
                atnam = tmpnam
              endif
              WRITE(IOUTB,9118) IATOM,atnam,JSEQ(NNSEQ),IRES,X,Y,Z
  440       CONTINUE
C
C           ----- GENERATE THE BASE OR POM ETC COORDINATES -----
C
  320       CONTINUE
            IDOI = IR1(NNSEQ)
            IDOF = IR2(NNSEQ)
            DO 460 IG = IDOI,IDOF
                 YYR = (HXMUL*PHI(IG)+ROT)*CONV
                 XRAD = R(IG)
                 X = XRAD*COS(YYR)
                 Y = XRAD*SIN(YYR)
                 Z = HXMUL*ZZ(IG)+HIGH
                 IATOM = IATOM+1
c
c                --- 3 char names start col. 14, 
c                    4 char names start col 13. ---
c
                 atnam = ' '
                 tmpnam = namat(ig)
                 if (tmpnam(4:4) .eq. ' ') then
                    atnam(2:4) = tmpnam
                 else
                    atnam = tmpnam
                 endif
                 WRITE(IOUTB,9118) IATOM,atnam,JSEQ(NNSEQ),IRES,X,Y,Z
  460       CONTINUE
C
C           ----- DO NOT INCREMENT THE HELIX IF HB ETC -----
C
            IF(NNSEQ.LT.4) GO TO 420
            ROT = ROT+UROT
            HIGH = HIGH+UHIGH
C
C           ----- END OF MOLECULE GENERATION -----
C
  420     CONTINUE
C
C         ----- REVERSE THE UNIT HEIGHT AND TWIST IN ORDER TO
C               DESCEND THE CHAIN -----
C
          UROT = -UROT
          ROT = ROT+UROT
          UHIGH = -UHIGH
          HIGH = HIGH+UHIGH
C
C         ----- END OF CHAIN GENERATION -----
C
  400   CONTINUE
C
      ELSE
C   
C       ------- PART USED FOR NEW DNA/RNA NOMENCLATURE ---------
C
        DO 4400 ICHAIN = 1,2
          IF(ICHAIN.EQ.2) HXMUL = 1.0E0
C
          DO 4420 I = 1,NRESH
            IF(ICHAIN.EQ.1) NNSEQ = IPMOL(I)
            IF(ICHAIN.EQ.2) NNSEQ = IPMOL(I+NRESH)
            IRES = IRES+1
C
C           ----- GENERATE NUCLEOSIDE COORDINATES -----
C
            DO 4470 K = 1,5
              IDOWN = INDEXD(K,NNSEQ)
              IUP   = INDEXU(K,NNSEQ)
              IF(IDOWN.EQ.0) GO TO 4470
C
C             --- correct for the last HE atom:
              if(k.eq.5.and.i.eq.nresh) then
c               write(6,1011) k,i
 1011           format(1x,'k,i',2i5)
                high = high + uhigh
                rot = rot + urot
              endif
C
              DO 4460 IG = IDOWN,IUP
                 YYR = (HXMUL*PHI(IG)+ROT)*CONV
                 XRAD = R(IG)
                 X = XRAD*COS(YYR)
                 Y = XRAD*SIN(YYR)
                 Z = HXMUL*ZZ(IG)+HIGH
                 IATOM = IATOM+1
c
c                -- terminal hydrogens (HB HE residues in old ff)
c                   have special goofy names
c
                 if (i.eq.1 .and. k.eq.1 .and. ig.eq.idown) then
                   atnam = ' H5T'
                 else if (i.eq.nresh .and. k.eq.5 .and. ig.eq.iup) then
                   atnam = ' H3T'
                 else
c
c                  --- adjust phosphate & sugar names
c
                   if (namat(ig).eq.katn(1)) then
                     namat(ig) = natn(1)
                   else if (namat(ig).eq.katn(2)) then
                     namat(ig) = natn(2)
                   else if (namat(ig).eq.katn(3)) then
                     namat(ig) = natn(3)
                   endif
c
c                  --- 3 char names start col. 14, 
c                      4 char names start col 13. ---
c
                   atnam = ' '
                   tmpnam = namat(ig)
                   if (tmpnam(4:4) .eq. ' ') then
                      atnam(2:4) = tmpnam
                   else
                      atnam = tmpnam
                   endif
                 endif
                 WRITE(IOUTB,9118) IATOM,atnam,KSEQ(NNSEQ),IRES,X,Y,Z
 4460         CONTINUE
 4470       CONTINUE
C
            ROT = ROT+UROT
            HIGH = HIGH+UHIGH
C
C           ----- END OF MOLECULE GENERATION -----
C
 4420     CONTINUE
C
          ROT = ROT - UROT
          HIGH = HIGH - UHIGH
C
C         ----- REVERSE THE UNIT HEIGHT AND TWIST IN ORDER TO
C               DESCEND THE CHAIN -----
C
          UROT = -UROT
          ROT = ROT+UROT
          UHIGH = -UHIGH
          HIGH = HIGH+UHIGH
C
C         ----- END OF CHAIN GENERATION -----
C
 4400   CONTINUE
      ENDIF
C
      CLOSE(UNIT=NF)
      RETURN
  900 WRITE(IOUT,9408)
      call mexit(iout, 1)
 9008 FORMAT(5A8)
 9018 FORMAT(3F10.4)
 9028 FORMAT(A4,3F10.4)
 9108 FORMAT(/5X,'SPECIAL HELICAL VALUES', 'UROT =',F8.3,
     +        ' UHIGH =',F8.3,/)
 9118 FORMAT('ATOM',2X,I5,1X,A4,1X,A3,2X,I4,4X,3F8.3)
 9218 FORMAT(/10X,'GENERATING RIGHT HANDED A-RNA(ARNOTT)',/)
 9228 FORMAT(/10X,'GENERATING RIGHT HANDED AP-RNA',/)
 9238 FORMAT(/10X,'GENERATING RIGHT HANDED BDNA (LANGRIDGE)',/)
 9248 FORMAT(/10X,'GENERATING RIGHT HANDED BDNA (ARNOTT)',/)
 9258 FORMAT(/10X,'GENERATING LEFT HANDED BDNA (SASISEKHARAN)',/)
 9268 FORMAT(/10X,'GENERATING RIGHT HANDED A-DNA',/)
 9408 FORMAT(/10X,'ERROR IN INPUT ... GENNUC')
      END
C
C************************************************************************
C                              AMBER                                   **
C                                                                      **
C                  Copyright (c) 1986, 1991, 1995                      **
C             Regents of the University of California                  **
C                       All Rights Reserved.                           ** 
C                                                                      **
C  This software provided pursuant to a license agreement containing   **
C  restrictions on its disclosure, duplication, and use. This software **
C  contains confidential and proprietary information, and may not be   **
C  extracted or distributed, in whole or in part, for any purpose      **
C  whatsoever, without the express written permission of the authors.  **
C  This notice, and the associated author list, must be attached to    **
C  all copies, or extracts, of this software. Any additional           **
C  restrictions set forth in the license agreement also apply to this  **
C  software.                                                           **
C************************************************************************
C
      SUBROUTINE SETSEQ(NRES,LBRES,IPMOL,JSEQ,KSEQ,inew)
      COMMON/IOFILE/IN,IOUT,IOUTB
      DIMENSION LBRES(*),IPMOL(*),JSEQ(*),KSEQ(*)

C
C     ---- TRY THE OLD NOMENCLATURE STYLE ----
C
      IERR = 0
      DO 100 I = 1,NRES
        DO 120 J = 1,9
          IF(LBRES(I).EQ.JSEQ(J)) GO TO 140
  120   CONTINUE
        IERR = 2 
  140   CONTINUE
        IPMOL(I) = J
  100 CONTINUE
C
      IF(IERR.EQ.2) THEN
C
C       ---- TRY THE NEW NOMENCLATURE STYLE ----
C
        ierr = 1
        inew = 1
        DO 1000 I = 1,NRES
          DO 1200 J = 1,15
            IF(LBRES(I).EQ.KSEQ(J)) GO TO 1400
 1200     CONTINUE
          IERR = 2 
          write(6,'(a,a)') ' Unknown residue: ', LBRES(I)
 1400     CONTINUE
          IPMOL(I) = J
 1000   CONTINUE
C
      ENDIF
C
      IF(IERR.lt.2) RETURN
  900 CONTINUE
      WRITE(IOUT,*) 
     .   ' RESIDUE NAMES ARE NOT UNIFORMLY OLD OR 1994 NUCLEIC ACIDS'
      call mexit(iout, 1)
      END
