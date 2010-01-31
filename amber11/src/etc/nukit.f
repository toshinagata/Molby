C         PROGRAM TO INTERROGATE USER & GENERATE INPUT FILES
C         NUC.DAT, LIN.DAT FOR NUCGEN & LINK ROUTINES 
C
C                                BILL ROSS, TINOCO LAB, BERKELEY,1987
C
      program nukit
C
C
      CHARACTER*80 INLINE, OUTLINE
      CHARACTER CH, BLANK
      CHARACTER*3 RES
      CHARACTER FLAG
      BLANK = ' '
C
      write(6,*) 'Residue naming convention? (O = pre-94, N = 94) O/N:'
      READ(*,200)INLINE
      if (inline .eq. 'N') then
         inew = 1
      elseif (inline .eq. 'O') then
         inew = 0
      else
         write(6,*) 'unknown option: ', inline
         stop
      endif
      WRITE(6,420)
      READ(*,200)INLINE
C
      OPEN(UNIT=9,FILE='lin.in',STATUS='NEW')
      OPEN(UNIT=8,FILE='nuc.in',STATUS='NEW')
      WRITE(9,200)INLINE
      WRITE(9,*) ' '
      WRITE(9,500)' '
      OUTLINE(1:80) = 'DU' 
      WRITE(9,200)OUTLINE
      OUTLINE(1:80) = '     0    0    0    0    0  '
      WRITE(9,200)OUTLINE
      WRITE(6,700)
      WRITE(6,800)
c     WRITE(6,810)
c
c     -- do 2 strands
c
      DO 4 I = 1,2
          WRITE(9,400)I
          WRITE(8,400)I
    1     INLINE(1:80) = ' '
          WRITE(6,900)
          READ(*,200)INLINE
          OUTLINE(1:80) = '   '
          if (inew.eq.0) then
              OUTLINE(1:80) = 'HB '
              K = 1
          else
              k = 0
          endif
          WRITE(6,920)
          READ(*,500)FLAG
          IF (I .EQ. 1) THEN
              WRITE(9,500)FLAG
          ELSE
              WRITE(9,500)FLAG,1,3
          ENDIF
          WRITE(8,500)FLAG
          J = 1
C
C         primitive while() loop over residues in strand:
C
   77     IF (INLINE(J:J) .eq. BLANK) goto 88
              IF (K .GT. 10) THEN
                  WRITE(8,200)OUTLINE
                  OUTLINE(5:5) = '2'
                  WRITE(9,200)OUTLINE
                  OUTLINE(1:80) = ' '
                  K = 0
              END IF
              CH = INLINE(J:J)
              if (inew.eq.0) then
                  IF (CH .EQ. 'C') THEN
                       RES = 'CYT'
                  ELSEIF (CH .EQ. 'G') THEN
                       RES = 'GUA'
                  ELSEIF (CH .EQ. 'T') THEN
                       RES = 'THY'
                  ELSEIF (CH .EQ. 'U') THEN
                       RES = 'URA'
                  ELSEIF (CH .EQ. 'A') THEN        
                       RES = 'ADE'
                  ELSE 
                       WRITE(6,600)J
                       GOTO 1
                  ENDIF
              else
                  IF (CH .EQ. 'C') THEN
                       RES = 'C  '
                  ELSEIF (CH .EQ. 'G') THEN
                       RES = 'G  '
                  ELSEIF (CH .EQ. 'T') THEN
                       RES = 'T  '
                  ELSEIF (CH .EQ. 'U') THEN
                       RES = 'U  '
                  ELSEIF (CH .EQ. 'A') THEN
                       RES = 'A  '
                  ELSE
                       WRITE(6,600)J
                       GOTO 1
                  ENDIF
                  if (j.eq.1) res(2:2) = '5'
              endif
              J = J + 1
              OUTLINE((5*K+1):(5*K+3)) = RES
              K = K + 1
              IF (inew.eq.0 .and. INLINE(J:J) .NE. BLANK) THEN
                   OUTLINE((5*K+1):(5*K+3)) = 'POM'
                   K = K + 1
              ENDIF
          goto 77
   88     continue
          if (inew.eq.0) then
              OUTLINE((5*K+1):(5*K+2)) = 'HE'
          else
              outline((5*(k-1)+2):(5*(k-1)+2)) = '3'
          endif
          WRITE(8,200)OUTLINE
          OUTLINE(5:5) = '2'
          WRITE(9,200)OUTLINE
          OUTLINE(1:80) = ' '
          WRITE(8,200)OUTLINE
          WRITE(9,200)OUTLINE
    4 CONTINUE

      OUTLINE(1:4) = 'QUIT'
      WRITE(9,200)OUTLINE
      OUTLINE(1:4) = 'END '
      WRITE(8,200)OUTLINE
      WRITE(6,901)
      WRITE(6,902)
      WRITE(6,903)
      WRITE(6,904)
      WRITE(6,905)
      WRITE(6,906)
      WRITE(6,907)
      WRITE(6,908)
      WRITE(6,909)
      write(6,930)
      write(6,932)
      WRITE(6,910)
      INLINE(1:80) = ' '
      READ(5,200) INLINE
      WRITE(8,300) INLINE
      CLOSE(9)
      CLOSE(8)
C
  200 FORMAT(A80)
  300 FORMAT(t1, '$',a8)
  400 FORMAT('  NUC ',I3,I8,I4)
  420 FORMAT(1X,'JOB NAME?  ')
  500 FORMAT(A,I14,I5)
  600 FORMAT('***COULDNT MATCH BASE ',I2,' TRY AGAIN:'/ /)
  700 FORMAT('--------(from here on, USE CAPITALS)       ---------')
  800 FORMAT('---------e.g. CGCGATAT                     ---------')
  810 FORMAT('--------(Note: "8" =brom8gua, "5" =brom5cyt)--------')
  900 FORMAT(/ /1X,'ENTER SEQUENCE {5-prime to 3-prime}:   ')
  901 FORMAT('             CONFORMATIONS:')
  902 FORMAT('     ARNA   right handed A rna (arnott)')
  903 FORMAT('     APRNA  right handed A-prime rna (arnott)')
  904 FORMAT('     LBDNA  right handed B dna (langridge)')
  905 FORMAT('     ABDNA  right handed B dna (arnott)')
  906 FORMAT('     SBDNA  left handed  B dna (sasisekharan)')
  907 FORMAT('     ADNA   right handed A dna (arnott)')
  908 FORMAT('     NIXON  none of above - nucgen.pdb user-supplied')
  909 FORMAT(' ')
  930 format('  A-forms may need work to place H1-primes properly. ')
  932 format(' ')
  910 FORMAT(1X,'CONFORMATION?  ')
  920 FORMAT(1X,'DNA OR RNA? (D/R): ')
      END
