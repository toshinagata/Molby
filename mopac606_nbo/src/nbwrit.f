C*****************************************************************************
      SUBROUTINE NBWRIT(IX,NX,IDAR)
C*****************************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (LENGTH = 256)
      PARAMETER (NBDAR = 100)
C
      COMMON/NBODAF/INBO,NAV,IONBO(NBDAR)
      COMMON/NBONAV/IXDNBO(LENGTH),NBNAV,ISINGL
C
      DIMENSION IX(1),IXSNBO(LENGTH/2)
C
      EQUIVALENCE (IXSNBO(1),IXDNBO(1))
C
      MAXIX = LENGTH * ISINGL / 2
      LDAR  = NX * ISINGL
      IF(IONBO(IDAR).NE.0) GO TO 100
C
C  If this is the first write to the NBO DAF:
C
      IONBO(IDAR) = NAV
      NBNAV = NAV
C
      MAX = 0
   10 MIN = MAX + 1
      MAX = MAX + MAXIX
      IF(MAX.GT.LDAR) MAX = LDAR
      DO 20 I = MIN,MAX
   20 IXDNBO(I-MIN+1) = IX(I)
      IF(ISINGL.EQ.1) WRITE(INBO,REC=NBNAV) IXSNBO
      IF(ISINGL.EQ.2) WRITE(INBO,REC=NBNAV) IXDNBO
      NBNAV = NBNAV + 1
      IF(MAX.LT.LDAR) GO TO 10
      NAV = NBNAV
      RETURN
C
C  Or if this is a rewrite:
C
  100 CONTINUE
      NBNAV = IONBO(IDAR)
      MAX = 0
  110 MIN = MAX + 1
      MAX = MAX + MAXIX
      IF(MAX.GT.LDAR) MAX = LDAR
      DO 120 I = MIN,MAX
  120 IXDNBO(I-MIN+1) = IX(I)
      IF(ISINGL.EQ.1) WRITE(INBO,REC=NBNAV) IXSNBO
      IF(ISINGL.EQ.2) WRITE(INBO,REC=NBNAV) IXDNBO
      NBNAV = NBNAV + 1
      IF(MAX.LT.LDAR) GO TO 110
      RETURN
      END
