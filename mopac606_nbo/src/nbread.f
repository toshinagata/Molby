C*****************************************************************************
      SUBROUTINE NBREAD(IX,NX,IDAR)
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
      NBNAV = IONBO(IDAR)
      MAXIX = LENGTH * ISINGL / 2
      LDAR  = NX * ISINGL
C
      MAX = 0
   10 MIN = MAX + 1
      MAX = MAX + MAXIX
      IF(MAX.GT.LDAR) MAX = LDAR
      IF(ISINGL.EQ.1) READ(INBO,REC=NBNAV) IXSNBO
      IF(ISINGL.EQ.2) READ(INBO,REC=NBNAV) IXDNBO
      DO 20 I = MIN,MAX
   20 IX(I) = IXDNBO(I-MIN+1)
      NBNAV = NBNAV + 1
      IF(MAX.LT.LDAR) GO TO 10
      RETURN
      END
