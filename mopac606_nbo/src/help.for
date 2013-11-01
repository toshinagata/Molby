      PROGRAM MOHELP
************************************************************************
*                                                                      *
*                        PROGRAM  MOHELP                               *
*                                                                      *
*    WRITTEN BY CPTN DONN STORCH, USAF, AND JAMES J. P. STEWART        *
*                                                                      *
* USAGE: A HELP FILE, WRITTEN IN VAX-FORMAT IS USED AS DATA.           *
*        HELP IS FOR ON-LINE USE; DATA ARE SUPPLIED IN THE SAME FORMAT *
*        AS FOR THE VAX HELP COMMAND, FOR EXAMPLE THE DATA LINE        *
*                                                                      *
*        KEYS                                                          *
*                                                                      *
*        WOULD GIVE A BRIEF DEFINITION OF THE KEY-WORDS, AND A FULL    *
*        LIST OF ALL OF THE KEY-WORDS                                  *
*                                                                      *
*        WRITTEN IN 1984                                               *
************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*80 COMAND
C
C  IF HELP FILES FOR MORE THAN ONE PROGRAM EXIST, UNCOMMENT THE
C  FOLLOWING LINES
C
C#  10  OPEN ( UNIT=9, FILE= 'MOH', STATUS='OLD')
C#      REWIND 9
C
C  FILE MOH SHOULD CONTAIN THE NAMES OF THE PROGRAMS FOR WHICH HELP
C  CAN BE PROVIDED.  FILES SHOULD BE THREE LETTER ABBREVIATIONS,
C  EXAMPLE, FOR MOPAC 'MOP'
C Help is available on the following topics
C
C      Topic         Acceptable Abbreviation
C
C      DENSITY             DEN
C      DRAW                DRA
C      MOPAC               MOP
C      HELP                HEL
C
C#  20  READ(9,'(A)',END=30,ERR=30)COMAND
C#      WRITE(6,'(A)')COMAND
C#      GOTO 20
C#  30  WRITE(6,'('' WHICH PROGRAM DO YOU WANT HELP WITH'')')
C#      READ(5,'(A)')COMAND
C#      CALL LCLEAN(COMAND,COMAND)
C#      COMAND=COMAND(1:3)
C#      IF(COMAND(1:1).EQ.' ') STOP
C
C  IF THE HELP FILE DOES NOT RESIDE IN YOUR HOME DIRECTORY
C  CHANGE THE NEXT LINE
      COMAND='MOHELP'
      OPEN ( UNIT=9, FILE= COMAND, STATUS='OLD',ERR=10)
      REWIND 9
      COMAND='MOP'
      CALL HELP(COMAND)
      STOP
   10 WRITE(6,'(A)')' MOPAC.HLP FILE MISSING'
      END
      SUBROUTINE HELP(COMAND)
************************************************************************
*                                                                      *
*                        SUBROUTINE HELP                               *
*                                                                      *
*    WRITTEN BY CPTN DONN STORCH, USAF, AND JAMES J. P. STEWART        *
*                                                                      *
* USAGE: A HELP FILE, WRITTEN IN VAX-FORMAT IS USED AS DATA.           *
*        HELP IS FOR ON-LINE USE; DATA ARE SUPPLIED IN THE SAME FORMAT *
*        AS FOR THE VAX HELP COMMAND, FOR EXAMPLE THE DATA LINE        *
*                                                                      *
*        KEYS                                                          *
*                                                                      *
*        WOULD GIVE A BRIEF DEFINITION OF THE KEY-WORDS, AND A FULL    *
*        LIST OF ALL OF THE KEY-WORDS                                  *
*                                                                      *
*        WRITTEN IN 1984                                               *
************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /DEBCOM/ DEBUG, DEBUGL, DEBUGN, DEBUGO, DEBUGP
      CHARACTER*80 LINE, COMAND, TLINE, SUBCOM
      CHARACTER*20 CWORDS(10), SPACE
      LOGICAL DEBUGL, DEBUGN, DEBUGO, BEGIN
      LOGICAL ERROR, FOUND, FIRST, MULTI, DEBUG, MORE, FEXIT, DEBUGP
      DIMENSION IRECNO( 10), LENCW( 10)
      DATA SPACE/'                    '/, BEGIN /.TRUE./
      IF( BEGIN )THEN
         IREAD=5
         I=INDEX(COMAND,'READ')
         IF(I.NE.0)IREAD=READA(COMAND,I)+0.5
         BEGIN=.FALSE.
      ENDIF
      FIRST = .TRUE.
      LEVEL = 0
      IRCD = 0
   10 LCMD = INDEX ( COMAND, ' ')-1
      FOUND = .FALSE.
      DO 20 LTOTL=80,1,-1
         IF ( COMAND(LTOTL:LTOTL) .NE. ' ') GOTO 30
   20 CONTINUE
   30 CONTINUE
* if COMAND ends in '-' then exit after listing only subcomands
      FEXIT = .FALSE.
      IF ( COMAND( LTOTL:LTOTL) .EQ. '-') THEN
         FEXIT = .TRUE.
         COMAND( LTOTL: LTOTL) = ' '
         LTOTL = LTOTL - 1
      ENDIF
      SUBCOM = COMAND(1:LCMD)
      LEVEL = LEVEL + 1
      IRECNO( LEVEL) = 0
      IF ( DEBUG) WRITE (*,*) 'SUBCOM='//SUBCOM(1:LCMD)//
     1      '; LCMD=',LCMD,', LTOTL=',LTOTL
      IF ( LCMD .LT. LTOTL) THEN
         COMAND = COMAND(LCMD+2:LTOTL)
         MULTI = .TRUE.
      ELSE
         COMAND = '          '
         MULTI = .FALSE.
      ENDIF
* here we are searching for the topic requested
   40 READ ( 9, '( A )', END=270 ) LINE
      IRCD = IRCD + 1
*
      IF ( LINE(1:1) .GE. '0' .AND. LINE(1:1) .LE. '9' 
     + .AND. LINE(2:2).EQ.' ') THEN
         ITEMP = READA ( LINE, 1, ERROR)
         IF ( ITEMP .LT. LEVEL ) GOTO 40
         DO 50 I=2,80
            IF ( LINE(I:I) .NE. ' ') GOTO 60
   50    CONTINUE
   60    CONTINUE
         IF ( LINE(I:I+LCMD-1) .EQ. SUBCOM(:LCMD) .AND.
     1        ITEMP .EQ. LEVEL ) THEN
            FOUND = .TRUE.
            CWORDS(LEVEL) = LINE(I:INDEX(LINE(I:),' ')+I-1)
            IRECNO( LEVEL) = IRCD
            DO 70 KTEMP = 80, I, -1
               IF ( LINE( KTEMP: KTEMP) .NE. ' ') GOTO 80
   70       CONTINUE
   80       LENCW( LEVEL) = KTEMP - I + 1
            CWORDS( LEVEL) = LINE( I: KTEMP)
            IF ( DEBUG ) WRITE (*,*) 'LABEL FOUND on RECORD',
     1                     IRECNO( LEVEL)
         ELSEIF ( FOUND ) THEN
            GOTO 150
         ELSE
            IF ( FOUND ) GOTO 150
            GOTO 40
         ENDIF
      ENDIF
      IF ( .NOT. FOUND ) GOTO 40
      IF ( MULTI ) GOTO 10
* topic found, print all info on it until sub topic found
   90 CONTINUE
      WRITE ( *, *)
      DO 100 I=1,LEVEL
         WRITE ( *, *) SPACE(:I*2)//CWORDS(I)(:LENCW( I))
  100 CONTINUE
      WRITE ( *, *)
  110 READ ( 9, '( A )', END=230) LINE
      IRCD = IRCD + 1
      IF ( LINE(1:1) .LT. '0' .OR. LINE(1:1) .GT. '9'
     +.OR. LINE(2:2).NE.' ' ) THEN
         DO 120 I=80,1,-1
            IF ( LINE(I:I) .NE. ' ') GOTO 130
  120    CONTINUE
         GOTO 110
  130    WRITE ( *, *) LINE(:I)
         GOTO 110
      ENDIF
* now search and keep all sub-topics
  140 CONTINUE
      MORE = .FALSE.
      ITEMP = READA ( LINE, 1, ERROR)
      IF ( ITEMP .LE. LEVEL ) GOTO 230
  150 WRITE ( *, '( 1X/,'' Other commands available:'',/ )' )
      MORE = .TRUE.
      INDX = 1
  160 DO 170 I=2,80
         IF ( LINE(I:I) .NE. ' ') GOTO 180
  170 CONTINUE
  180 CONTINUE
      DO 190 ITEMP=80,2,-1
         IF ( LINE(ITEMP:ITEMP) .NE. ' ') GOTO 200
  190 CONTINUE
  200 CONTINUE
      IF ( DEBUG ) WRITE ( *, *) 'GATHERING KEYS:'//LINE(I:ITEMP)
      TLINE = TLINE(:INDX)//LINE(I:ITEMP)
      INDX = INDX + ITEMP - I + 1
      IF ( INDX .GT. 60 ) THEN
         WRITE ( *, '( 1X, A )' ) TLINE(:INDX)
         INDX = 1
      ELSE
         ITEMP = 15 - MOD( INDX, 15)
         TLINE = TLINE(:INDX)//SPACE(:ITEMP)
         INDX = INDX + ITEMP
      ENDIF
  210 READ ( 9, '( A )', END=220) LINE
      IRCD = IRCD + 1
      IF ( LINE(1:1) .LT. '0' .OR. LINE(1:1) .GT. '9'
     +.OR.LINE(2:2).NE.' ') GOTO 210
      ITEMP = READA ( LINE, 1, ERROR)
      IF ( ITEMP .LE. LEVEL) GOTO 220
      IF ( ITEMP-1 .GT. LEVEL ) GOTO 210
      GOTO 160
*
  220 WRITE ( *, '( 1X, A )' ) TLINE(:INDX)
      IF ( FEXIT) GOTO 310
  230 CONTINUE
      WRITE ( *, *)
      IF ( .NOT. MORE ) THEN
         LEVEL = LEVEL - 1
         IF ( LEVEL .LT. 1 ) GOTO 310
      ENDIF
      INDX = 1
      DO 240 I= 1, LEVEL
         TLINE = TLINE(1:INDX)//CWORDS( I)(:LENCW( I))//': '
         INDX = INDX + LENCW( I) + 2
  240 CONTINUE
      CALL UPROMP( 'HELP: '// TLINE(:INDX) )
      READ ( IREAD, '( A )', END=250 ) COMAND
      CALL LCLEAN( COMAND, COMAND)
  250 IF ( COMAND(1:1) .EQ. ' ') THEN
         IF ( MORE ) LEVEL = LEVEL - 1
         IF ( LEVEL .LT. 1 ) GOTO 310
         GOTO 230
      ELSEIF ( COMAND( 1: 1) .EQ. '?') THEN
         GOTO 280
      ELSE
         REWIND ( UNIT= 9)
         DO 260 I= 1, IRECNO( LEVEL)
            READ ( 9, '( A )' ) LINE
  260    CONTINUE
         IRCD = IRECNO( LEVEL)
         GOTO 10
      ENDIF
*
  270 WRITE ( *, *) 'No information available on '//SUBCOM(:LCMD)
      LEVEL = LEVEL - 1
      IF ( LEVEL .LT. 1) GOTO 310
  280 REWIND ( UNIT= 9)
      DO 290 I= 1, IRECNO( LEVEL)
         READ ( 9, '( A )' ) LINE
  290 CONTINUE
      IRCD = IRECNO( LEVEL)
      IF ( DEBUG ) WRITE (*,*) 'LEVEL=',LEVEL
  300 READ ( 9, '( A )', END=230 ) LINE
      IRCD = IRCD + 1
      IF ( LINE( 1: 1) .GE. '0' .AND. LINE( 1: 1) .LE. '9'
     +.AND. LINE(2:2).EQ.' ' ) THEN
         ITEMP = READA ( LINE, 1, ERROR)
         IF ( ITEMP .LT. LEVEL ) GOTO 230
         IF ( ITEMP .GT. LEVEL ) GOTO 140
         IF ( ITEMP .EQ. LEVEL ) GOTO 300
      ENDIF
      GOTO 300
*
  310 CONTINUE
      CLOSE ( UNIT= 9)
      RETURN
      END
      SUBROUTINE UPROMP( TEXT )
      CHARACTER*(*) TEXT
*************************************************
*
*   USE OF A DOLLAR SIGN IS NON-STANDARD FORTRAN
*   BUT VAX USERS MIGHT WANT TO PUT ONE AFTER
*   THE 'A' FORMAT SPECIFIER
*
************************************************
      IF ( TEXT(1:1) .EQ. '+' ) THEN
         WRITE ( *, '(A)' ) TEXT
      ELSE
         WRITE ( *, '( 1X, A)') TEXT
      ENDIF
*
      RETURN
      END
      SUBROUTINE LCLEAN (LINE1, LINE2)
      LOGICAL DEBUGL, DEBUGN, DEBUGO
      LOGICAL DEBUG,DEBUGP
      COMMON /DEBCOM/ DEBUG, DEBUGL, DEBUGN, DEBUGO, DEBUGP
      CHARACTER*1 TAB, ASCII
      CHARACTER*(*) LINE1, LINE2
      COMMON /ASCIIC/ ASCII( 0: 255)
      EQUIVALENCE ( TAB, ASCII(9))
*
      I = 0
      ILOWA = ICHAR('a')
      ILOWZ = ICHAR('z')
      ICAPA = ICHAR('A')
************************************************************************
      LINE2=LINE1
      DO 10 I=1,80
         ILINE=ICHAR(LINE2(I:I))
         IF(ILINE.GE.ILOWA.AND.ILINE.LE.ILOWZ) THEN
            LINE2(I:I)=CHAR(ILINE+ICAPA-ILOWA)
         ENDIF
   10 CONTINUE
************************************************************************
   20 IF ( LINE2(1:1) .EQ. ' ') THEN
         LINE2 = LINE2( 2:)
         I = I + 1
         IF ( I .GT. 10) RETURN
         GOTO 20
      ENDIF
   30 I = INDEX(LINE2,',')
      IF (I .GT. 0) THEN
         LINE2(I:I)=' '
         GO TO 30
      ENDIF
   40 I = INDEX(LINE2,';')
      IF (I .GT. 0) THEN
         LINE2(I:I)=' '
         GO TO 40
      ENDIF
   50 I = INDEX(LINE2,TAB)
      IF (I .GT. 0) THEN
         LINE2(I:I)=' '
         GO TO 50
      ENDIF
   60 I = INDEX(LINE2,'/')
      IF (I .GT. 0) THEN
         LINE2(I:I)=' '
         GO TO 60
      ENDIF
      RETURN
      END
      DOUBLE PRECISION FUNCTION READA(STRING,ISTART)
C     FORTRAN FUNCTION TO EXTRACT NUMBER FROM STRING
C
      CHARACTER STRING*(*)
      DOUBLE PRECISION DIGIT
      LOGICAL DEFALT,EXPNNT
C
C     DEFINE ASCII VALUES OF NUMERIC FIELD CHARACTERS
      I0=ICHAR('0')
      I9=ICHAR('9')
      IDOT=ICHAR('.')
      INEG=ICHAR('-')
      IPOS=ICHAR('+')
      ICAPD=ICHAR('D')
      ICAPE=ICHAR('E')
      ISMLD=ICHAR('d')
      ISMLE=ICHAR('e')
C
      L=LEN(STRING)
C
C     FIND THE START OF THE NUMERIC FIELD
      DO 10 I=ISTART,L
         IADD=0
         N=ICHAR(STRING(I:I))
C
C       SIGNAL START OF NUMERIC FIELD IF DIGIT FOUND
         IF(N.GE.I0.AND.N.LE.I9)GOTO 20
C
C       ACCOUNT FOR CONSECUTIVE SIGNS [- AND(OR) +]
         IF(N.EQ.INEG.OR.N.EQ.IPOS)THEN
            IADD=IADD+1
            IF(I+IADD.GT.L)GOTO 50
            N=ICHAR(STRING(I+IADD:I+IADD))
            IF(N.GE.I0.AND.N.LE.I9)GOTO 20
         ENDIF
C
C       ACCOUNT FOR CONSECUTIVE DECIMAL POINTS (.)
         IF(N.EQ.IDOT)THEN
            IADD=IADD+1
            IF(I+IADD.GT.L)GOTO 50
            N=ICHAR(STRING(I+IADD:I+IADD))
            IF(N.GE.I0.AND.N.LE.I9)GOTO 20
         ENDIF
   10 CONTINUE
      GOTO 50
C
C     FIND THE END OF THE NUMERIC FIELD
   20 EXPNNT=.FALSE.
      DO 30 J=I+1,L
         IADD=0
         N=ICHAR(STRING(J:J))
C
C       CONTINUE SEARCH FOR END IF DIGIT FOUND
         IF(N.GE.I0.AND.N.LE.I9)GOTO 30
C
C       CONTINUE SEARCH FOR END IF SIGN FOUND AND EXPNNT TRUE
         IF(N.EQ.INEG.OR.N.EQ.IPOS)THEN
            IF(.NOT.EXPNNT)GOTO 40
            IADD=IADD+1
            IF(J+IADD.GT.L)GOTO 40
            N=ICHAR(STRING(J+IADD:J+IADD))
            IF(N.GE.I0.AND.N.LE.I9)GOTO 30
         ENDIF
         IF(N.EQ.IDOT)THEN
            IADD=IADD+1
            IF(J+IADD.GT.L)GOTO 40
            N=ICHAR(STRING(J+IADD:J+IADD))
            IF(N.GE.I0.AND.N.LE.I9)GOTO 30
            IF(N.EQ.ICAPE.OR.N.EQ.ISMLE.OR.N.EQ.ICAPD.OR.N.EQ.ISMLD)
     1    GOTO 30
         ENDIF
         IF(N.EQ.ICAPE.OR.N.EQ.ISMLE.OR.N.EQ.ICAPD.OR.N.EQ.ISMLD)THEN
            IF(EXPNNT)GOTO 40
            EXPNNT=.TRUE.
            GOTO 30
         ENDIF
         GOTO 40
   30 CONTINUE
      J=L+1
   40 N=ICHAR(STRING(J-1:J-1))
      IF(N.EQ.ICAPE.OR.N.EQ.ISMLE.OR.N.EQ.ICAPD.OR.N.EQ.ISMLD)J=J-1
C
C     FOUND THE END OF THE NUMERIC FIELD (IT RUNS 'I' THRU 'J-1')
      N=0
      N=N+INDEX(STRING(I:J-1),'e')
      N=N+INDEX(STRING(I:J-1),'E')
      N=N+INDEX(STRING(I:J-1),'d')
      N=N+INDEX(STRING(I:J-1),'D')
      IF(N.EQ.0)THEN
         READA=DIGIT(STRING(I:J-1),1)
      ELSE
         READA=DIGIT(STRING(:I+N-2),I)*1.D1**DIGIT(STRING(:J-1),I+N)
      ENDIF
      DEFALT=.FALSE.
      RETURN
C
C     DEFAULT VALUE RETURNED BECAUSE NO NUMERIC FIELD FOUND
   50 READA=0.D0
      DEFALT=.TRUE.
      RETURN
      END
C     ******************************************************************
      DOUBLE PRECISION FUNCTION DIGIT(STRING,ISTART)
C     FORTRAN FUNCTION TO CONVERT NUMERIC FIELD TO DOUBLE PRECISION
C     NUMBER.  THE STRING IS ASSUMED TO BE CLEAN (NO INVALID DIGIT
C     OR CHARACTER COMBINATIONS FROM ISTART TO THE FIRST NONSPACE,
C     NONDIGIT, NONSIGN, AND NONDECIMAL POINT CHARACTER).
C
      CHARACTER STRING*(*)
      DOUBLE PRECISION C1,C2,DECIML
      LOGICAL SIGN
C
C     DEFINE ASCII VALUES OF NUMERIC FIELD CHARACTERS
      I0=ICHAR('0')
      I9=ICHAR('9')
      INEG=ICHAR('-')
      IPOS=ICHAR('+')
      IDOT=ICHAR('.')
      ISPC=ICHAR(' ')
C
      C1=0.D0
      C2=0.D0
      SIGN=.TRUE.
      L=LEN(STRING)
C
C     DETERMINE THE CONTRIBUTION TO THE NUMBER GREATER THAN ONE
      IDIG=0
      DO 10 I=ISTART,L
         N=ICHAR(STRING(I:I))
         IF(N.GE.I0.AND.N.LE.I9)THEN
            IDIG=IDIG+1
            C1=C1*1.D1+N-I0
         ELSEIF(N.EQ.INEG.OR.N.EQ.IPOS.OR.N.EQ.ISPC)THEN
            IF(N.EQ.INEG)SIGN=.FALSE.
         ELSEIF(N.EQ.IDOT)THEN
            GOTO 20
         ELSE
            GOTO 40
         ENDIF
   10 CONTINUE
C
C     DETERMINE THE CONTRIBUTION TO THE NUMBER LESS THAN THAN ONE
   20 DECIML=1.D0
      DO 30 J=I+1,L
         N=ICHAR(STRING(J:J))
         IF(N.GE.I0.AND.N.LE.I9)THEN
            DECIML=DECIML/1.D1
            C2=C2+(N-I0)*DECIML
         ELSEIF(N.NE.ISPC)THEN
            GOTO 40
         ENDIF
   30 CONTINUE
C
C     PUT THE PIECES TOGETHER
   40 DIGIT=C1+C2
      IF(.NOT.SIGN)DIGIT=-DIGIT
      RETURN
      END
