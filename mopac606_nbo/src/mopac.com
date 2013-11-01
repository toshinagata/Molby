$! COMMAND FILE TO SUBMIT A MOPAC JOB TO A BATCH QUEUE. 
$! THIS COMMAND CAN BE USED 'MANUALLY' BY ENTERING THE INSTRUCTION
$!    "$ @MOPAC <filename> <queue> <priority>" or
$!    "$ @MOPAC" and be prompted for the other arguments,
$! A RECOMMENDED APPROACH WOULD BE TO INSERT INTO THE LOGIN COMMAND THE LINE
$!    "$ MOPAC :== @MOPACDIRECTORY:MOPAC" This would allow the command
$!
$!		"$MOPAC <filename>  <queue>  <priority>"   to be used
$!
$!	Fetch parameters
$!
$	IF P1.NES."" THEN GOTO H
$ G:
$	INQUIRE P1 "What file? "
$ H:
$	LEN = 'F$LEN(P1)'
$	DOT = 'F$LOC(".",P1)
$	IF LEN.EQ.DOT THEN GOTO CHECK
$		LEXT = LEN - DOT
$		EXT := 'F$EXT(DOT,LEXT,P1)'
$		P1 := 'F$EXT(0,DOT,P1)'
$CHECK:
$!    check to see if file is there, if not send error message
$!
$! DEFINITIONS OF FILENAME EXTENSIONS.  FEEL FREE TO CHANGE THESE AS NECESSARY
$! IF YOU DO CHANGE THEM, THEN ALSO CHANGE RMOPAC.COM AND MOPACCOM.COM
$!
$       P05=".DAT"   !  DATA FILE
$! END OF DEFINITIONS
$	OPEN /ERR=NOFILE DUMMY 'P1''P05' 	! see if data file is there
$       READ DUMMY LINE
$       YLEN = 'F$LEN(LINE)'
$       YY = 'F$LOC("RESTART",LINE)'
$	CLOSE DUMMY				! Yes
$	GOTO OKAY
$NOFILE:
$	WRITE SYS$OUTPUT -
	" error opening ''P1'''P05'"
$	P1 := ""
$	GOTO G
$OKAY:
$ C:
$ K:
$!    SUBSTITUTE THE ACTUAL NAMES OF THE QUEUES FOR THE DUMMY NAMES
$!    QUEUE1 AND QUEUE2 IN THE FOLLOWING LINES. QUEUE3 IS INTENDED 
$!    FOR VERY LONG JOBS, QUEUE2 FOR JOBS OF LESS THAN 2 HOURS, AND
$!    QUEUE1 FOR JOBS OF LESS THAN HALF AN HOUR.
$                                   P4 = "QUEUE3"
$		IF P2.LT. 121  THEN P4 = "QUEUE2"
$		IF P2.LT.  31  THEN P4 = "QUEUE1"
$   P2=P4
$ M:
$!    IF YOU WANT THE DEFAULT PRIORITY TO BE DIFFERENT FROM 3 THEN
$!    ALTER THE FOLLOWING LINE.
$       IF P3.EQS."" THEN P3 := "3"
$!
$!	Submit
$!
$ TEMP1 := 'F$PARSE (P1+P05,,,"DEVICE")'
$ TEMP2 := 'F$PARSE (P1+P05,,,"DEVICE")''F$PARSE (P1+P05,,,"DIRECTORY")'
$ DIRECTT := 'F$DIRECTORY()'
$ LDIRECTT = F$LENGTH(DIRECTT)
$ LOG_FILE = DIRECTT
$!
$!  IF YOU WANT THE LOG FILE TO GO TO YOUR ROOT DIRECTORY, REMOVE THE
$!  EXCLAMATION POINTS ON THE FOLLOWING TWO LINES.
$! PPER = 'F$LOCATE(".",DIRECTT)'
$! IF PPER .LT. LDIRECTT THEN LOG_FILE := 'F$EXTRACT(0,PPER,DIRECTT)']
$!  IF YOU WANT THE LOG FILE TO GO TO YOUR ROOT DIRECTORY, REMOVE THE
$!  EXCLAMATION POINTS ON THE PREVIOUS TWO LINES.
$!
$SUBMIT/notify /NAME='P1' MOPACDIRECTORY:RMOPAC /QUEUE='P2' -
/PARA=("''P1'","''TEMP2'")/PRIO='P3' -
/LOG_FILE='TEMP1''LOG_FILE''P1'.LOG/NOPRINT
$WRITE SYS$OUTPUT "  "
$WRITE SYS$OUTPUT "  "
$WRITE SYS$OUTPUT "         Notice of Public Domain nature of this Program  "  
$WRITE SYS$OUTPUT "  "
$WRITE SYS$OUTPUT "      'This computer program is a work of the United States   "
$WRITE SYS$OUTPUT "       Government and as such is not subject to protection by   "
$WRITE SYS$OUTPUT "       copyright (17 U.S.C. # 105.)  Any person who fraudulently   "
$WRITE SYS$OUTPUT "       places a copyright notice or does any other act contrary   "
$WRITE SYS$OUTPUT "       to the provisions of 17 U.S. Code 506(c) shall be subject   "
$WRITE SYS$OUTPUT "       to the penalties provided therein.  This notice shall not   "
$WRITE SYS$OUTPUT "       be altered or removed from this software and is to be on   "
$WRITE SYS$OUTPUT "       all reproductions.'  "
$WRITE SYS$OUTPUT "  "
$EXIT
