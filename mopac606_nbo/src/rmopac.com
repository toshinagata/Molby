$!   COMMAND FILE TO RUN A MOPAC JOB. THIS FILE SHOULD RESIDE IN THE
$!   MOPAC DIRECTORY. IT SHOULD BE ACCESSED FROM MOPAC.COM, BUT CAN
$!   STAND ALONE, IF NECESSARY.
$!
$!   THE CALL IS	$RMOPAC filename  directory  
$!
$	SHOW TIME
$!
$	SET VERIFY
$!
$!	Make assignments
$!
$DEL*ETE :== DELETE
$	IF P2 .NES. ""	THEN SET DEFAULT 'P2'
$	OPEN /ERR=NOEND1 DUMMY 'P1'.END  ! see if shutdown file is there
$       DELETE/NOCONFIRM 'P1'.END;*
$NOEND1:
$	ASSIGN 'P1'.DAT FOR005
$	ASSIGN 'P1'.OUT FOR006
$	ASSIGN 'P1'.RES FOR009
$	ASSIGN 'P1'.DEN FOR010
$       ASSIGN SYS$OUTPUT FOR011
$       ASSIGN SETUP.DAT SETUP
$	ASSIGN 'P1'.ARC FOR012
$	ASSIGN 'P1'.GPT FOR013
$	ASSIGN 'P1'.SYB FOR016
$	ASSIGN 'P1'.UMP FOR020
$       ASSIGN 'P1'.END SHUTDOWN
$!
$	ON ERROR	THEN GOTO AA
$	ON CONTROL_Y	THEN GOTO AA	! Cleanup if ^Y
$       TIM=F$TIME()
$       STARTTIME=F$CVTIME(TIM)
$	RUN/NODEBUG MOPACDIRECTORY:MOPAC
$       STOPTIME=F$CVTIME(F$TIME())
$!	SET NOVERIFY
$!
$	SHOW TIME
$!
$!      Get job information to mail to user at end of job.
$!      This will only happen if the elapsed time of the MOPAC
$!      run is greater than 2 hours.
$!
$!      Start the calculation of the time difference
$!
$ DAYS1=(F$INTEGER(F$EXTRACT(5,2,STARTTIME))*31) -
        +F$INTEGER(F$EXTRACT(8,2,STARTTIME))
$ DAYS2=(F$INTEGER(F$EXTRACT(5,2,STOPTIME))*31) -
        +F$INTEGER(F$EXTRACT(8,2,STOPTIME))
$ HRS1=F$INTEGER(F$EXTRACT(11,2,STARTTIME))
$ HRS2=F$INTEGER(F$EXTRACT(11,2,STOPTIME))
$ MIN1=F$INTEGER(F$EXTRACT(14,2,STARTTIME))
$ MIN2=F$INTEGER(F$EXTRACT(14,2,STOPTIME))
$ TDIFF=(((24*DAYS2)+HRS2)-((24*DAYS1)+HRS1))*60 +MIN2-MIN1
$!
$!      TDIFF now contains the number of minutess between the start and end
$!      of the job.  If this is greater than 5 then send mail to the user 
$!      (in case the user forgot what was being run)
$!
$  IF TDIFF .LE. 60 THEN GOTO AA
$!
$!   The job lasted for more than 60 minutes by the wall-clock
$!
$!   Open the .DAT file to get the job data to write into the message file
$!
$ OPEN/READ FILE1 'P2''P1'.DAT
$ READ FILE1 REC1
$ READ FILE1 REC2
$ READ FILE1 REC3
$ CLOSE FILE1
$!
$!  Write data out to file to be mailed to the user.
$!
$ STARTDAT=F$EXTRACT(0,17,TIM)
$ OPEN/WRITE FILE2 'P2''P1'.SUM
$ WRITE FILE2 "Job ''p1' run from directory ''p2' has finished."
$ WRITE FILE2 ""
$ WRITE FILE2 "Run started: ''startdat'"
$ WRITE FILE2 "Run time approximately: ''tdiff' minutes"
$ write file2 "Information extracted from ''p1'.DAT file follows:"
$ WRITE FILE2 ""
$ WRITE FILE2 rec1
$ WRITE FILE2 rec2
$ WRITE FILE2 rec3
$ CLOSE FILE2
$!
$!  Get name of user to whom to mail file
$!
$!   TO INSTALLER OF MOPAC:  ONE OF THE FOLLOWING SETS OF COMMANDS SHOULD
$!   WORK
$!                           SET 1
$ USENAM == F$EXTRACT(1,F$LENGTH(F$USER())-2,F$USER())
$!                           SET 2
$ DIRECTT := 'F$DIRECTORY()'
$ LDIRECTT = F$LENGTH(DIRECTT)
$ USENAM = F$EXTRACT(1,LDIRECTT-2,DIRECTT)
$ PPER = 'F$LOCATE(".",DIRECTT)'
$ IF PPER .LT. LDIRECTT THEN USENAM := 'F$EXTRACT(1,PPER-1,DIRECTT)'
$!                           END OF SET 2
$!
$!   Mail the file to user notifying end-of-job.
$!
$ MAIL 'P1'.SUM 'USENAM'
$ DELETE/NOCONFIRM 'P1'.SUM;
$!
$ AA:
$	SET NOCONTROL_Y			 ! Continue cleanup if ^Y
$       OPEN /ERR=NOEND DUMMY 'P1''PEND' ! see if shutdown file is there
$       DELETE/NOCONFIRM 'P1'.END;*
$NOEND: CLOSE DUMMY
$	SET CONTROL_Y
$	SET NOVERIFY
$!
$!	END
$!
