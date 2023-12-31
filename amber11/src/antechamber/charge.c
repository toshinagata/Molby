# define CONVERG 0.00001
# define GASMAXITER 500
# define DAMPFACTOR 0.5
GASTEIGER gas[MAXGAS];

char *amberhome;

void wpdb_optimized(char *filename, int atomnum, ATOM atom[] ,int flag) {
ATOM *atom_tmp;
int i;
/* when read in a mopout file, always output a pdb file that has the optimized coordinates */
                atom_tmp = (ATOM *) malloc(sizeof(ATOM) * atomnum);
                if (atom_tmp == NULL) {
                        fprintf(stderr, "memory allocation error for *atom_tmp\n");
                        exit(1);
                }
                for(i=0;i<atomnum;i++) {
                        atom_tmp[i] = atom[i];
			atom_tmp[i].connum = 0;
		}
		if(flag == 0) {
                	rmopout_coord(filename, atom_tmp);
                	wpdb("mopac.pdb", atomnum, atom_tmp);
		}
		if(flag == 1) {
                	rdivout_coord(filename, atom_tmp);
                	wpdb("divcon.pdb", atomnum, atom_tmp);
		}
                free(atom_tmp);
}

void rsqmcharge(char *filename, int atomnum, ATOM atom[], MOLINFO *minfo)
{

	int index;
	int number = 0;
	int tmpint;
	FILE *fpout;
	char tmpchar[MAXCHAR];
	char tmpchar1[MAXCHAR];
	char tmpchar2[MAXCHAR];
	char tmpchar3[MAXCHAR];
	char tmpchar4[MAXCHAR];
	char tmpchar5[MAXCHAR];
	char tmpchar6[MAXCHAR];
	char line[MAXCHAR];
	double tmpfloat1;

	if ((fpout = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "Cannot open sqm output file: %s in rmopacharge(), exit\n", filename);
		exit(1);
	}
	index = 0;
	for (;;) {
		if (fgets(line, MAXCHAR, fpout) == NULL)
			break;

		sscanf(line, "%s%s%s%s%s%s", tmpchar1, tmpchar2, tmpchar3,
			   tmpchar4, tmpchar5, tmpchar6);
		if (strcmp("Mulliken", tmpchar3) == 0 ){
			index = 1;
			continue;
		}
		if (index == 1) {
			sscanf(line, "%d%s%lf", &tmpint, tmpchar, &tmpfloat1);
		 	atom[number].charge = tmpfloat1;
			number++;
		}
	}
	fclose(fpout);
	if (number == 0) {
		fprintf(stderr, "Error: unable to find sqm charges in %s\n", filename);
		fprintf(stderr, "       Examine that file for evidence of errors\n" );
		exit(1);
	}
	if ((*minfo).usercharge < -9990)
		(*minfo).icharge = intcharge(atomnum, atom);
}

void rmopcharge(char *filename, int atomnum, ATOM atom[], MOLINFO *minfo)
{

	int index;
	int number = 0;
	int tmpint;
	FILE *fpout;
	char tmpchar[MAXCHAR];
	char tmpchar1[MAXCHAR];
	char tmpchar2[MAXCHAR];
	char tmpchar3[MAXCHAR];
	char tmpchar4[MAXCHAR];
	char tmpchar5[MAXCHAR];
	char tmpchar6[MAXCHAR];
	char line[MAXCHAR];
	double tmpfloat1;
	double tmpfloat2;

	if ((fpout = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "Cannot open mopac output file: %s in rmopacharge(), exit\n", filename);
		exit(1);
	}
	index = 0;
	for (;;) {
		if (fgets(line, MAXCHAR, fpout) == NULL)
			break;

		sscanf(line, "%s%s%s%s%s%s", tmpchar1, tmpchar2, tmpchar3,
			   tmpchar4, tmpchar5, tmpchar6);
		if (strcmp("NET", tmpchar1) == 0 && strcmp("ATOMIC", tmpchar2) == 0
			&& strcmp("CHARGES", tmpchar3) == 0
			&& strcmp("AND", tmpchar4) == 0
			&& strcmp("DIPOLE", tmpchar5) == 0
			&& strcmp("CONTRIBUTIONS", tmpchar6) == 0) {
			index = 1;
			continue;
		}
		if (index == 1 && strcmp("ATOM", tmpchar1) == 0
			&& strcmp("NO.", tmpchar2) == 0) {
			index = 2;
			continue;
		}
		if (index == 2 && strcmp("DIPOLE", tmpchar1) == 0
			&& strcmp("X", tmpchar2) == 0 && strcmp("Y", tmpchar3) == 0
			&& strcmp("Z", tmpchar4) == 0
			&& strcmp("TOTAL", tmpchar5) == 0) {
			index = 3;
			break;
		}
		if (index == 2) {
			sscanf(line, "%d%s%lf%lf", &tmpint, tmpchar, &tmpfloat1,
			   	&tmpfloat2);
		 	atom[number].charge = tmpfloat1;
			number++;
		}
	}
	fclose(fpout);
	if (number == 0) {
		fprintf(stderr, "Error: unable to find mopac charges in %s\n", filename);
		fprintf(stderr, "       Examine that file for evidence of errors\n" );
		exit(1);
	}
	if ((*minfo).usercharge < -9990)
		(*minfo).icharge = intcharge(atomnum, atom);
}


void rdivcharge(char *filename, int atomnum, ATOM atom[], MOLINFO *minfo, int flag)
{

	/* now modified for divcon.out output */

	int index;
	int number = 0;
	int tmpint;
	FILE *fpout;
	char tmpchar[MAXCHAR];
	char tmpchar1[MAXCHAR];
	char tmpchar2[MAXCHAR];
	char tmpchar3[MAXCHAR];
	char tmpchar4[MAXCHAR];
	char tmpchar5[MAXCHAR];
	char tmpchar6[MAXCHAR];
	char line[MAXCHAR];
	double tmpfloat1;
	double tmpfloat2;
	double tmpfloat3;

	if ((fpout = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "Cannot open divcon output file: %s in rdivcharge() , exit\n", filename);
		exit(1);
	}
	index = 0;
	for (;;) {
		if (fgets(line, MAXCHAR, fpout) == NULL)
			break;

		sscanf(line, "%s%s%s%s%s%s", tmpchar1, tmpchar2, tmpchar3,
			   tmpchar4, tmpchar5, tmpchar6);
		/* parse DIVCON output  */
		if (strcmp("NO.", tmpchar1) == 0 && strcmp("SYMBOL", tmpchar2) == 0
			&& strcmp("MULLIKEN", tmpchar3) == 0
			&& strcmp("CM1", tmpchar4) == 0
			&& strcmp("CM2", tmpchar5) == 0 ){
			index = 1;
			continue;
		}
		if (index == 1 && strcmp("CHARGE", tmpchar1) == 0
			&& strcmp("CHARGE", tmpchar2) == 0) {
			index = 2;
			fgets(line, MAXCHAR, fpout);
			continue;
		}
		if (index == 2 && strcmp("\n", line) == 0 ) {
			index = 3;
			break;
		}
		if (index == 2) {
			sscanf(line, "%d%s%lf%lf%lf", &tmpint, tmpchar, &tmpfloat1,
			   	&tmpfloat2, &tmpfloat3);
			if(flag == 1) 
				atom[number].charge = tmpfloat1;
			if(flag == 2) 
				atom[number].charge = tmpfloat2;
			if(flag == 3) 
				atom[number].charge = tmpfloat3;
			number++;
		}
	}
	fclose(fpout);
	if (number == 0) {
		fprintf(stderr, "Error: unable to find divcon charges in %s\n", filename);
		exit(1);
	}
	if ((*minfo).usercharge < -9990)
		(*minfo).icharge = intcharge(atomnum, atom);
}

void rgaucharge(char *filename, char *chargemethod, int atomnum,
				ATOM atom[], MOLINFO *minfo)
{
	int chargeindex;
	int num;
	int Found_Stationary = 0;
	int index = 0;
	char line[MAXCHAR];
	FILE *fpin;

	if ((fpin = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "Cannot open the gaussian output file: %s in rgaucharge(), exit\n", filename);
		return;
	}
	if (strcmp(chargemethod, "mul") == 0)
		chargeindex = 1;
	else
		chargeindex = 2;
	num = 0;
	for (;;) {
		if (fgets(line, MAXCHAR, fpin) == NULL)
			break;
		if (strncmp("-- Stationary point found", &line[4], 25) == 0) {
			Found_Stationary = 1; 
			continue; 
		}
		if (chargeindex == 1 && Found_Stationary == 1 && index == 0) {
			if(strncmp("Total atomic charges:", &line[1], 21) == 0) index = 1;
			if(strncmp("Mulliken atomic charges:", &line[1], 24) == 0) index = 1;
			continue;
		}
		if (chargeindex == 2 && Found_Stationary == 1 && index == 0 && 
			strncmp("Fitting point charges", &line[1], 21) == 0)
			index = 1;
		if (index ==1 && strncmp("Atomic charges with hydrogens summed into heavy atoms", 
			&line[1], 52) == 0)
			index = 2;
		if (Found_Stationary == 1 && index == 1) {
			if(line[11]=='.' ) {
				sscanf(&line[8], "%lf", &atom[num].charge);
				num = num + 1;
			}
			if(line[14]=='.') {
				sscanf(&line[11], "%lf", &atom[num].charge);
				num = num + 1;
			}
		}
		if (num > atomnum)
			break;
	}
	fclose(fpin);
	if ((*minfo).usercharge < -9990)
		(*minfo).icharge = intcharge(atomnum, atom);
}



/*CHARGE METHOD : READ CHARGE */
void rcharge(char *filename, int atomnum, ATOM atom[], CONTROLINFO cinfo,
			 MOLINFO *minfo)
{
	FILE *fpcharge;
	char line[MAXCHAR];
	int i;
	double tmpf1, tmpf2, tmpf3, tmpf4, tmpf5, tmpf6, tmpf7, tmpf8;

	if ((fpcharge = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "Cannot open charge file to read: %s , exit\n", filename);
		exit(1);
	}
	i = 0;
	for (;;) {
		if (fgets(line, MAXCHAR, fpcharge) == NULL) {
/*       printf("\nFinished reading file %s", filename); */
			break;
		}
		sscanf(line, "%lf%lf%lf%lf%lf%lf%lf%lf", &tmpf1, &tmpf2, 
			   &tmpf3, &tmpf4, &tmpf5, &tmpf6, &tmpf7, &tmpf8);
		if(i==atomnum) break;
		atom[i].charge = tmpf1;
		if(i==atomnum-1) break;
		atom[i+1].charge = tmpf2;
		if(i==atomnum-2) break;
		atom[i+2].charge = tmpf3;
		if(i==atomnum-3) break;
		atom[i+3].charge = tmpf4;
		if(i==atomnum-4) break;
		atom[i+4].charge = tmpf5;
		if(i==atomnum-5) break;
		atom[i+5].charge = tmpf6;
		if(i==atomnum-6) break;
		atom[i+6].charge = tmpf7;
		if(i==atomnum-7) break;
		atom[i+7].charge = tmpf8;
		i = i + 8;
	}
	fclose(fpcharge);
	if ((*minfo).usercharge < -9990.)
		(*minfo).icharge = intcharge(atomnum, atom);
}

/*CHARGE METHOD : WRITE CHARGE */

void wcharge(char *filename, int atomnum, ATOM atom[], CONTROLINFO cinfo,
			 MOLINFO minfo)
{
	FILE *fpcharge;
	int i;
	int num = 0;

	if ((fpcharge = fopen(filename, "w")) == NULL) {
		fprintf(stderr, "Cannot open charge file to write: %s , exit\n", filename);
		exit(1);
	}
	for (i = 0; i < atomnum; i++) {
		fprintf(fpcharge, "%10.6lf", atom[i].charge);
		num++;
		if (num == 8) {
			num = 0;
			fprintf(fpcharge, "\n");
		}
	}
	fclose(fpcharge);
}



/*CHARGE METHOD : RESP CHARGE */

void resp(char *filename, int atomnum, ATOM * atom, int bondnum,
		  BOND * bond, CONTROLINFO cinfo, MOLINFO minfo)
{
	char tmpchar1[MAXCHAR];
	char tmpchar2[MAXCHAR];
	char tmpchar3[MAXCHAR];
	char tmpchar4[MAXCHAR];
	char tmpchar5[MAXCHAR];
	size_t copied_size;
	int status = 0;

	if (strcmp(cinfo.intype, "gout") != 0
		&& strcmp(cinfo.intype, "11") != 0) {
		printf("\n Sorry, RESP charge needs a Gaussian output file (gout)");
		return;
	}
	wac("ANTECHAMBER_RESP.AC", atomnum, atom, bondnum, bond, cinfo, minfo);

        copied_size = build_exe_path(tmpchar1,
                "espgen -o ANTECHAMBER.ESP -i ", MAXCHAR );
        strncat(tmpchar1, filename, MAXCHAR - copied_size );
 
        copied_size = build_exe_path(tmpchar2, "respgen"
                " -i ANTECHAMBER_RESP.AC -o ANTECHAMBER_RESP1.IN -f resp1"
                , MAXCHAR );
 
        copied_size = build_exe_path(tmpchar3, "respgen"
                " -i ANTECHAMBER_RESP.AC -o ANTECHAMBER_RESP2.IN -f resp2"
                , MAXCHAR );
 
        copied_size = build_exe_path(tmpchar4, "resp", MAXCHAR );
 
        copied_size = build_exe_path(tmpchar5, "resp", MAXCHAR );

        /*  run espgen:  */
        if (cinfo.intstatus == 2) fprintf(stderr, "\nRunning: %s\n\n", tmpchar1);

	status = system(tmpchar1);
        if(status != 0) {
        	fprintf(stderr, "Error: cannot run \"%s\" in resp() of charge.c properly, exit\n", tmpchar1);
                exit(1);
        }

        /*  run  respgen to make ANTECHAMBER_RESP1.IN:   */
        if (cinfo.intstatus == 2) fprintf(stderr, "\nRunning: %s\n\n", tmpchar2);
        status = system(tmpchar2);
    	if(status != 0) {
        	fprintf(stderr, "Error: cannot run \"%s\" in resp() of charge.c properly, exit\n", tmpchar2);
        	exit(1);
    	}

        /*  run  respgen to make ANTECHAMBER_RESP2.IN:   */
        if (cinfo.intstatus == 2) fprintf(stderr, "\nRunning: %s\n\n", tmpchar3);
        status = system(tmpchar3);
    	if(status != 0) {
        	fprintf(stderr, "Error: cannot run \"%s\" in resp() of charge.c properly, exit\n", tmpchar3);
        	exit(1);
    	}
        /* run first stage of resp:  */
    	strcat(tmpchar4,
        	" -O -i ANTECHAMBER_RESP1.IN -o ANTECHAMBER_RESP1.OUT -e ANTECHAMBER.ESP -t qout");
        if (cinfo.intstatus == 2) fprintf(stderr, "\nRunning: %s\n\n", tmpchar4);
    	status = system(tmpchar4);
    	if(status != 0) {
        	fprintf(stderr, "Error: cannot run \"%s\" in resp() of charge.c properly, exit\n", tmpchar4 );
        	exit(1);
    	}

        /* run second stage of resp:  */
    	strcat(tmpchar5,
        	" -O -i ANTECHAMBER_RESP2.IN -o ANTECHAMBER_RESP2.OUT -e ANTECHAMBER.ESP -q qout -t QOUT");
        if (cinfo.intstatus == 2) fprintf(stderr, "\nRunning: %s\n\n", tmpchar5);
    	status = system(tmpchar5);
    	if(status != 0) {
        	fprintf(stderr, "Error: cannot run \"%s\" in resp() of charge.c properly, exit\n", tmpchar5 );
        	exit(1);
    	}
        rcharge("QOUT", atomnum, atom, cinfo, &minfo);
}

/* CHARGE METHOD : BCC-AM1 */
void bccharge(int atomnum, ATOM atom[], int bondnum, BOND bond[],
			  AROM arom[], CONTROLINFO *cinfo, MOLINFO *minfo)
{
	int i;
	int atomnum_tmp;
	int bondnum_tmp;
	int status = 0;
	ATOM *atom_tmp;
	BOND *bond_tmp;
	char tmpchar[MAXCHAR];
	char tmpchar1[MAXCHAR];
	char tmpchar2[MAXCHAR];
	size_t copied_size;

	wac("ANTECHAMBER_AM1BCC_PRE.AC", atomnum, atom, bondnum, bond, *cinfo,
		*minfo);
	copied_size = build_exe_path(tmpchar1, "am1bcc"
        	" -i ANTECHAMBER_AM1BCC_PRE.AC -o ANTECHAMBER_AM1BCC.AC -f ac -p " , MAXCHAR );
       	/* tmpchar2[0] = '\0';
       	strcpy(tmpchar2, amberhome);
        strcat(tmpchar2, "/dat/antechamber/BCCPARM.DAT"); */
        build_exe_path(tmpchar2, "../dat/antechamber/BCCPARM.DAT", MAXCHAR);
       	strncat(tmpchar1, tmpchar2, MAXCHAR - copied_size );
	strcat(tmpchar1, " -s ");
	sprintf(tmpchar2, "%d", (*cinfo).intstatus);
/*      newitoa((*cinfo).intstatus, tmpchar2); */
	strcat(tmpchar1, tmpchar2);
/*
	if (strcmp((*cinfo).intype, "ac") == 0 || strcmp((*cinfo).intype, "1") == 0
		|| strcmp((*cinfo).intype, "mol2") == 0
		|| strcmp((*cinfo).intype, "2") == 0) {
		if ((*cinfo).prediction_index == -1)
			(*cinfo).prediction_index = 0;
		sprintf(tmpchar, "%d", (*cinfo).prediction_index); 
		strcat(tmpchar1, " -j ");
		strcat(tmpchar1, tmpchar);
	} else
		strcat(tmpchar1, " -j 4");
*/
	if ((*cinfo).prediction_index == -1)
		(*cinfo).prediction_index = 0;
/*	newitoa((*cinfo).prediction_index, tmpchar); */
	sprintf(tmpchar, "%d", (*cinfo).prediction_index);
	strcat(tmpchar1, " -j ");
	strcat(tmpchar1, tmpchar);
    quote_full_line(tmpchar1, MAXCHAR);

	if ((*cinfo).intstatus == 2)
		fprintf(stderr, "\nRunning: %s\n", tmpchar1);
	status = system(tmpchar1);
        if(status != 0) {
                fprintf(stderr, "Error: cannot run \"%s\" in bccharge() of charge.c properly, exit\n", tmpchar1);
                exit(1);
        }
	(*cinfo).maxatom = atomnum + 10;
	(*cinfo).maxbond = bondnum + 10;
	atom_tmp = (ATOM *) malloc(sizeof(ATOM) * (*cinfo).maxatom);
	if (atom_tmp == NULL) {
		fprintf(stderr, "memory allocation error for *atom_tmp\n");
		exit(1);
	}
	bond_tmp = (BOND *) malloc(sizeof(BOND) * (*cinfo).maxbond);
	if (bond_tmp == NULL) {
		fprintf(stderr, "memory allocation error for *bond_tmp\n");
		exit(1);
	}
	for (i = 0; i < (*cinfo).maxbond; ++i) {
		bond_tmp[i].jflag = -1; /* bond type has not been assigned */
	}
	rac("ANTECHAMBER_AM1BCC.AC", &atomnum_tmp, atom_tmp, &bondnum_tmp,
		bond_tmp, cinfo, minfo);

	for (i = 0; i < atomnum; i++)
		atom[i].charge = atom_tmp[i].charge;
	free(atom_tmp);
	free(bond_tmp);

}

/* CHARGE METHOD : BCC-AM1 */
void bcc(char *filename, int atomnum, ATOM atom[], int bondnum,
		 BOND bond[], AROM arom[], CONTROLINFO *cinfo, MOLINFO *minfo)
{
	char tmpchar[MAXCHAR];
	int status = 0;	

	if (strcmp((*cinfo).intype, "mopout") == 0
		|| strcmp((*cinfo).intype, "12") == 0)
		rmopcharge(filename, atomnum, atom, minfo);
	else if (strcmp((*cinfo).intype, "divout") == 0
		|| strcmp((*cinfo).intype, "22") == 0)
		rdivcharge(filename, atomnum, atom, minfo,1);
/*	else if (strcmp((*cinfo).intype, "sqmout") == 0
 *		|| strcmp((*cinfo).intype, "32") == 0)
 *		rsqmcharge(filename, atomnum, atom, minfo);
 */
	else {
		if((*minfo).divcon == 0) {
			wmopcrt("mopac.in", atomnum, atom, *minfo);
		/*	strcpy(tmpchar, amberhome);
            strcat(tmpchar, "/bin/mopac.sh"); */
			build_exe_path(tmpchar, "mopac.sh", MAXCHAR);
		}
		if((*minfo).divcon == 1) {
			wdivcrt("divcon.in", atomnum, atom, *minfo);
		/*	strcpy(tmpchar, amberhome);
			strcat(tmpchar, "/bin/divcon"); */
			build_exe_path(tmpchar, "divcon", MAXCHAR);
		}
		if((*minfo).divcon == 2) {
			wsqmin("sqm.in", atomnum, atom, *minfo);
		/*	strcpy(tmpchar, amberhome);
			strcat(tmpchar, "/bin/sqm -O -i sqm.in -o sqm.out"); */
			build_exe_path(tmpchar, "sqm -O -i sqm.in -o sqm.out", MAXCHAR);
		}
		if ((*cinfo).intstatus == 2 || (*cinfo).intstatus == 1)
			fprintf(stderr, "\nRunning: %s\n", tmpchar);
		status = system(tmpchar);
		if(status != 0) {
			fprintf(stderr, "Error: cannot run \"%s\" of bcc() in charge.c properly, exit\n", tmpchar);
			exit(1);
		}
		if((*minfo).divcon == 0) {
			rmopcharge("mopac.out", atomnum, atom, minfo);
 			wpdb_optimized("mopac.out",atomnum,atom,0);
		}
		if((*minfo).divcon == 1) {
			rdivcharge("divcon.out", atomnum, atom, minfo, 1);
			wpdb_optimized("divcon.out",atomnum,atom,1);
		}
		if((*minfo).divcon == 2) {
			rsqmcharge("sqm.out", atomnum, atom, minfo);
			/* wpdb_optimized("divcon.out",atomnum,atom,2); */
		}
	}
	bccharge(atomnum, atom, bondnum, bond, arom, cinfo, minfo);
}

/* CHARGE METHOD : CM1 */
void cm1(int atomnum, ATOM atom[], CONTROLINFO *cinfo, MOLINFO *minfo)
{
        char tmpchar[MAXCHAR];
	int status = 0;

	if((*minfo).divcon == 0) {
        	wmopcrt("mopac.in", atomnum, atom, *minfo);
            /*  strcpy(tmpchar, amberhome);
                strcat(tmpchar, "/bin/mopac.sh"); */
                build_exe_path(tmpchar, "mopac.sh", MAXCHAR);
	}
	if((*minfo).divcon == 1) {
        	wdivcrt("divcon.in", atomnum, atom, *minfo);
             /* strcpy(tmpchar, amberhome);
                strcat(tmpchar, "/bin/divcon"); */
                build_exe_path(tmpchar, "divcon", MAXCHAR);
	}
         if ((*cinfo).intstatus == 2 || (*cinfo).intstatus == 1)
         	fprintf(stderr, "\nRunning: %s\n", tmpchar);
        status = system(tmpchar);
	if(status != 0) {
		fprintf(stderr, "Error: cannot run \"%s\" of cm1() in charge.c properly, exit\n", tmpchar);
		exit(1);
	}
	if((*minfo).divcon == 0) {
		rmopcharge("mopac.out", atomnum, atom, minfo);
		fprintf(stderr, "Warning: if the mopac program cannot generate the CM1 charges, Mulliken charges are generated\n");
 		wpdb_optimized("mopac.out",atomnum,atom,0);
	}

	if((*minfo).divcon == 1) {
		rdivcharge("divcon.out", atomnum, atom, minfo, 2);
		wpdb_optimized("divcon.out",atomnum,atom,1);
	}
}

/* CHARGE METHOD : CM2 */
void cm2(int atomnum, ATOM atom[], CONTROLINFO *cinfo, MOLINFO *minfo)
{
        char tmpchar[MAXCHAR];
	int status = 0;

	if((*minfo).divcon == 0) {
        	wmopcrt("mopac.in", atomnum, atom, *minfo);
             /* strcpy(tmpchar, amberhome);
                strcat(tmpchar, "/bin/mopac.sh"); */
                build_exe_path(tmpchar, "mopac.sh", MAXCHAR);
	}
	if((*minfo).divcon == 1) {
        	wdivcrt("divcon.in", atomnum, atom, *minfo);
             /* strcpy(tmpchar, amberhome);
                strcat(tmpchar, "/bin/divcon"); */
                build_exe_path(tmpchar, "divcon", MAXCHAR);
	}
         if ((*cinfo).intstatus == 2 || (*cinfo).intstatus == 1)
         	fprintf(stderr, "\nRunning: %s\n", tmpchar);
        status = system(tmpchar);
	if(status != 0) {
		fprintf(stderr, "Error: cannot run \"%s\" of cm2() in charge.c properly, exit\n", tmpchar);
		exit(1);
	}
	if((*minfo).divcon == 0) {
		rmopcharge("mopac.out", atomnum, atom, minfo);
		fprintf(stderr, "Warning: if the mopac program cannot generate the CM2 charges, Mulliken charges are generated\n");
 		wpdb_optimized("mopac.out",atomnum,atom,0);
	}
	if((*minfo).divcon == 1) {
		rdivcharge("divcon.out", atomnum, atom, minfo, 3);
		wpdb_optimized("divcon.out",atomnum,atom,1);
	}
}

/* CHARGE METHOD : Mulliken  */
void mul(char *filename, int atomnum, ATOM atom[], CONTROLINFO *cinfo,
		 MOLINFO *minfo)
{
	char tmpchar[MAXCHAR];
	int status = 0;

	if (strcmp((*cinfo).intype, "gout") == 0
		|| strcmp((*cinfo).intype, "11") == 0)
		rgaucharge(filename, "mul", atomnum, atom, minfo);
	else if (strcmp((*cinfo).intype, "mopout") == 0
			 || strcmp((*cinfo).intype, "12") == 0)
		rmopcharge(filename, atomnum, atom, minfo);
	else if (strcmp((*cinfo).intype, "divout") == 0
			 || strcmp((*cinfo).intype, "22") == 0)
		rdivcharge(filename, atomnum, atom, minfo,1);
	else {
		if((*minfo).divcon == 0) {
                        wmopcrt("mopac.in", atomnum, atom, *minfo);
                     /* strcpy(tmpchar, amberhome);
                        strcat(tmpchar, "/bin/mopac.sh"); */
                        build_exe_path(tmpchar, "mopac.sh", MAXCHAR);
                }
		if((*minfo).divcon == 1) {
                        wdivcrt("divcon.in", atomnum, atom, *minfo);
                     /* strcpy(tmpchar, amberhome);
                        strcat(tmpchar, "/bin/divcon"); */
                        build_exe_path(tmpchar, "divcon", MAXCHAR);
                }
                if ((*cinfo).intstatus == 2 || (*cinfo).intstatus == 1)
                        fprintf(stderr, "\nRunning: %s\n", tmpchar);
                status = system(tmpchar);
		if(status != 0) {
			fprintf(stderr, "Error: cannot run \"%s\" of mul() in charge.c properly, exit\n", tmpchar);
			exit(1);
		}
		if((*minfo).divcon == 0) {
                        rmopcharge("mopac.out", atomnum, atom, minfo);
 			wpdb_optimized("mopac.out",atomnum,atom,0);
		}
		if((*minfo).divcon == 1) {
                        rdivcharge("divcon.out", atomnum, atom, minfo, 1);
			wpdb_optimized("divcon.out",atomnum,atom,1);
		}
        }
}

/* CHARGE METHOD : ESP(Kollman) */
void esp(char *filename, int atomnum, ATOM atom[], CONTROLINFO cinfo,
		 MOLINFO minfo)
{
	if (strcmp(cinfo.intype, "gout") == 0)
		rgaucharge(filename, "esp", atomnum, atom, &minfo);
	else
		printf("\n Kollman ESP charges IS Only Used for Gaussian Output");
}

/* CHARGE METHOD : Gasteiger charge */
void rgasparm(char *filename, int *gasparmnum, GASTEIGER gas[])
{
	FILE *fp;
	int num = 0;
	char line[MAXCHAR];

	if ((fp = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "Cannot open the Gasteiger parameter file: %s, exit\n", filename);
		return;
	}
	for (;;) {
		if (fgets(line, MAXCHAR, fp) == NULL)
			break;
		if (strncmp(line, "GASPARM", 7) == 0) {
			gas[num].name[0] = '\0';
			sscanf(&line[8], "%s%lf%lf%lf%lf%lf", gas[num].name, &gas[num].a,
				   &gas[num].b, &gas[num].c, &gas[num].d, &gas[num].charge);
			num++;
		}
	}
	*gasparmnum = num;
	fclose(fp);
}
void assign(int atomnum, ATOM atom[], int gasparmnum, GASTEIGER gas[],
			int gasparmindex[])
{
	int i, j;
	int flag;
	for (i = 0; i < atomnum; i++) {
		flag = 0;
		for (j = 0; j < gasparmnum; j++)
			if (strcmp(gas[j].name, atom[i].ambername) == 0) {
				gasparmindex[i] = j;
				flag = 1;
				break;
			}
		if (flag == 0) {
			fprintf(stderr,
					"\nNo Gasteiger parameter for atom[%d]:%s:%s, exit\n",
					i, atom[i].name, atom[i].ambername);
			exit(1);
		}
	}
}

double rmscal(int atomnum, double gaschargep[], double gaschargea[])
{
	double total = 0.0;
	double rms;
	int i;

	for (i = 0; i < atomnum; i++) {
		total +=
			(gaschargep[i] - gaschargea[i]) * (gaschargep[i] -
											   gaschargea[i]);
		gaschargep[i] = gaschargea[i];
	}
	rms = total / atomnum;
	rms = sqrt(rms);
	return rms;
}

void gasiter(int atomnum, ATOM atom[], int gasparmnum, GASTEIGER gas[], int netcharge)
{
	int i, j;
	int iteration = 0;
	double q;
	double xx;
	double rmsd;
	int *gasparmindex;
	double *x;
	double *gaschargep;
	double *gaschargea;
	double TotalGasCharge = 0.0;

	gasparmindex = (int *) malloc(sizeof(int) * atomnum);
	if (gasparmindex == NULL) {
		fprintf(stderr,
				"memory allocation error for gasparmindex in gasiter(), exit\n");
		exit(1);
	}
	x = (double *) malloc(sizeof(double) * atomnum);
	if (x == NULL) {
		fprintf(stderr,
				"memory allocation error for x in gasiter(), exit\n");
		exit(1);
	}
	gaschargep = (double *) malloc(sizeof(double) * atomnum);
	if (gaschargep == NULL) {
		fprintf(stderr,
				"memory allocation error for gaschargep in gasiter(), exit\n");
		exit(1);
	}
	gaschargea = (double *) malloc(sizeof(double) * atomnum);
	if (gaschargea == NULL) {
		fprintf(stderr,
				"memory allocation error for gaschargea in gasiter(), exit\n");
		exit(1);
	}
	assign(atomnum, atom, gasparmnum, gas, gasparmindex);
	for (i = 0; i < atomnum; i++) {
		gaschargep[i] = gas[gasparmindex[i]].charge; 
		gaschargea[i] = gas[gasparmindex[i]].charge;
		TotalGasCharge +=  gas[gasparmindex[i]].charge;
	}
	if(TotalGasCharge != netcharge) {
		printf("\nThe net charge of the molecule (%d) does not equal to the total charge (%-5.2lf) based on Gasteiger atom type, exit\n", netcharge, TotalGasCharge);
		exit(1);
	}
	do {
		for (i = 0; i < atomnum; i++) {
			x[i] =
				gas[gasparmindex[i]].a +
				gas[gasparmindex[i]].b * gaschargep[i];
			x[i] += gas[gasparmindex[i]].c * gaschargep[i] * gaschargep[i];
			if (x[i] == 0.0)
				x[i] = 0.0000000001;
		}
		for (i = 0; i < atomnum; i++)
			for (j = i + 1; j < atomnum; j++)
				if (atom[i].con[0] == j || atom[i].con[1] == j
					|| atom[i].con[2] == j || atom[i].con[3] == j
					|| atom[i].con[4] == j || atom[i].con[5] == j) {
					if (x[i] <= x[j]) {
						xx = gas[gasparmindex[i]].d ;
						q = (x[j] - x[i]) / xx * pow(DAMPFACTOR,
													 iteration + 1);
						gaschargea[i] += q;
						gaschargea[j] -= q;
					}
					if (x[i] > x[j]) {
						xx = gas[gasparmindex[j]].d; 
						q = (x[i] - x[j]) / xx * pow(DAMPFACTOR,
													 iteration + 1);
						gaschargea[i] -= q;
						gaschargea[j] += q;
					}
				}
		iteration++;
/*
		printf("\nIteration %5d", iteration);
		for (i = 0; i < atomnum; i++) 
			printf("\n%5d %5s %8.4lf %8.4lf", i+1, atom[i].name, gaschargep[i], gaschargea[i]);
*/
		rmsd = rmscal(atomnum, gaschargep, gaschargea);
	} while (rmsd > CONVERG && iteration < GASMAXITER);

	for (i = 0; i < atomnum; i++)
		atom[i].charge = gaschargea[i];
	free(gasparmindex);
	free(x);
	free(gaschargep);
	free(gaschargea);
}

void gascharge(int atomnum, ATOM atom[], int bondnum, BOND bond[],
			   CONTROLINFO cinfo, MOLINFO *minfo)
{
	int i;
	int atomnum_tmp;
	int bondnum_tmp;
	int status = 0;
	ATOM *atom_tmp;
	BOND *bond_tmp;

	int gasparmnum = 0;
	char tmpchar1[MAXCHAR];
	char tmpchar2[MAXCHAR];
	size_t copied_size;	
	GASTEIGER *gas;

	gas = (GASTEIGER *) malloc(sizeof(GASTEIGER) * MAXGAS);
	if (gas == NULL) {
		fprintf(stderr,
				"memory allocation error for gas in gascharge(),increase GASMAX in define.h exit\n");
		exit(1);
	}
	(*minfo).gfilename[0] = '\0';
	strcpy((*minfo).gfilename, amberhome);
	strcat((*minfo).gfilename, "/dat/antechamber/GASPARM.DAT");
	rgasparm((*minfo).gfilename, &gasparmnum, gas);	/*principle parameter file */
	wac("ANTECHAMBER_GAS.AC", atomnum, atom, bondnum, bond, cinfo, (*minfo));
        copied_size = build_exe_path(tmpchar1, "atomtype"
                " -i ANTECHAMBER_GAS.AC -o ANTECHAMBER_GAS_AT.AC -d "
                , MAXCHAR );
        tmpchar2[0] = '\0';
     /*   strcpy(tmpchar2, amberhome);
        strcat(tmpchar2, "/dat/antechamber/ATOMTYPE_GAS.DEF"); */
        build_exe_path(tmpchar2, "../dat/antechamber/ATOMTYPE_GAS.DEF", MAXCHAR);
        strncat(tmpchar1, tmpchar2, MAXCHAR - copied_size );
        quote_full_line(tmpchar1, MAXCHAR);

	if (cinfo.intstatus == 2)
		fprintf(stderr, "Running: %s\n", tmpchar1);
	status = system(tmpchar1);
	if(status != 0) {
        	fprintf(stderr, "Error: cannot run \"%s\" of gascharge() in charge.c properly, exit\n", tmpchar1);
                exit(1);
	}
	cinfo.maxatom = atomnum + 10;
	cinfo.maxbond = bondnum + 10;
	atom_tmp = (ATOM *) malloc(sizeof(ATOM) * cinfo.maxatom);
	if (atom_tmp == NULL) {
		fprintf(stderr, "memory allocation error for *atom_tmp\n");
		exit(1);
	}
	bond_tmp = (BOND *) malloc(sizeof(BOND) * cinfo.maxbond);
	if (bond_tmp == NULL) {
		fprintf(stderr, "memory allocation error for *bond_tmp\n");
		exit(1);
	}
	for (i = 0; i < cinfo.maxbond; ++i) {
		bond_tmp[i].jflag = -1; /* bond type has not been assigned */
	}

	rac("ANTECHAMBER_GAS_AT.AC", &atomnum_tmp, atom_tmp, &bondnum_tmp,
		bond_tmp, &cinfo, minfo);
	gasiter(atomnum_tmp, atom_tmp, gasparmnum, gas, (*minfo).usercharge);
	for (i = 0; i < atomnum; i++)
		atom[i].charge = atom_tmp[i].charge;
	free(atom_tmp);
	free(bond_tmp);
	free(gas);
}



void write_sybyl_bat(char *str)
{
	FILE *fpout;
	char tmpchar[MAXCHAR];
	int status = 0;

	if ((fpout = fopen("antechamber_sybyl.bat", "w")) == NULL) {
		fprintf(stderr, "Cannot open antechamber_sybyl.bat , exit\n");
		exit(1);
	}
	strcpy(tmpchar, amberhome);
	strcat(tmpchar, "/dat/antechamber/charge.spl");

	fprintf(fpout, "%s", "#!/bin/csh");
	fprintf(fpout, "\n%s", "sybyl << @");
	fprintf(fpout, "\n%s", "take ");
	fprintf(fpout, "%s", tmpchar);
	if (strcmp(str, "gas1") == 0)
		fprintf(fpout, "\n%s", "CHARGE GASTEIGER");
	if (strcmp(str, "gas2") == 0)
		fprintf(fpout, "\n%s", "CHARGE GAST_HUCK");
	if (strcmp(str, "del") == 0)
		fprintf(fpout, "\n%s", "CHARGE DELRE");
	if (strcmp(str, "pull") == 0)
		fprintf(fpout, "\n%s", "CHARGE PULLMAN ");
	if (strcmp(str, "huc") == 0)
		fprintf(fpout, "\n%s", "CHARGE HUCKEL");
	if (strcmp(str, "mmff") == 0)
		fprintf(fpout, "\n%s", "CHARGE MMFF94");
	fprintf(fpout, "\n%s", "QUIT YES");
	fprintf(fpout, "\n%s\n\n", "@");
	fclose(fpout);
	status = system("chmod +x antechamber_sybyl.bat");
	if(status != 0) {
                fprintf(stderr, "Error: cannot run \"%s\" of write_sybyl_bat() in charge.c properly, exit\n", "chmod +x antechamber_sybyl.bat");
                exit(1);
        }
	status = system("antechamber_sybyl.bat");
	if(status != 0) {
                fprintf(stderr, "Error: cannot run \"%s\" of write_sybyl_bat() in charge.c properly, exit\n", "antechamber_sybyl.bat");
                exit(1);
        }

}


/* CHARGE METHOD : Gasteiger, using sybyl */
void gas1(int atomnum, ATOM atom[], int bondnum, BOND bond[], AROM arom[],
		  CONTROLINFO cinfo, MOLINFO minfo)
{
	wmol2("ANTECHAMBER_INPUT.mol2", atomnum, atom, bondnum, bond, arom,
		  cinfo, minfo);
	write_sybyl_bat("gas1");
	rmol2("ANTECHAMBER_INPUT.mol2", &atomnum, atom, &bondnum, bond, &cinfo,
		  &minfo, 0);
}

/* CHARGE METHOD : Del Re, using sybyl*/
void del(int atomnum, ATOM atom[], int bondnum, BOND bond[], AROM arom[],
		 CONTROLINFO cinfo, MOLINFO minfo)
{
	wmol2("ANTECHAMBER_INPUT.mol2", atomnum, atom, bondnum, bond, arom,
		  cinfo, minfo);
	write_sybyl_bat("del");
	rmol2("ANTECHAMBER_INPUT.mol2", &atomnum, atom, &bondnum, bond, &cinfo,
		  &minfo, 0);
}

/* CHARGE METHOD : Pullman, using sybyl */
void pull(int atomnum, ATOM atom[], int bondnum, BOND bond[], AROM arom[],
		  CONTROLINFO cinfo, MOLINFO minfo)
{
	wmol2("ANTECHAMBER_INPUT.mol2", atomnum, atom, bondnum, bond, arom,
		  cinfo, minfo);
	write_sybyl_bat("pull");
	rmol2("ANTECHAMBER_INPUT.mol2", &atomnum, atom, &bondnum, bond, &cinfo,
		  &minfo, 0);
}

/* CHARGE METHOD :Gasteiger-Huckel, using sybyl*/
void gas2(int atomnum, ATOM atom[], int bondnum, BOND bond[], AROM arom[],
		  CONTROLINFO cinfo, MOLINFO minfo)
{
	wmol2("ANTECHAMBER_INPUT.mol2", atomnum, atom, bondnum, bond, arom,
		  cinfo, minfo);
	write_sybyl_bat("gas2");
	rmol2("ANTECHAMBER_INPUT.mol2", &atomnum, atom, &bondnum, bond, &cinfo,
		  &minfo, 0);
}

/* CHARGE METHOD : Huckel, using sybyl */
void huc(int atomnum, ATOM atom[], int bondnum, BOND bond[], AROM arom[],
		 CONTROLINFO cinfo, MOLINFO minfo)
{

	wmol2("ANTECHAMBER_INPUT.mol2", atomnum, atom, bondnum, bond, arom,
		  cinfo, minfo);
	write_sybyl_bat("huc");
	rmol2("ANTECHAMBER_INPUT.mol2", &atomnum, atom, &bondnum, bond, &cinfo,
		  &minfo, 0);
}

/* CHARGE METHOD :MMFF94, using sybyl*/
void mmff(int atomnum, ATOM atom[], int bondnum, BOND bond[], AROM arom[],
		  CONTROLINFO cinfo, MOLINFO minfo)
{
	wmol2("ANTECHAMBER_INPUT.mol2", atomnum, atom, bondnum, bond, arom,
		  cinfo, minfo);
	write_sybyl_bat("mmff");
	rmol2("ANTECHAMBER_INPUT.mol2", &atomnum, atom, &bondnum, bond, &cinfo,
		  &minfo, 0);
}
