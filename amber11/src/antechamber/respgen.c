/*
************************************************************************
*           All Copyright Reserved!                                    *
*                                                                      *
*  Prog:    respgen                                                    *
*  Version: version 1.0                                                *
*  Author:  Junmei Wang                                                *
*                                                                      *
*  Department of Pharmaceutical Chemistry                              *
*  School of Pharmacy                                                  *
*  University of California                                            *
*  San Francisco   CA 94143                                            *
*  Octomber, 2001                                                      *
************************************************************************
*/
char *amberhome;
# include "common.h"
# include "define.h"
# include "atom.h"
# include "utility.c"
# include "common.c"
# include "ac.c"
# define MAXPATHATOMNUM 5000
# define MAX_RESP_ATOM 1000

ATOM *atom;
BOND *bond;
int atomnum = 0;
int bondnum = 0;
CONTROLINFO cinfo;
MOLINFO minfo;
int *selectindex;
int *equatomno;
int *pathnum;
int *selectelement;
int *pathatomnum;
int maxlength = -1;
double *pathscore[MAX_RESP_ATOM] ;
int i, j, k, l;
FILE *fpin;
FILE *fpout;
FILE *fpcharge;
FILE *fpaddinfo;
char line[MAXCHAR];
char ifilename[MAXCHAR];
char ofilename[MAXCHAR];
char afilename[MAXCHAR];

int method = 0;
int selectnum = 0;
int pathnumindex = 0;
int atomindex = 0;
double charge = 0.0;

int iaddinfo = 0;
double *pcharge; /* predefined partial charge */
int numcharge = 0;
int *icharge; /*flag of predefined partial charge */
int nconf = 1;

void readinfo() {
int i,j;
int tmpint;
int num;
double tmpf;
numcharge = 0;
for (;;) {
	if (fgets(line, MAXCHAR, fpaddinfo) == NULL) break;
        if (strncmp("CHARGE", line, 6) == 0) {
		sscanf(&line[6], "%lf%d", &tmpf, &tmpint); 
		if(tmpint < 1) {
			fprintf(stderr, "\nAtom ID should not smaller than 1, exit");
			exit(1);	
		}
		pcharge[tmpint-1] = tmpf;
		icharge[tmpint-1] = 1;
		numcharge ++;
        }
}
num = 0;
for (i = 0; i < nconf; i++)
	for (j = 0; j < atomnum; j++) {
        	fprintf(fpcharge, "%10.6lf", pcharge[j]);
        	num++;
        	if (num == 8) {
                	num = 0;
                	fprintf(fpcharge, "\n");
        	}
	}

}

void scorepath(ATOM atm[], int selectnum, int startnum)
{
	int i, j, k;
	int start;
	double score;
	int resetindex;
	start = -1;
	resetindex = -1;
	selectindex[startnum] = selectnum;
	selectelement[selectnum++] = startnum;
	if(maxlength != -1 && selectnum > maxlength)
		return;
	for (i = 0; i < 6; i++) {
		if (atm[startnum].con[i] == -1) {
			score = 0.0;
			for (j = 0; j < selectnum; j++) {
				/* printf("%5d", selectelement[j]); */
				score +=
					(j + 1) * 0.11 +
					atom[selectelement[j]].atomicnum * 0.08;
			}
			pathscore[atomindex][pathnumindex++] = score;
			if (pathnumindex >= pathatomnum[atomindex]) {
				pathatomnum[atomindex] += MAXPATHATOMNUM;
				fprintf
					(stderr, "\nInfo: the number of the path atoms exceeds MAXPATHATOMNUM(%d) for atom[%d],extend the size and reallocate the memory automatically",
					 pathatomnum[atomindex], atomindex);
				pathscore[atomindex] =
					(double *) realloc(pathscore[atomindex],
									   pathatomnum[atomindex] *
									   sizeof(double));
				if (pathscore[atomindex] == NULL) {
					fprintf(stderr,
							" reallocate memory for pathscore[%d] failed\n",
							atomindex);
					exit(1);
				}
			}
			/*    printf("\n %5d%8.4lf", selectnum, score); */
			return;
		}
		start = atm[startnum].con[i];
		for (k = 0; k < selectnum; k++)
			if (start == selectelement[k]) {
				resetindex = 1;
				break;
			}
		if (resetindex == 1) {
			resetindex = -1;
			continue; 
		}
		if (start == -1)
			return;
		scorepath(atm, selectnum, start);
		/* we have already visited this atom */
	}
}
void sort(double array[], int elemnum)
{
	int i, j;
	double tmp;
	for (i = 0; i < elemnum; i++)
		for (j = i + 1; j < elemnum; j++) {
			if (array[j] < array[i]) {
				tmp = array[i];
				array[i] = array[j];
				array[j] = tmp;
			}
		}
/* printf("\n"); */
/*
  printf("%8d", elemnum);  
  for(i=0;i<elemnum;i++)
  printf("%8.4lf\n", array[i]);  
*/
}

void equatom(void)
{
	int i, j, k;
	int equindex;
	double sum;
	for (i = 0; i < atomnum; i++)
		equatomno[i] = -1;
	for (i = 0; i < atomnum; i++) {
/* if(i!=1&&i!=3) continue;*/
		selectnum = 0;
		pathnumindex = 0;
		atomindex = i;
		for (j = 0; j < atomnum; j++) {
			selectindex[j] = -1;
			selectelement[i] = -1;
		}
		scorepath(atom, 0, i);
		pathnum[i] = pathnumindex;
	}
	selectnum = 0;
	for (i = 0; i < atomnum; i++) {
/*if(i!=2&&i!=4) continue;*/
		sum = 0.0;
/* printf("%d", pathnum[i]);*/
		for (j = 0; j < pathnum[i]; j++) {
			sum += pathscore[i][j];
			/* printf("%7.3lf",pathscore[i][j]); */
		}
/*  printf("%7.3lf\n", sum); */
	}
	for (i = 0; i < atomnum; i++) {
/*
printf("\n%s", atom[i].name);
*/
		sort(pathscore[i], pathnum[i]);
	}
	for (i = 0; i < atomnum; i++) {
		for (j = i + 1; j < atomnum; j++) {
			if (equatomno[j] != -1)
				continue;
			equindex = 1;
			if (pathnum[i] != pathnum[j])
				continue;
			for (k = 0; k < pathnum[i]; k++)
				if (pathscore[i][k] != pathscore[j][k]) {
					equindex = -1;
					break;
				}
			if (equindex == 1)
				equatomno[j] = i;
		}
	}
/*     for(i=0;i<atomnum;i++) printf("\n %5d%5d", i+1, equatomno[i]+1); */
}
void respin(int method)
{
	int i, j;
	int cindex;
	int *chindex;
	int times = 0;
	int tmpint;
	int count;
	int tcount;
	int tmpnum;
	int igroup = 0;

	double tmpcharge;
	FILE *fpout;

	chindex = (int *) malloc(sizeof(int) * (atomnum + 10));
        if (chindex == NULL) {
                fprintf(stderr, "memory allocation error for *chindex in respin()\n");
                exit(1);
        }
	
	if ((fpout = fopen(ofilename, "w")) == NULL) {
		printf("\n Cannot open file %s to write in respin(), exit", ofilename);
		return;
	}
	for (i = 0; i < atomnum; i++)
		chindex[i] = -1;
	for (i = 0; i < atomnum; i++) {
		cindex = 0;
		if (atom[i].name[0] == 'C')
			for (j = 0; j < 6; j++)
				if (atom[i].con[j] >= 0)
					if (atom[atom[i].con[j]].atomicnum == 1)
						cindex++;
		if (cindex >= 2) {
			chindex[i] = 1;
			for (j = 0; j < 6; j++)
				if (atom[i].con[j] >= 0)
					if (atom[atom[i].con[j]].atomicnum == 1) 
						chindex[atom[i].con[j]] = 1;
		}
	}
/* for(i=0;i<atomnum;i++)
 printf("\n %5d%5s%5d%5d", i+1, atom[i].name, chindex[i], equatomno[i]);
 */
	fprintf(fpout, "%s\n\n", "Resp charges for organic molecule");
	fprintf(fpout, "%s\n\n", " &cntrl");
	fprintf(fpout, " nmol = %d,\n", nconf);
	fprintf(fpout, "%s\n", " ihfree = 1,");
	fprintf(fpout, "%s\n", " ioutopt = 1,");
	if (method == 1 && numcharge > 0) 
		fprintf(fpout, "%s\n", " iqopt = 2,");
	if (method == 2) {
		fprintf(fpout, "%s\n", " iqopt = 2,");
		fprintf(fpout, "%s\n", " qwt = 0.001,");
	}
	if (method == 0) {
		fprintf(fpout, "%s\n", " idqrtrnt = 0,");
		fprintf(fpout, "%s\n", " dmscale = 1,");
	}
	fprintf(fpout, "\n%s\n", " &end");
	while(times < nconf) {
		fprintf(fpout, "%s\n", "    1.0");
		fprintf(fpout, "%s\n", "Resp charges for organic molecule");
		fprintf(fpout, "%5d%5d\n", (int) charge, atomnum);
		if (method == 0) {
			for (i = 0; i < atomnum; i++) {
				if(iaddinfo == 1 && icharge[i] == 1)
					fprintf(fpout, "%5d%5d", atom[i].atomicnum, -99);
				else {
					if (equatomno[i] == -1)
						fprintf(fpout, "%5d%5d", atom[i].atomicnum, 0);
					else
						fprintf(fpout, "%5d%5d", atom[i].atomicnum,
								equatomno[i] + 1);
				}
				if (chindex[i] == 1 && atom[i].atomicnum == 6)
					fprintf(fpout, "%8.4lf\n", 0.001);
				else
					fprintf(fpout, "\n");
			}
		}
		if (method == 1) {
			for (i = 0; i < atomnum; i++) {
				if(iaddinfo == 1 && icharge[i] == 1)
					fprintf(fpout, "%5d%5d\n", atom[i].atomicnum, -99);
				else {
					if (equatomno[i] == -1)
						fprintf(fpout, "%5d%5d\n", atom[i].atomicnum, 0);
					else if (chindex[i] == 1)
						fprintf(fpout, "%5d%5d\n", atom[i].atomicnum, 0);
					else
						fprintf(fpout, "%5d%5d\n", atom[i].atomicnum,
								equatomno[i] + 1);
				}
			}
		}

		if (method == 2) {
			for (i = 0; i < atomnum; i++) {
				if(iaddinfo == 1 && icharge[i] == 1)
					fprintf(fpout, "%5d%5d\n", atom[i].atomicnum, -99);
				else {
					if (chindex[i] != 1)
						fprintf(fpout, "%5d%5d\n", atom[i].atomicnum, -99);
					else if (equatomno[i] == -1)
						fprintf(fpout, "%5d%5d\n", atom[i].atomicnum, 0);
					else
						fprintf(fpout, "%5d%5d\n", atom[i].atomicnum,
								equatomno[i] + 1);
				}
			}
		}
		fprintf(fpout, "\n");
		times++;
   	}
	if(iaddinfo == 1) {
/* group */
		rewind(fpaddinfo);
		tmpnum = -1;
		for (;;) {
        		if (fgets(line, MAXCHAR, fpaddinfo) == NULL) break;
        		if (strncmp("GROUP", line, 5) == 0) {
				if(tmpnum > 0) {
					if(tmpnum != tcount) {
						fprintf(stderr,"The number of ATOM lines does not equal to the number defined in the GROUP line\n"); 
						exit(1);
					}
					else
						fprintf(fpout, "\n"); 
				} 
                		sscanf(&line[5], "%d%lf", &tmpnum, &tmpcharge);
				fprintf(fpout, "%5d%8.3lf\n", tmpnum, tmpcharge); 
				count = 0;
				tcount = 0;
				igroup++;
				continue;
			}
        		if (strncmp("ATOM", line, 4) == 0) {
                		sscanf(&line[4], "%d", &tmpint);
				fprintf(fpout, "%5d%5d", 1, tmpint); 
				count ++;
				tcount ++;
				if(count >= 8) {
					count = 0;
					fprintf(fpout, "\n");
				}	
			}
        	}

		if(tmpnum > 0) {
			if(tmpnum != tcount) {
				fprintf(stderr,"The number of ATOM lines does not equal to the number defined in the GROUP line\n"); 
				exit(1);
			}
			else
				fprintf(fpout, "\n\n"); 
		} 
/* equaliation*/
/* need one more blank line when no group is defined */
		if(igroup == 0) fprintf(fpout, "\n");
		for(i=1; i< nconf; i++)
			for(j=0; j< atomnum ; j++) {
				fprintf(fpout, "%5d\n", 2); 
				fprintf(fpout, "%5d%5d%5d%5d\n", 1, j+1, i+1, j+1); 
			}
	}
	if(iaddinfo == 0 && nconf > 1) {
		fprintf(fpout, "\n\n");
		for(i=1; i< nconf; i++)
			for(j=0; j< atomnum ; j++) {
				fprintf(fpout, "%5d\n", 2); 
				fprintf(fpout, "%5d%5d%5d%5d\n", 1, j+1, i+1, j+1); 
			}
	}
	free(chindex);
	fprintf(fpout, "\n");
	fclose(fpout);
}

int main(int argc, char *argv[])
{

	int overflow_flag = 0;			/*if overflow_flag ==1, reallocate memory */

    amberhome = (char *) getenv("AMBERHOME");
    if( amberhome == NULL ){
       fprintf( stderr, "AMBERHOME is not set!\n" );
       exit(1);
    }
	if (strcmp(COLORTEXT, "YES") == 0 || strcmp(COLORTEXT, "yes") == 0) {
		if (argc == 2
			&& (strcmp(argv[1], "-h") == 0
				|| strcmp(argv[1], "-H") == 0)) {
			printf("[31mUsage: respgen -i[0m input file name(ac)\n"
				   "[31m               -o[0m output file name\n"
				   "[31m               -l[0m maximum path length (default is -1, only recommand if the default setting\n"
				   "                  taking long time and/or cause core dump. If applied, a value of 8 to 10 should good)\n"
				   "[31m               -f[0m output file format (resp1 or resp2) \n"
				   "[32m                  resp1[0m - first stage resp fitting \n"
				   "[32m                  resp2[0m - second stage resp fitting\n"
				   "[31m               -a[0m additional input data (predefined charges, atom groups etc).)\n"
				   "[31m               -n[0m number of conformations (default is 1)\n");
			exit(1);
		}
		if (argc != 7 && argc != 9 && argc != 11 && argc != 13) {
			printf("[31mUsage: respgen -i[0m input file name(ac)\n"
				   "[31m               -o[0m output file name\n"
				   "[31m               -l[0m maximum path length (default is -1, only recommand if the default setting\n"
				   "                  taking long time and/or cause core dump. If applied, a value of 8 to 10 should good)\n"
				   "[31m               -f[0m output file format (resp1 or resp2) \n"
				   "[32m                  resp1[0m - first stage resp fitting \n"
				   "[32m                  resp2[0m - second stage resp fitting\n"
				   "[31m               -a[0m additional input data (predefined charges, atom groups etc).)\n"
				   "[31m               -n[0m number of conformations (default is 1)\n");
			exit(1);
		}
	} else {
		if (argc == 2
			&& (strcmp(argv[1], "-h") == 0
				|| strcmp(argv[1], "-H") == 0)) {
			printf("Usage respgen -i input file name(ac)\n");
			printf("              -o output file name\n");
			printf("	      -l maximum path length (default is -1, only recommand if the default setting\n");
			printf("                 taking long time and/or cause core dump. If applied, a value of 8 to 10 should good)\n");
			printf("              -f output file format (resp1 or resp2)\n");
			printf("                 resp1 - first stage resp fitting\n");
			printf("                 resp2 - second stage resp fitting \n");
		        printf("              -a additional input data (predefined charges, atom groups etc).)\n");
		        printf("              -n number of conformations (default is 1)\n");
			exit(1);
		}
		if (argc != 7 && argc != 9 && argc != 11) {
			printf("Usage respgen -i input file name(ac)\n");
			printf("              -o output file name\n");
			printf("	      -l maximum path length (default is -1, only recommand if the default setting\n");
			printf("                 taking long time and/or cause core dump. If applied, a value of 8 to 10 should good)\n");
			printf("              -f output file format (resp1 or resp2)\n");
			printf("                 resp1 - first stage resp fitting\n");
			printf("                 resp2 - second stage resp fitting \n");
		        printf("              -a additional input data (predefined charges, atom groups etc).)\n");
			exit(1);
		}
	}

	method = 0;
	maxlength = -1;
	for (i = 1; i < argc; i += 2) {
		if (strcmp(argv[i], "-i") == 0)
			strcpy(ifilename, argv[i + 1]);
		if (strcmp(argv[i], "-o") == 0)
			strcpy(ofilename, argv[i + 1]);
		if (strcmp(argv[i], "-a") == 0) {
			strcpy(afilename, argv[i + 1]);
			iaddinfo = 1;
		}
		if (strcmp(argv[i], "-n") == 0)
			nconf = atoi(argv[i+1]);
		if (strcmp(argv[i], "-l") == 0)
			maxlength = atoi(argv[i+1]); 
		if (strcmp(argv[i], "-f") == 0) {
			if (strcmp("resp", argv[i + 1]) == 0)
				method = 0;
			if (strcmp("resp1", argv[i + 1]) == 0)
				method = 1;
			if (strcmp("resp2", argv[i + 1]) == 0)
				method = 2;
		}
	}
	if(nconf < 1) {
		printf("\nNumber of conformations must be equal to or larger than 1"); 
		exit(1);
	}

	default_minfo(&minfo);
	default_cinfo(&cinfo);

	atom = (ATOM *) malloc(sizeof(ATOM) * cinfo.maxatom);
	if (atom == NULL) {
		fprintf(stderr, "memory allocation error for *atom\n");
		exit(1);
	}

	bond = (BOND *) malloc(sizeof(BOND) * cinfo.maxbond);
	if (bond == NULL) {
		fprintf(stderr, "memory allocation error for *bond\n");
		exit(1);
	}
	int i;
	for (i = 0; i < cinfo.maxbond; ++i) {
		bond[i].jflag = -1; /* bond type has not been assigned */
	}

	overflow_flag =
		rac(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);

	if (overflow_flag) {
		cinfo.maxatom = atomnum + 10;
		cinfo.maxbond = bondnum + 10;
		free(atom);
		free(bond);
		atom = (ATOM *) malloc(sizeof(ATOM) * cinfo.maxatom);
		if (atom == NULL) {
			fprintf(stderr, "memory allocation error for *atom\n");
			exit(1);
		}
		bond = (BOND *) malloc(sizeof(BOND) * cinfo.maxbond);
		if (bond == NULL) {
			fprintf(stderr, "memory allocation error for *bond\n");
			exit(1);
		}
		int i;
		for (i = 0; i < cinfo.maxbond; ++i) {
			bond[i].jflag = -1; /* bond type has not been assigned */
		}
		overflow_flag =
			rac(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
	}
	atomicnum(atomnum, atom);
	adjustatomname(atomnum, atom, 1);
	if(minfo.dcharge >= -9990) charge = minfo.dcharge;

	selectindex = (int *) malloc(sizeof(int) * atomnum);
	if (selectindex == NULL) {
		fprintf(stderr, "memory allocation error for *selectindex\n");
		exit(1);
	}
	equatomno = (int *) malloc(sizeof(int) * atomnum);
	if (equatomno == NULL) {
		fprintf(stderr, "memory allocation error for *equatomno\n");
		exit(1);
	}
	pathnum = (int *) malloc(sizeof(int) * atomnum);
	if (pathnum == NULL) {
		fprintf(stderr, "memory allocation error for *pathnum\n");
		exit(1);
	}
	selectelement = (int *) malloc(sizeof(int) * atomnum);
	if (selectelement == NULL) {
		fprintf(stderr, "memory allocation error for *selectelement\n");
		exit(1);
	}
	pathatomnum = (int *) malloc(sizeof(int) * atomnum);
	if (pathatomnum == NULL) {
		fprintf(stderr, "memory allocation error for *pathatomnum\n");
		exit(1);
	}
	if(atomnum > MAX_RESP_ATOM) {
		fprintf(stderr, "The number of atoms (%d) exceed the MAX_RESP_ATOM (%d) defined in respgen.c, extend MAX_RESP_ATOM and recompile the program\n", 
             atomnum, MAX_RESP_ATOM);
		exit(1);
	}
	for (i = 0; i < atomnum; i++) {
		pathatomnum[i] = MAXPATHATOMNUM;
		pathscore[i] = (double *) calloc(pathatomnum[i], sizeof(double));
		if (pathscore == NULL) {
			fprintf(stderr, "memory allocation error for *pathscore[%d]\n",
					i + 1);
			exit(1);
		}
	}
	equatom();

	icharge = (int *) malloc(sizeof(int) * atomnum);
       	if (icharge == NULL) {
               	fprintf(stderr, "memory allocation error for *icharge in respin()\n");
               	exit(1);
       	}
	pcharge = (double *) malloc(sizeof(double) * atomnum);
       	if (pcharge == NULL) {
               	fprintf(stderr, "memory allocation error for *pcharge in respin()\n");
               	exit(1);
       	}
	for(i=0;i<atomnum;i++) {
		pcharge[i] = 0.0;	
		icharge[i] = 0;
	}

	if(iaddinfo == 1) {
		if ((fpaddinfo = fopen(afilename, "r")) == NULL) {
        		fprintf(stderr, "Cannot open the additional file %s to read in main(), exit\n", afilename);
        		exit(1);
		}
		if ((fpcharge = fopen("QIN", "w")) == NULL) {
        		fprintf(stderr, "Cannot open file QIN, exit\n");
        		exit(1);
		}
		readinfo();
	}
	respin(method);
	if(iaddinfo == 1) {
		fclose(fpcharge);
		fclose(fpaddinfo);
	}
	printf("\n");
/*
	 free(atom);
	 free(selectindex);
	 free(selectelement);
	 free(equatomno);
	 free(pathnum);
	 free(pathatomnum);
	 for (i =0 ;i <atomnum; i++) free(pathscore[i]);
*/
	return (0);

}
