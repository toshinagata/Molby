char *amberhome;
# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <ctype.h>
# include <stdlib.h>
# include <string.h>
# include "common.h"
# include "define.h"
# include "atom.h"
# include "utility.c"
# include "common.c"
# include "ring.c"
# include "rotate.c"
# include "ac.c"
# include "charmm.c"
# include "mol2.c"
# include "mopcrt.c"
# include "divcrt.c"
# include "sqmcrt.c"
# include "mopint.c"
# include "mopout.c"
# include "divout.c"
# include "gcrt.c"
# include "gzmat.c"
# include "gout.c"
# include "pdb.c"
# include "csd.c"
# include "mdl.c"
# include "alc.c"
# include "hin.c"
# include "prep.c"
# include "rst.c"
# include "jzmat.c"
# include "jcrt.c"
# include "jout.c"
# include "charge.c"

int ao_flag = 0;
/*For addtional file 
 1, only readin coordinates
 2, only readin charge
 3, only readin atom names
 4, only readin atom types
 5, only readin bond types 
*/
int atomtype_flag = 0;				/*judge atom type? */
int bondtype_flag = 0;				/*judge bond type? */
int default_flag = 0;				/*assign default information? */
int atomname_flag = 0;				/*assign atom name? */
int atomicnum_flag = 0;				/*judge atomic number according to atom name ? */
int adjustatomname_flag = 0;		/*adjust atom name? */
int duplicatedname_flag = 0;		/*check atom name duplication? */
int cartcoord_flag = 0;				/*generate coordinate from internal coordinate ? */
int connect_flag = 0;				/*judge atom connectivity and generate bond ? */
int divcon_flag = 2;  
int ek_flag = 0; 				/* read empircal calculation keyword or not */

int atomnum = 0;
int bondnum = 0;
int ringnum = 0;
ATOM *atom;
BOND *bond;
RING *ring;
AROM *arom;
MOLINFO minfo;
CONTROLINFO cinfo;

int atomnum_tmp = 0;
int bondnum_tmp = 0;
ATOM *atom_tmp;
BOND *bond_tmp;
char line[MAXCHAR];
char ifilename[MAXCHAR];
char ofilename[MAXCHAR];
char afilename[MAXCHAR] = "";
char cfilename[MAXCHAR];

/*The following four functions, read_at(), write_at(), read_bt() and write_bt() are
used in amber sybyl interface development and are unreachable to the users
*/
int ra_flag = 0;
int rb_flag = 0;
int wa_flag = 0;
int wb_flag = 0;
char at_filename[MAXCHAR];
char bt_filename[MAXCHAR];

FILE *fpin;
FILE *fpout;

size_t copied_size;

void usage()
{
	if (strcmp(COLORTEXT, "YES") == 0 || strcmp(COLORTEXT, "yes") == 0) {
		printf("[31mUsage: antechamber -i  [0m input file name\n"
			   "[31m                   -fi [0m input file format\n"
			   "[31m                   -o  [0m output file name\n"
			   "[31m                   -fo [0m output file format\n"
			   "[31m                   -c  [0m charge method\n"
			   "[31m                   -cf [0m charge file name\n"
			   "[31m                   -nc [0m net molecular charge (int)\n"
			   "[31m                   -a  [0m additional file name\n"
			   "[31m                   -fa [0m additional file format\n"
			   "[31m                   -ao [0m additional file operation\n"
			   "[35m                        crd [0m: only read in coordinate\n"
			   "[35m                        crg[0m: only read in charge\n"
			   "[35m                        name  [0m: only read in atom name\n"
			   "[35m                        type  [0m: only read in atom type\n"
			   "[35m                        bond  [0m: only read in bond type \n"
			   "[31m                   -m  [0m multiplicity (2S+1), default is 1\n"
			   "[31m                   -rn [0m residue name, overrides input file, default is MOL\n"
			   "[31m                   -rf [0m residue toplogy file name in prep input file,\n"
               "                                   default is molecule.res\n"
			   "[31m                   -ch [0m check file name for gaussian, default is molecule\n"
			   "[31m                   -ek [0m mopac or sqm keyword, inside quotes\n"
			   "[31m                   -gk [0m gaussian keyword, inside quotes\n"
			   "[31m                   -df [0m am1-bcc flag, 2 - use sqm(default); 0 - use mopac\n"
			   "[31m                   -at [0m atom type, can be gaff (default), amber, bcc and sybyl\n"
			   "[31m                   -du [0m fix duplicate atom names: yes(y)[default] or no(n)\n"
			   "[31m                   -j  [0m atom type and bond type prediction index, default is 4 \n"
			   "[35m                        0    [0m: no assignment\n"
			   "[35m                        1    [0m: atom type \n"
			   "[35m                        2    [0m: full  bond types \n"
			   "[35m                        3    [0m: part  bond types \n"
			   "[35m                        4    [0m: atom and full bond type \n"
			   "[35m                        5    [0m: atom and part bond type \n"
			   "[31m                   -s  [0m status information: 0(brief), 1(default) or 2(verbose)\n"
			   "[31m                   -pf [0m remove intermediate files: yes(y) or no(n)[default]\n"
			   "                   -i -o -fi and -fo must appear; others are optional[0m");
		printf
			("\n\n	         	    [31m List of the File Formats [0m \n");
		printf
			("\n	 	file format type  abbre. index | file format type abbre. index");
		printf
			("\n		--------------------------------------------------------------- ");
		printf
			("\n		Antechamber        ac       1  | Sybyl Mol2         mol2    2 ");
                printf
                        ("\n		PDB                pdb      3  | Modified PDB       mpdb    4 ");
		printf
			("\n		AMBER PREP (int)   prepi    5  | AMBER PREP (car)   prepc   6 ");
		printf
			("\n		Gaussian Z-Matrix  gzmat    7  | Gaussian Cartesian gcrt    8 ");
		printf
			("\n		Mopac Internal     mopint   9  | Mopac Cartesian    mopcrt 10 ");
		printf
			("\n		Gaussian Output    gout    11  | Mopac Output       mopout 12 ");
		printf
			("\n		Alchemy            alc     13  | CSD                csd    14 ");
		printf
			("\n		MDL                mdl     15  | Hyper              hin    16 ");
		printf
			("\n		AMBER Restart      rst     17  | Jaguar Cartesian   jcrt   18 ");
		printf
			("\n		Jaguar Z-Matrix    jzmat   19  | Jaguar Output      jout   20");
		printf
			("\n		Divcon Input       divcrt  21  | Divcon Output      divout 22");
		printf
			("\n		Charmm             charmm  23  | SQM Output         sqmout 24");
		printf
			("\n		--------------------------------------------------------------\n");
		printf
			("\n                AMBER restart file can only be read in as additional file.\n");
		printf
			("\n	         	    [31m List of the Charge Methods [0m \n");
		printf
			("\n		charge method     abbre.  index | charge method    abbre. index");
		printf
			("\n		----------------------------------------------------------------  ");
		printf
			("\n		RESP               resp     1  |  AM1-BCC            bcc     2");
		printf
			("\n		CM1                cm1      3  |  CM2                cm2     4");
		printf
			("\n		ESP (Kollman)      esp      5  |  Mulliken           mul     6");
		printf
			("\n		Gasteiger          gas      7  |  Read in charge     rc      8");
		printf
			("\n		Write out charge   wc       9  |  Delete Charge      dc     10");
		printf
			("\n		----------------------------------------------------------------\n");
	} else {
	}
		printf("Usage: antechamber -i  input file name\n"
			   "                   -fi input file format\n"
			   "                   -o  output file name\n"
			   "                   -fo output file format\n"
			   "                   -c  charge method\n"
			   "                   -cf charge file name\n"
			   "                   -nc net molecular charge (int)\n"
			   "                   -a  additional file name\n"
			   "                   -fa additional file format\n"
			   "                   -ao additional file operation\n"
			   "                        crd  only read in coordinate\n"
			   "                        crg only read in charge\n"
			   "                        name   only read in atom name\n"
			   "                        type   only read in atom type\n"
			   "                        bond   only read in bond type \n"
			   "                   -m  multiplicity (2S+1), default is 1\n"
			   "                   -rn residue name, overrides input file, default is MOL\n"
			   "                   -rf residue toplogy file name in prep input file,\n"
               "                              default is molecule.res\n"
			   "                   -ch check file name for gaussian, default is molecule\n"
			   "                   -ek mopac or sqm keyword, inside quotes\n"
			   "                   -gk gaussian keyword, inside quotes\n"
			   "                   -df am1-bcc flag, 2 - use sqm(default); 0 - use mopac\n"
			   "                   -at atom type, can be gaff (default), amber, bcc and sybyl\n"
			   "                   -du fix duplicate atom names: yes(y)[default] or no(n)\n"
			   "                   -j  atom type and bond type prediction index, default is 4 \n"
			   "                        0     no assignment\n"
			   "                        1     atom type \n"
			   "                        2     full  bond types \n"
			   "                        3     part  bond types \n"
			   "                        4     atom and full bond type \n"
			   "                        5     atom and part bond type \n"
			   "                   -s  status information: 0(brief), 1(default) or 2(verbose)\n"
			   "                   -pf remove intermediate files: yes(y) or no(n)[default]\n"
			   "                   -i -o -fi and -fo must appear; others are optional");
		printf
			("\n\n	         	     List of the File Formats \n");
		printf
			("\n	 	file format type  abbre. index | file format type abbre. index");
		printf
			("\n		--------------------------------------------------------------- ");
		printf
			("\n		Antechamber        ac       1  | Sybyl Mol2         mol2    2 ");
                printf
                        ("\n		PDB                pdb      3  | Modified PDB       mpdb    4 ");
		printf
			("\n		AMBER PREP (int)   prepi    5  | AMBER PREP (car)   prepc   6 ");
		printf
			("\n		Gaussian Z-Matrix  gzmat    7  | Gaussian Cartesian gcrt    8 ");
		printf
			("\n		Mopac Internal     mopint   9  | Mopac Cartesian    mopcrt 10 ");
		printf
			("\n		Gaussian Output    gout    11  | Mopac Output       mopout 12 ");
		printf
			("\n		Alchemy            alc     13  | CSD                csd    14 ");
		printf
			("\n		MDL                mdl     15  | Hyper              hin    16 ");
		printf
			("\n		AMBER Restart      rst     17  | Jaguar Cartesian   jcrt   18 ");
		printf
			("\n		Jaguar Z-Matrix    jzmat   19  | Jaguar Output      jout   20");
		printf
			("\n		Divcon Input       divcrt  21  | Divcon Output      divout 22");
		printf
			("\n		Charmm             charmm  23  | SQM Output         sqmout 24");
		printf
			("\n		--------------------------------------------------------------\n");
		printf
			("\n                AMBER restart file can only be read in as additional file.\n");
		printf
			("\n	         	     List of the Charge Methods \n");
		printf
			("\n		charge method     abbre.  index | charge method    abbre. index");
		printf
			("\n		----------------------------------------------------------------  ");
		printf
			("\n		RESP               resp     1  |  AM1-BCC            bcc     2");
		printf
			("\n		CM1                cm1      3  |  CM2                cm2     4");
		printf
			("\n		ESP (Kollman)      esp      5  |  Mulliken           mul     6");
		printf
			("\n		Gasteiger          gas      7  |  Read in charge     rc      8");
		printf
			("\n		Write out charge   wc       9  |  Delete Charge      dc     10");
		printf
			("\n		----------------------------------------------------------------\n");
}

void memory(int flag, int maxatom, int maxbond, int maxring)
{
	if (flag == 0) {
		atom = (ATOM *) malloc(sizeof(ATOM) * maxatom);
		if (atom == NULL) {
			fprintf(stderr, "memory allocation error for *atom\n");
			exit(1);
		}
		arom = (AROM *) malloc(sizeof(AROM) * maxatom);
		if (arom == NULL) {
			fprintf(stderr, "memory allocation error for *arom\n");
			exit(1);
		}
		bond = (BOND *) malloc(sizeof(BOND) * maxbond);
		if (bond == NULL) {
			fprintf(stderr, "memory allocation error for *bond\n");
			exit(1);
		}
		int i;
		for (i = 0; i < maxbond; ++i) {
			bond[i].jflag = -1; /* bond type has not been assigned */
			bond[i].type = 0; /* 20110707 Toshi Nagata */
		}
	}
/*flag = 1  <->atom
       = 2  <->bond
       = 3  <->arom
       = 4  <->atom + bond
       = 5  <->atom + arom 
       = 6  <->bond + arom
       = 7  <->atom + arom +bond
*/
	if (flag == 1 || flag == 4 || flag == 5 || flag == 7) {
		free(atom);
		atom = (ATOM *) malloc(sizeof(ATOM) * maxatom);
		if (atom == NULL) {
			fprintf(stderr, "memory allocation error for *atom\n");
			exit(1);
		}
	}
	if (flag == 2 || flag == 4 || flag == 6 || flag == 7) {
		free(bond);
		bond = (BOND *) malloc(sizeof(BOND) * maxbond);
		if (bond == NULL) {
			fprintf(stderr, "memory allocation error for *bond\n");
			exit(1);
		}
		int i;
		for (i = 0; i < maxbond; ++i) {
			bond[i].jflag = -1; /* bond type has not been assigned */
			bond[i].type = 0; /* 20110707 Toshi Nagata */
		}
	}
	if (flag == 3 || flag == 5 || flag == 6 || flag == 7) {
		free(arom);
		arom = (AROM *) malloc(sizeof(AROM) * maxatom);
		if (arom == NULL) {
			fprintf(stderr, "memory allocation error for *arom\n");
			exit(1);
		}
	}
	if (flag == 8) {
		free(ring);
		ring = (RING *) malloc(sizeof(RING) * maxring);
		if (ring == NULL) {
			fprintf(stderr, "memory allocation error for *ring\n");
			exit(1);
		}
	}
}

void judgebondtype(int atomnum, ATOM * atom, int bondnum, BOND * bond,
				   CONTROLINFO cinfo, MOLINFO minfo, int bondtype_flag)
{
	char tmpchar[MAXCHAR];
	char *options;
	int status = 0;
	wac("ANTECHAMBER_BOND_TYPE.AC0", atomnum, atom, bondnum, bond, cinfo,
		minfo);
	copied_size = build_exe_path(tmpchar, "bondtype", MAXCHAR );
        if (bondtype_flag == 1)
                options = " -j part -i ANTECHAMBER_BOND_TYPE.AC0"
                          " -o ANTECHAMBER_BOND_TYPE.AC -f ac" ;
        else
                options = " -j full -i ANTECHAMBER_BOND_TYPE.AC0"
                          " -o ANTECHAMBER_BOND_TYPE.AC -f ac" ;
        strncat(tmpchar, options, MAXCHAR - copied_size );

	if (cinfo.intstatus == 2)
		fprintf(stdout, "Running: %s\n", tmpchar);
	status = system(tmpchar);
        if(status != 0) {
                fprintf(stderr, "Error: cannot run \"%s\" in judgebondtype() of antechamber.c properly, exit\n", tmpchar);
                exit(1);
        }
	rac("ANTECHAMBER_BOND_TYPE.AC", &atomnum, atom, &bondnum, bond, &cinfo,
		&minfo);
}

int read_at(char *filename)  {
	FILE *fpin;
	int num = 0;
	char tmpchar1[MAXCHAR];
	char tmpchar2[MAXCHAR];
	char tmpchar3[MAXCHAR];
	char line[MAXCHAR];

        if ((fpin = fopen(filename, "r")) == NULL) {
                fprintf(stderr, "Cannot open file %s in read_at(), exit\n", filename);
                return 0;
        }
        for (;;) {
                if (fgets(line, MAXCHAR, fpin) == NULL)
                        break;
                        sscanf(line, "%s%s%s", tmpchar1, tmpchar2, tmpchar3);
                        strcpy(atom[num++].ambername, tmpchar3);
        }
	fclose(fpin);
	return 0;
}

int read_bt(char *filename)  {
	FILE *fpin;
	int num = 0;
	char tmpchar1[MAXCHAR];
	char tmpchar2[MAXCHAR];
	char tmpchar3[MAXCHAR];
	char line[MAXCHAR];
        if ((fpin = fopen(filename, "r")) == NULL) {
                fprintf(stderr, "Cannot open file %s in read_bt(), exit\n", filename);
                return 0;
        }
        for (;;) {
                if (fgets(line, MAXCHAR, fpin) == NULL)
                        break;
                        sscanf(line, "%s%s%s%d", tmpchar1, tmpchar2, tmpchar3, &bond[num++].type);
        }
	fclose(fpin);
	return 0;

}

int write_at(char *filename)  {
	FILE *fpout;
	int i;
        if ((fpout = fopen(filename, "w")) == NULL) {
                fprintf(stderr, "Cannot open file %s to write in write_at(), exit\n", filename);
                return 0;
        }
	for(i=0;i<atomnum;i++)
		fprintf(fpout, "%5d %5s %5s\n", i+1, atom[i].name, atom[i].ambername);
	fclose(fpout);
	return 0;

}

int write_bt(char *filename)  {
	FILE *fpout;
	int i;
        if ((fpout = fopen(filename, "w")) == NULL) {
                fprintf(stderr, "Cannot open file %s to write in write_bt(), exit\n", filename);
                return 0;
        }
	for(i=0;i<bondnum;i++)
		fprintf(fpout, "%5d %5d %5d %5d\n", i+1, bond[i].bondi+1, bond[i].bondj, bond[i].type);
	fclose(fpout);
	return 0;

}

int main(int argc, char *argv[])
{
	int i,j,k;
	int index;
	int status = 0;
	char tmpchar[MAXCHAR];
	double fraction;
	double tmpf;
	int overflow_flag = 0;			/*if overflow_flag ==1, reallocate memory */

    amberhome = (char *) getenv("AMBERHOME");
    if( amberhome == NULL ){
       fprintf( stderr, "AMBERHOME is not set!\n" );
       exit(1);
    }
	if (argc == 2)
		if (strncmp(argv[1], "-h", 2) == 0
			|| strncmp(argv[1], "-H", 2) == 0) {
			usage();
			exit(1);
		}
	if (argc == 1) {
		usage();
		exit(1);
	}

/* 	set defaults information */
	default_minfo(&minfo);
	default_cinfo(&cinfo);
	atomtype_flag = 0;
	bondtype_flag = 0;
	default_flag = 0;
	atomname_flag = 0;
	atomicnum_flag = 0;
	adjustatomname_flag = 0;
	duplicatedname_flag = 1;
	cartcoord_flag = 0;
	connect_flag = 0;


	index = 0;
	for (i = 1; i < argc - 1; i += 2) {
		if (strcmp(argv[i], "-i") == 0) {
			index++;
			strcpy(ifilename, argv[i + 1]);
			continue;
		} else if (strcmp(argv[i], "-o") == 0) {
			index++;
			strcpy(ofilename, argv[i + 1]);
			continue;
		} else if (strcmp(argv[i], "-fi") == 0) {
			index++;
			strcpy(cinfo.intype, argv[i + 1]);
			continue;
		} else if (strcmp(argv[i], "-j") == 0) {
			cinfo.prediction_index = atoi(argv[i + 1]);
			continue;
		} else if (strcmp(argv[i], "-fo") == 0) {
			index++;
			strcpy(cinfo.outtype, argv[i + 1]);
			continue;
		} else if (strcmp(argv[i], "-a") == 0) {
			strcpy(afilename, argv[i + 1]);
			continue;
		} else if (strcmp(argv[i], "-fa") == 0) {
			strcpy(cinfo.atype, argv[i + 1]);
			continue;
		} else if (strcmp(argv[i], "-c") == 0) {
			strcpy(cinfo.chargetype, argv[i + 1]);
			continue;
		} else if (strcmp(argv[i], "-m") == 0) {
			minfo.multiplicity = atoi(argv[i + 1]);
			continue;
		} else if (strcmp(argv[i], "-s") == 0) {
			cinfo.intstatus = atoi(argv[i + 1]);
			continue;
		} else if (strcmp(argv[i], "-ra") == 0) {
			strcpy(at_filename, argv[i + 1]);
			ra_flag = 1;
			continue;
		} else if (strcmp(argv[i], "-rb") == 0) {
			strcpy(bt_filename, argv[i + 1]);
			rb_flag = 1;
			continue;
		} else if (strcmp(argv[i], "-wa") == 0) {
			strcpy(at_filename, argv[i + 1]);
			wa_flag = 1;
			continue;
		} else if (strcmp(argv[i], "-wb") == 0) {
			strcpy(bt_filename, argv[i + 1]);
			wb_flag = 1;
			continue;
		} else if (strcmp(argv[i], "-ao") == 0) {
			if (strcmp(argv[i + 1], "crd") == 0
				|| strcmp(argv[i + 1], "CRD") == 0) {
				ao_flag = 1;
				continue;
			}
			if (strcmp(argv[i + 1], "crg") == 0
				|| strcmp(argv[i + 1], "CRG") == 0) {
				ao_flag = 2;
				continue;
			}
			if (strcmp(argv[i + 1], "name") == 0
				|| strcmp(argv[i + 1], "NAME") == 0) {
				ao_flag = 3;
				continue;
			}
			if (strcmp(argv[i + 1], "type") == 0
				|| strcmp(argv[i + 1], "NAME") == 0) {
				ao_flag = 4;
				continue;
			}
			if (strcmp(argv[i + 1], "bond") == 0
				|| strcmp(argv[i + 1], "BOND") == 0) {
				ao_flag = 5;
				continue;
			}
		} else if (strcmp(argv[i], "-cf") == 0) {
			strcpy(cfilename, argv[i + 1]);
			continue;
		} else if (strcmp(argv[i], "-nc") == 0) {
			minfo.usercharge = atoi(argv[i + 1]);
			continue;
		} else if (strcmp(argv[i], "-rn") == 0) {
			strcpy(minfo.longresname, argv[i + 1]);
/*			strncpy(minfo.resname, argv[i + 1], 3); */
			strcpy(minfo.resname, argv[i + 1]); 
			cinfo.rnindex = 1;
			continue;
		} else if (strcmp(argv[i], "-ch") == 0) {
			strcpy(minfo.chkfile, argv[i + 1]);
			continue;
		} else if (strcmp(argv[i], "-at") == 0) {
			strcpy(minfo.atom_type_def, argv[i + 1]);
			continue;
		} else if (strcmp(argv[i], "-du") == 0) {
			if (strcmp(argv[i + 1], "Yes") == 0
				|| strcmp(argv[i + 1], "Y") == 0
				|| strcmp(argv[i + 1], "yes") == 0
				|| strcmp(argv[i + 1], "y") == 0) 
				duplicatedname_flag = 1;
			if (strcmp(argv[i + 1], "NO") == 0
				|| strcmp(argv[i + 1], "No") == 0
				|| strcmp(argv[i + 1], "N") == 0
				|| strcmp(argv[i + 1], "no") == 0
				|| strcmp(argv[i + 1], "n") == 0) 
				duplicatedname_flag = 0;
			continue;
		} else if (strcmp(argv[i], "-rf") == 0) {
			strcpy(minfo.resfilename, argv[i + 1]);
			continue;
		} else if (strcmp(argv[i], "-pf") == 0) {
			if (strcmp(argv[i + 1], "yes") == 0
				|| strcmp(argv[i + 1], "YES") == 0
				|| strcmp(argv[i + 1], "Y") == 0
				|| strcmp(argv[i + 1], "y") == 0
				|| strcmp(argv[i + 1], "Yes") == 0)
				cinfo.pfindex = 1;
			if (strcmp(argv[i + 1], "no") == 0
				|| strcmp(argv[i + 1], "NO") == 0
				|| strcmp(argv[i + 1], "N") == 0
				|| strcmp(argv[i + 1], "n") == 0
				|| strcmp(argv[i + 1], "No") == 0)
				cinfo.pfindex = 0;
			continue;
		} else if (strcmp(argv[i], "-ek") == 0) {
			strcpy(minfo.ekeyword, argv[i + 1]);
			ek_flag = 1;
			continue;
		} else if (strcmp(argv[i], "-gk") == 0) {
			strcpy(minfo.gkeyword, argv[i + 1]);
			continue;
		} else if (strcmp(argv[i], "-df") == 0) {
			divcon_flag = atoi(argv[i + 1]);
			if(divcon_flag != 0 && divcon_flag != 2) {
				fprintf(stderr, "divcon_flag can only be 2 (sqm) or 0 (mopac), exit\n");
				exit(1);
			}
			minfo.divcon = divcon_flag;
			continue;
		} else {
			fprintf(stderr, "Flag not recognized: %s\n", argv[i]);
			fprintf(stderr,
					"Use antechamber -h for command-line syntax\n");
			exit(1);
		}
	}

	if (index != 4 && (strcmp(cinfo.chargetype, "wc") != 0)) {
		fprintf(stderr, "Need both input and output files & formats\n");
		fprintf(stderr, "Use antechamber -h for command-line syntax\n");
		exit(1);
	}

/* 	set ekeyword if it is not read in*/
	if(ek_flag == 0) {
		if(divcon_flag == 0)
			strcpy(minfo.ekeyword, minfo.mkeyword);		
		else if(divcon_flag == 1) 
			strcpy(minfo.ekeyword, minfo.dkeyword);		
		else if(divcon_flag == 2) 
			strcpy(minfo.ekeyword, minfo.skeyword);		
	}

        create_output_file_comment( argc, argv );
/* 	for connect.tpl and radius parameter files */
        minfo.connect_file[0] = '\0';
        strcpy(minfo.connect_file, amberhome);
        strcat(minfo.connect_file, "/dat/antechamber/CONNECT.TPL");
        minfo.radius_file[0] = '\0';
        strcpy(minfo.radius_file, amberhome);
        strcat(minfo.radius_file, "/dat/antechamber/RADIUS.DAT");
/*      allocate memory using default parameters MAXATOM and MAXBOND */
	memory(0, MAXATOM, MAXBOND, MAXRING);

/******************************************/
/* 	The following codes readin input file */
/******************************************/

	if (strcmp("ac", cinfo.intype) == 0 || strcmp("1", cinfo.intype) == 0) {
		overflow_flag =
			rac(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag =
				rac(ifilename, &atomnum, atom, &bondnum, bond, &cinfo,
					&minfo);
		}
		adjustatomname_flag = 1;
		atomicnum_flag = 1;
	}

	if (strcmp("mol2", cinfo.intype) == 0
		|| strcmp("2", cinfo.intype) == 0) {
		
		overflow_flag =
			rmol2(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo, 0);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag =
				rmol2(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo, 0);
		}
		default_flag = 1;
		atomicnum_flag = 1;
		adjustatomname_flag = 1;
	}
	if (strcmp("mopint", cinfo.intype) == 0
		|| strcmp("9", cinfo.intype) == 0) {
		overflow_flag = rmopint(ifilename, &atomnum, atom, cinfo, minfo);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag =
				rmopint(ifilename, &atomnum, atom, cinfo, minfo);
		}
		cartcoord_flag = 1;
		atomicnum_flag = 1;
		atomname_flag = 1;
		default_flag = 1;
		connect_flag = 1;
		bondtype_flag = 2;

	}

	if (strcmp("mopcrt", cinfo.intype) == 0
		|| strcmp("10", cinfo.intype) == 0) {
		overflow_flag = rmopcrt(ifilename, &atomnum, atom, cinfo, minfo);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag =
				rmopcrt(ifilename, &atomnum, atom, cinfo, minfo);
		}
		atomicnum_flag = 1;
		atomname_flag = 1;
		default_flag = 1;
		connect_flag = 1;
		bondtype_flag = 2;
	}

	if (strcmp("mopout", cinfo.intype) == 0
		|| strcmp("12", cinfo.intype) == 0) {
		overflow_flag = rmopout(ifilename, &atomnum, atom, &cinfo, &minfo);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag =
				rmopout(ifilename, &atomnum, atom, &cinfo, &minfo);
		}
/* when read in a mopout file, always output a pdb file that has the optimized coordinates */
                atom_tmp = (ATOM *) malloc(sizeof(ATOM) * cinfo.maxatom);
                if (atom_tmp == NULL) {
                        fprintf(stderr, "memory allocation error for *atom_tmp\n");
                        exit(1);
                }
		for(i=0;i<atomnum;i++) {
			atom_tmp[i] = atom[i];
			atom_tmp[i].connum = 0;
		}
		rmopout_coord(ifilename, atom_tmp);
		wpdb("mopac.pdb", atomnum, atom_tmp);
		free(atom_tmp);
		atomicnum_flag = 1;
		atomname_flag = 1;
		default_flag = 1;
		connect_flag = 1;
		bondtype_flag = 2;
	}

	if (strcmp("gcrt", cinfo.intype) == 0
		|| strcmp("8", cinfo.intype) == 0) {
		overflow_flag = rgcrt(ifilename, &atomnum, atom, cinfo, minfo);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag = rgcrt(ifilename, &atomnum, atom, cinfo, minfo);
		}
		atomname_flag = 1;
		default_flag = 1;
		connect_flag = 1;
		bondtype_flag = 2;
	}
	if (strcmp("gzmat", cinfo.intype) == 0
		|| strcmp("7", cinfo.intype) == 0) {
		overflow_flag = rgzmat(ifilename, &atomnum, atom, cinfo, minfo);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag =
				rgzmat(ifilename, &atomnum, atom, cinfo, minfo);
		}
		cartcoord_flag = 1;
		atomicnum_flag = 1;
		atomname_flag = 1;
		default_flag = 1;
		connect_flag = 1;
		bondtype_flag = 2;
	}

	if (strcmp("gout", cinfo.intype) == 0
		|| strcmp("11", cinfo.intype) == 0) {
		overflow_flag = rgout(ifilename, &atomnum, atom, cinfo, &minfo);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag = rgout(ifilename, &atomnum, atom, cinfo, &minfo);
		}
		atomname_flag = 1;
		default_flag = 1;
		connect_flag = 1;
		bondtype_flag = 2;
	}

	if (strcmp("jcrt", cinfo.intype) == 0
		|| strcmp("18", cinfo.intype) == 0) {
		overflow_flag = rjcrt(ifilename, &atomnum, atom, cinfo, minfo);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag = rjcrt(ifilename, &atomnum, atom, cinfo, minfo);
		}
		default_flag = 1;
		connect_flag = 1;
		bondtype_flag = 2;
	}

        if (strcmp("jzmat", cinfo.intype) == 0
                || strcmp("19", cinfo.intype) == 0) {
                overflow_flag = rjzmat(ifilename, &atomnum, atom, cinfo, minfo);
                if (overflow_flag) {
                        cinfo.maxatom = atomnum + 10;
                        cinfo.maxbond = bondnum + 10;
                        memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
                        overflow_flag =
                                rjzmat(ifilename, &atomnum, atom, cinfo, minfo);
                }
                cartcoord_flag = 1;
                atomicnum_flag = 1;
                atomname_flag = 0;
                default_flag = 1;
                connect_flag = 1;
                bondtype_flag = 2;
        }

        if (strcmp("jout", cinfo.intype) == 0
                || strcmp("20", cinfo.intype) == 0) {
                overflow_flag = rjout(ifilename, &atomnum, atom, cinfo, &minfo);
                if (overflow_flag) {
                        cinfo.maxatom = atomnum + 10;
                        cinfo.maxbond = bondnum + 10;
                        memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
                        overflow_flag = rjout(ifilename, &atomnum, atom, cinfo, &minfo);
                }
                default_flag = 1;
                connect_flag = 1;
                bondtype_flag = 2;
        }

	if (strcmp("pdb", cinfo.intype) == 0 || strcmp("3", cinfo.intype) == 0) {
		overflow_flag = rpdb(ifilename, &atomnum, atom, cinfo, minfo, 0);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag =
				rpdb(ifilename, &atomnum, atom, cinfo, minfo, 0);
		}
		adjustatomname_flag = 1;
		atomicnum_flag = 1;
		default_flag = 1;
		bondtype_flag = 2;
		index = 0;
		for(i=0; i<atomnum; i++)
			if(atom[i].connum > 0) {
				index = 1;
				break;
			}
		if(index == 1) {
			bondnum = 0;
			for(i=0; i<atomnum-1; i++)
				for(j=i+1; j<atomnum; j++) 
					for(k=0; k<6; k++)  {
						if(atom[i].con[k] == -1) break;
						if(atom[i].con[k] == j) {
							bond[bondnum].bondi = i;
							bond[bondnum].bondj = j;
							bondnum++;
		  				}
					}
		}
		else 
			connect_flag = 1;
	}

	if (strcmp("mpdb", cinfo.intype) == 0
		|| strcmp("4", cinfo.intype) == 0) {
		overflow_flag = rpdb(ifilename, &atomnum, atom, cinfo, minfo, 1);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag =
				rpdb(ifilename, &atomnum, atom, cinfo, minfo, 1);
		}
		adjustatomname_flag = 1;
		atomicnum_flag = 1;
		default_flag = 0;
		bondtype_flag = 2;
		index = 0;
		for(i=0; i<atomnum; i++)
			if(atom[i].connum > 0) {
				index = 1;
				break;
			}
                if(index == 1) {
                        bondnum = 0;
                        for(i=0; i<atomnum-1; i++)
                                for(j=i+1; j<atomnum; j++) 
                                        for(k=0; k<6; k++)  {       
                                                if(atom[i].con[k] == -1) break;
                                                if(atom[i].con[k] == j) {
                                                        bond[bondnum].bondi = i;
                                                        bond[bondnum].bondj = j;
                                                        bondnum++;
                                                }
					}
                }
                else 
                        connect_flag = 1;

	}
	if (strcmp("csd", cinfo.intype) == 0
		|| strcmp("14", cinfo.intype) == 0) {
		overflow_flag = rcsd(ifilename, &atomnum, atom, cinfo, minfo);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag = rcsd(ifilename, &atomnum, atom, cinfo, minfo);
		}
		atomicnum_flag = 1;
		default_flag = 1;
		connect_flag = 1;
		bondtype_flag = 2;
	}

	if (strcmp("mdl", cinfo.intype) == 0
		|| strcmp("sd", cinfo.intype) == 0
		|| strcmp("sdf", cinfo.intype) == 0
		|| strcmp("15", cinfo.intype) == 0) {
		overflow_flag =
			rmdl(ifilename, &atomnum, atom, &bondnum, bond, cinfo, minfo);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag =
				rmdl(ifilename, &atomnum, atom, &bondnum, bond, cinfo,
					 minfo);
		}
		atomicnum_flag = 1;
		atomname_flag = 1;
		default_flag = 1;
		bondtype_flag = 1;
	}

	if (strcmp("alc", cinfo.intype) == 0
		|| strcmp("13", cinfo.intype) == 0) {
		overflow_flag =
			ralc(ifilename, &atomnum, atom, &bondnum, bond, cinfo, minfo);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag =
				ralc(ifilename, &atomnum, atom, &bondnum, bond, cinfo,
					 minfo);
		}
		atomicnum_flag = 1;
		atomname_flag = 1;
		bondtype_flag = 1;
		default_flag = 1;
	}

	if (strcmp("hin", cinfo.intype) == 0
		|| strcmp("16", cinfo.intype) == 0) {
		overflow_flag =
			rhin(ifilename, &atomnum, atom, &bondnum, bond, cinfo, minfo);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag =
				rhin(ifilename, &atomnum, atom, &bondnum, bond, cinfo,
					 minfo);
		}
		atomicnum_flag = 1;
		atomname_flag = 1;
		bondtype_flag = 1;
		default_flag = 1;
	}

	if (strcmp("prepi", cinfo.intype) == 0
		|| strcmp("5", cinfo.intype) == 0) {
		overflow_flag = rprepi(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag =
				rprepi(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
		}
		atomicnum_flag = 1;
		bondtype_flag = 2;
		connect_flag = 0;
	}

	if (strcmp("prepc", cinfo.intype) == 0
		|| strcmp("5", cinfo.intype) == 0) {
		overflow_flag = rprepc(ifilename, &atomnum, atom, &cinfo, &minfo);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag =
				rprepc(ifilename, &atomnum, atom, &cinfo, &minfo);
		}
		atomicnum_flag = 1;
		bondtype_flag = 2;
		connect_flag = 1;
	}

	if (strcmp("rst", cinfo.intype) == 0
		|| strcmp("17", cinfo.intype) == 0) {
		fprintf(stderr,
				"RST (17) file format can only be additional file because it only has coordinate information\n");
		exit(1);
	}

        if (strcmp("divcrt", cinfo.intype) == 0
                || strcmp("21", cinfo.intype) == 0) {
                overflow_flag = rdivcrt(ifilename, &atomnum, atom, cinfo, minfo);
                if (overflow_flag) {
                        cinfo.maxatom = atomnum + 10;
                        cinfo.maxbond = bondnum + 10;
                        memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
                        overflow_flag =
                                rdivcrt(ifilename, &atomnum, atom, cinfo, minfo);
                }
                atomicnum_flag = 1;
                atomname_flag = 1;
                default_flag = 1;
                connect_flag = 1;
                bondtype_flag = 2;
        }
                                                                                                                                                                                                           
        if (strcmp("divout", cinfo.intype) == 0
                || strcmp("22", cinfo.intype) == 0) {
                overflow_flag = rdivout(ifilename, &atomnum, atom, &cinfo, &minfo);
                if (overflow_flag) {
                        cinfo.maxatom = atomnum + 10;
                        cinfo.maxbond = bondnum + 10;
                        memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
                        overflow_flag =
                                rdivout(ifilename, &atomnum, atom, &cinfo, &minfo);
                }
/* when read in a divout file, always output a pdb file that has the optimized coordinates */
                atom_tmp = (ATOM *) malloc(sizeof(ATOM) * cinfo.maxatom);
                if (atom_tmp == NULL) {
                        fprintf(stderr, "memory allocation error for *atom_tmp\n");
                        exit(1);
                }
		for(i=0;i<atomnum;i++) {
			atom_tmp[i] = atom[i];
			atom_tmp[i].connum = 0;
		}
		rdivout_coord(ifilename, atom_tmp);
		wpdb("divcon.pdb", atomnum, atom_tmp);
		free(atom_tmp);
                atomicnum_flag = 1;
                atomname_flag = 1;
                default_flag = 1;
                connect_flag = 1;
                bondtype_flag = 2;
        }
	if (strcmp("charmm", cinfo.intype) == 0
		|| strcmp("23", cinfo.intype) == 0) {
		
		overflow_flag =
			rcharmm(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag =
				rcharmm(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
		}
		default_flag = 1;
                atomicnum_flag = 1;
		adjustatomname_flag = 1;
	}


/*****************************************************************************/
/*	assign atomtype_flag, bondtype_flag and charge_flag according to -j flag */
/*****************************************************************************/

	if (cinfo.prediction_index == 0) {
		atomtype_flag = 0;
		bondtype_flag = 0;
	}
	if (cinfo.prediction_index == 1) {
		atomtype_flag = 1;
		bondtype_flag = 0;
	}
	if (cinfo.prediction_index == 2) {
		atomtype_flag = 0;
		bondtype_flag = 2;
	}
	if (cinfo.prediction_index == 3) {
		atomtype_flag = 0;
		bondtype_flag = 1;
	}
	if (cinfo.prediction_index == 4) {
		atomtype_flag = 1;
		bondtype_flag = 2;
	}
	if (cinfo.prediction_index == 5) {
		atomtype_flag = 1;
		bondtype_flag = 1;
	}

/*	reassign the connect_flag according to the output types 	*/
        if (strcmp("mopcrd", cinfo.outtype) == 0
                || strcmp("mopout", cinfo.outtype) == 0 
                || strcmp("gcrt", cinfo.outtype) == 0 
                || strcmp("gout", cinfo.outtype) == 0 
                || strcmp("jcrt", cinfo.outtype) == 0 
                || strcmp("jout", cinfo.outtype) == 0 
                || strcmp("pdb", cinfo.outtype) == 0 
                || strcmp("rst", cinfo.outtype) == 0 
                || strcmp("3", cinfo.outtype) == 0
                || strcmp("8", cinfo.outtype) == 0
                || strcmp("10", cinfo.outtype) == 0
                || strcmp("11", cinfo.outtype) == 0
                || strcmp("12", cinfo.outtype) == 0
                || strcmp("17", cinfo.outtype) == 0
                || strcmp("18", cinfo.outtype) == 0
                || strcmp("20", cinfo.outtype) == 0) {
		connect_flag = 0;
		bondtype_flag = 0;
		atomtype_flag = 0;
	}

        if (strcmp("mopout", cinfo.outtype) == 0
                || strcmp("gout", cinfo.outtype) == 0
                || strcmp("jout", cinfo.outtype) == 0
                || strcmp("rst", cinfo.outtype) == 0
                || strcmp("11", cinfo.outtype) == 0
                || strcmp("12", cinfo.outtype) == 0
                || strcmp("17", cinfo.outtype) == 0
                || strcmp("20", cinfo.outtype) == 0)
                duplicatedname_flag = 0;

        if (strcmp("mopint", cinfo.outtype) == 0
                || strcmp("gzmat", cinfo.outtype) == 0 
                || strcmp("jzmat", cinfo.outtype) == 0 
                || strcmp("7", cinfo.outtype) == 0 
                || strcmp("9", cinfo.outtype) == 0 
                || strcmp("19", cinfo.outtype) == 0 ) {
		bondtype_flag = 0;
		atomtype_flag = 0;
	}
/*     	the following code judge or assign atom name, atom type, bond type etc according to flags */
	if (adjustatomname_flag) {
		if(strcmp(cinfo.intype, "mol2")==0 || strcmp(cinfo.intype, "2")==0 ||
		   strcmp(cinfo.intype, "ac")==0 || strcmp(cinfo.intype, "1")==0)
			adjustatomname(atomnum, atom, 1);
		else
			adjustatomname(atomnum, atom, 0);
	}
	if (atomicnum_flag) 
		atomicnum(atomnum, atom);
	if (atomname_flag) 
		atomname(atomnum, atom);
	if (default_flag)
		default_inf(atomnum, atom, default_flag);
	if (cartcoord_flag)
		cartcoord(atomnum, atom);
	if (connect_flag) {
		overflow_flag =
			connect(minfo.connect_file, atomnum, atom, &bondnum, bond,
					cinfo.maxbond);
		if (overflow_flag) {
			cinfo.maxbond = bondnum + 10;
			memory(2, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag =
				connect(minfo.connect_file, atomnum, atom, &bondnum, bond,
						cinfo.maxbond);
		}
	}
	if (bondtype_flag && bondnum > 0) {
		judgebondtype(atomnum, atom, bondnum, bond, cinfo, minfo,
					  bondtype_flag);
		if(cinfo.prediction_index ==2||cinfo.prediction_index ==3) {
			cinfo.prediction_index = 0;
			atomtype_flag = 0;
			bondtype_flag = 0;
		}
		if(cinfo.prediction_index ==4||cinfo.prediction_index ==5) {
			cinfo.prediction_index = 1;
			atomtype_flag = 1;
			bondtype_flag = 0;
		}
	}

	if (duplicatedname_flag)
		duplicatedname(atomnum, atom);
	if (atomtype_flag) {
		wac("ANTECHAMBER_AC.AC0", atomnum, atom, bondnum, bond, cinfo,
			minfo);
                copied_size = build_exe_path(tmpchar, "atomtype"
                        " -i ANTECHAMBER_AC.AC0 -o ANTECHAMBER_AC.AC -p "
                        , MAXCHAR );
                strncat(tmpchar, minfo.atom_type_def, MAXCHAR - copied_size );

		if (cinfo.intstatus == 2)
			fprintf(stdout, "\nRunning: %s\n", tmpchar);
		status = system(tmpchar);
	        if(status != 0) {
                	fprintf(stderr, "Error: cannot run \"%s\" in main() of antechamber.c properly, exit\n", tmpchar);
                	exit(1);
       		}
		rac("ANTECHAMBER_AC.AC", &atomnum, atom, &bondnum, bond, &cinfo,
			&minfo);
	}

/*     the following code readin or calculate charges */
/* 	usercharge info*/
	if (minfo.usercharge > -9999) {	/*read in charge with -nc flag */
		minfo.icharge = minfo.usercharge;
		minfo.dcharge = minfo.usercharge;
	}
	if(minfo.dcharge < -9990) {
		minfo.dcharge = 0.0;
        	for (i = 0; i < atomnum; i++)
                	minfo.dcharge += atom[i].charge;
        		fraction = modf(minfo.dcharge, &tmpf);
        		minfo.icharge = (int) tmpf;
        		if (fabs(fraction) >= 0.50) {
                		if (minfo.dcharge < 0)
                        		minfo.icharge--;
                		if (minfo.dcharge > 0)
                        		minfo.icharge++;
        		}
	}	
/*zero weird charges */
	if(minfo.usercharge < -9990 && (minfo.icharge <= -10 || minfo.icharge >= 10)) {
		fprintf(stderr, "Warning: Weird total charge: %d!"
            "The net charge is assumed to be 0.\n"
		      "         If the weird charge was correct, "
            "specify it via the -nc net charge flag.\n", minfo.icharge);
		minfo.icharge = 0;	
	} 
	if (strcmp("resp", cinfo.chargetype) == 0
		|| strcmp("1", cinfo.chargetype) == 0)
		resp(ifilename, atomnum, atom, bondnum, bond, cinfo, minfo);
	if (strcmp("bcc", cinfo.chargetype) == 0
		|| strcmp("2", cinfo.chargetype) == 0)
		bcc(ifilename, atomnum, atom, bondnum, bond, arom, &cinfo, &minfo);
	if (strcmp("cm1", cinfo.chargetype) == 0
		|| strcmp("3", cinfo.chargetype) == 0)
		cm1(atomnum, atom, &cinfo, &minfo);
	if (strcmp("cm2", cinfo.chargetype) == 0
		|| strcmp("4", cinfo.chargetype) == 0)
		cm2(atomnum, atom, &cinfo, &minfo);
	if (strcmp("esp", cinfo.chargetype) == 0
		|| strcmp("5", cinfo.chargetype) == 0)
		esp(ifilename, atomnum, atom, cinfo, minfo);
	if (strcmp("mul", cinfo.chargetype) == 0
		|| strcmp("6", cinfo.chargetype) == 0)
		mul(ifilename, atomnum, atom, &cinfo, &minfo);
	if (strcmp("gas", cinfo.chargetype) == 0
		|| strcmp("7", cinfo.chargetype) == 0)
		gascharge(atomnum, atom, bondnum, bond, cinfo, &minfo);
	if (strcmp("rc", cinfo.chargetype) == 0
		|| strcmp("8", cinfo.chargetype) == 0)
		rcharge(cfilename, atomnum, atom, cinfo, &minfo);
	if (strcmp("wc", cinfo.chargetype) == 0
		|| strcmp("9", cinfo.chargetype) == 0)
		wcharge(cfilename, atomnum, atom, cinfo, minfo);
	if (strcmp("dc", cinfo.chargetype) == 0
		|| strcmp("10", cinfo.chargetype) == 0) {
		for(i=0;i<atomnum;i++)
			atom[i].charge = 0.0;
	}

/*	read the radii*/
	if (strcmp("mpdb", cinfo.intype) != 0 && strcmp("4", cinfo.intype) != 0 &&
		(strcmp("mpdb", cinfo.outtype) == 0 || strcmp("4", cinfo.outtype) == 0))
		read_radius(minfo.radius_file, atomnum, atom);

/*    	read in additional files*/
/*	expand the array size a little bit; 
  	expand the array size a little bit 
  	since atom file type such as prepi may have additional atom records  
*/
	if (strlen(afilename) > 1 && strlen(cinfo.atype) > 1) {
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
			bond[i].jflag = -1; /* bond type has not been assigned */
		}

		if (strcmp("ac", cinfo.atype) == 0
			|| strcmp("1", cinfo.atype) == 0) {
			overflow_flag =
				rac(afilename, &atomnum_tmp, atom_tmp, &bondnum_tmp,
					bond_tmp, &cinfo, &minfo);
			if (overflow_flag) {
				fprintf(stderr,
						"Overflow happens for additional files, exit");
				exit(1);
			}
		}
		if (strcmp("mol2", cinfo.atype) == 0
			|| strcmp("2", cinfo.atype) == 0) {
			overflow_flag =
				rmol2(afilename, &atomnum_tmp, atom_tmp, &bondnum_tmp,
					  bond_tmp, &cinfo, &minfo, 0);
			if (overflow_flag) {
				fprintf(stderr,
						"Overflow happens for additional files, exit");
				exit(1);
			}
		}
		if (strcmp("mopint", cinfo.atype) == 0
			|| strcmp("9", cinfo.atype) == 0) {
			overflow_flag =
				rmopint(afilename, &atomnum_tmp, atom_tmp, cinfo, minfo);
			if (overflow_flag) {
				fprintf(stderr,
						"Overflow happens for additional files, exit");
				exit(1);
			}
		}
		if (strcmp("mopcrt", cinfo.atype) == 0
			|| strcmp("10", cinfo.atype) == 0) {
			overflow_flag =
				rmopcrt(afilename, &atomnum_tmp, atom_tmp, cinfo, minfo);
			if (overflow_flag) {
				fprintf(stderr,
						"Overflow happens for additional files, exit");
				exit(1);
			}
		}

		if (strcmp("mopout", cinfo.atype) == 0
			|| strcmp("12", cinfo.atype) == 0) {
			overflow_flag =
				rmopout(afilename, &atomnum_tmp, atom_tmp, &cinfo, &minfo);
			if (overflow_flag) {
				fprintf(stderr,
						"Overflow happens for additional files, exit");
				exit(1);
			}
		}
		if (strcmp("gcrt", cinfo.atype) == 0
			|| strcmp("8", cinfo.atype) == 0) {
			overflow_flag =
				rgcrt(afilename, &atomnum_tmp, atom_tmp, cinfo, minfo);
			if (overflow_flag) {
				fprintf(stderr,
						"Overflow happens for additional files, exit");
				exit(1);
			}
		}
		if (strcmp("gzmat", cinfo.atype) == 0
			|| strcmp("7", cinfo.atype) == 0) {
			overflow_flag =
				rgzmat(afilename, &atomnum_tmp, atom_tmp, cinfo, minfo);
			if (overflow_flag) {
				fprintf(stderr,
						"Overflow happens for additional files, exit");
				exit(1);
			}
		}

		if (strcmp("gout", cinfo.atype) == 0
			|| strcmp("11", cinfo.atype) == 0) {
			overflow_flag =
				rgout(afilename, &atomnum_tmp, atom_tmp, cinfo, &minfo);
			if (overflow_flag) {
				fprintf(stderr,
						"Overflow happens for additional files, exit");
				exit(1);
			}
		}

		if (strcmp("pdb", cinfo.atype) == 0
			|| strcmp("3", cinfo.atype) == 0) {
			overflow_flag =
				rpdb(afilename, &atomnum_tmp, atom_tmp, cinfo, minfo, 0);
			if (overflow_flag) {
				fprintf(stderr,
						"Overflow happens for additional files, exit");
				exit(1);
			}
		}

		if (strcmp("mpdb", cinfo.atype) == 0
			|| strcmp("4", cinfo.atype) == 0) {
			overflow_flag =
				rpdb(afilename, &atomnum_tmp, atom_tmp, cinfo, minfo, 1);
			if (overflow_flag) {
				fprintf(stderr,
						"Overflow happens for additional files, exit");
				exit(1);
			}
		}

		if (strcmp("csd", cinfo.atype) == 0
			|| strcmp("14", cinfo.atype) == 0) {
			overflow_flag =
				rcsd(afilename, &atomnum_tmp, atom_tmp, cinfo, minfo);
			if (overflow_flag) {
				fprintf(stderr,
						"Overflow happens for additional files, exit");
				exit(1);
			}
		}

		if (strcmp("mdl", cinfo.atype) == 0
			|| strcmp("sd", cinfo.atype) == 0
			|| strcmp("sdf", cinfo.atype) == 0
			|| strcmp("15", cinfo.atype) == 0) {
			overflow_flag =
				rmdl(afilename, &atomnum_tmp, atom_tmp, &bondnum_tmp,
					 bond_tmp, cinfo, minfo);
			if (overflow_flag) {
				fprintf(stderr,
						"Overflow happens for additional files, exit");
				exit(1);
			}
		}

		if (strcmp("alc", cinfo.atype) == 0
			|| strcmp("13", cinfo.atype) == 0) {
			overflow_flag =
				ralc(afilename, &atomnum_tmp, atom_tmp, &bondnum_tmp,
					 bond_tmp, cinfo, minfo);
			if (overflow_flag) {
				fprintf(stderr,
						"Overflow happens for additional files, exit");
				exit(1);
			}
		}

		if (strcmp("hin", cinfo.atype) == 0
			|| strcmp("16", cinfo.atype) == 0) {
			overflow_flag =
				rhin(afilename, &atomnum_tmp, atom_tmp, &bondnum_tmp,
					 bond_tmp, cinfo, minfo);
			if (overflow_flag) {
				fprintf(stderr,
						"Overflow happens for additional files, exit");
				exit(1);
			}
		}

		if (strcmp("prepi", cinfo.atype) == 0
			|| strcmp("5", cinfo.atype) == 0) {
			overflow_flag =
				rprepi(afilename, &atomnum_tmp, atom_tmp, &bondnum_tmp, bond_tmp, &cinfo, &minfo);
			if (overflow_flag) {
				fprintf(stderr,
						"Overflow happens for additional files, exit");
				exit(1);
			}
		}

		if (strcmp("prepc", cinfo.atype) == 0
			|| strcmp("6", cinfo.atype) == 0) {
			overflow_flag =
				rprepc(afilename, &atomnum_tmp, atom_tmp, &cinfo, &minfo);
			if (overflow_flag) {
				fprintf(stderr,
						"Overflow happens for additional files, exit");
				exit(1);
			}
		}
		if (strcmp("rst", cinfo.atype) == 0
			|| strcmp("17", cinfo.atype) == 0) {
			overflow_flag = rrst(afilename, &atomnum_tmp, atom_tmp, cinfo);
			if (overflow_flag) {
				fprintf(stderr,
						"Overflow happens for additional files, exit");
				exit(1);
			}
		}
                if (strcmp("divcrt", cinfo.atype) == 0
                        || strcmp("10", cinfo.atype) == 0) {
                        overflow_flag =
                                rdivcrt(afilename, &atomnum_tmp, atom_tmp, cinfo, minfo);
                        if (overflow_flag) {
                                fprintf(stderr,
                                                "Overflow happens for additional files, exit");
                                exit(1);
                        }
                }
                                                                                                                                                                                                           
                if (strcmp("divout", cinfo.atype) == 0
                        || strcmp("12", cinfo.atype) == 0) {
                        overflow_flag =
                                rdivout(afilename, &atomnum_tmp, atom_tmp, &cinfo, &minfo);
                        if (overflow_flag) {
                                fprintf(stderr,
                                                "Overflow happens for additional files, exit");
                                exit(1);
                        }
                }
		if (ao_flag == 1)
			for (i = 0; i < atomnum; i++) {
				atom[i].x = atom_tmp[i].x;
				atom[i].y = atom_tmp[i].y;
				atom[i].z = atom_tmp[i].z;
			}

		if (ao_flag == 2)
			for (i = 0; i < atomnum; i++)
				atom[i].charge = atom_tmp[i].charge;

		if (ao_flag == 3)
			for (i = 0; i < atomnum; i++)
				strcpy(atom[i].name, atom_tmp[i].name);
		if (ao_flag == 4)
			for (i = 0; i < atomnum; i++)
				strcpy(atom[i].ambername, atom_tmp[i].ambername);
		if (ao_flag == 5)
			for (i = 0; i < bondnum; i++)
				bond[i].type = bond_tmp[i].type;
		free(atom_tmp);
		free(bond_tmp);
	}

/*	This part may not necessary! commented by Junmei Wang
	if(atomtype_flag == 0 && (strcmp("charmm", cinfo.outtype) == 0|| strcmp("23", cinfo.outtype) == 0  ||
	                         strcmp("ac", cinfo.outtype) == 0 || strcmp("1", cinfo.outtype) == 0 ||
	                         strcmp("mol2", cinfo.outtype) == 0 || strcmp("2", cinfo.outtype) == 0 ||
	                         strcmp("csd", cinfo.outtype) == 0 || strcmp("14", cinfo.outtype) == 0 ||
	                         strcmp("mdl", cinfo.outtype) == 0 || strcmp("15", cinfo.outtype) == 0 ||
	                         strcmp("sdf", cinfo.outtype) == 0 || strcmp("sd", cinfo.outtype) == 0 ||
	                         strcmp("alc", cinfo.outtype) == 0 || strcmp("13", cinfo.outtype) == 0 ||
	                         strcmp("hin", cinfo.outtype) == 0 || strcmp("16", cinfo.outtype) == 0)) {
        	memory(8, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
        	overflow_flag =
			ringdetect(atomnum, atom, bondnum, bond, &ringnum, ring, arom,
			  	 cinfo.maxatom, cinfo.maxring, minfo.inf_filename, 0);
        	if (overflow_flag) {
			cinfo.maxring = ringnum + 10;
                	memory(8, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
        		overflow_flag =
				ringdetect(atomnum, atom, bondnum, bond, &ringnum, ring, arom,
			   		cinfo.maxatom, cinfo.maxring, minfo.inf_filename, 0);
		}
	}
*/

	if(ra_flag==1)	
		read_at(at_filename);
	if(rb_flag==1)	
		read_bt(bt_filename);

/*      write out files */
	if (strcmp("ac", cinfo.outtype) == 0
		|| strcmp("1", cinfo.outtype) == 0)
		wac(ofilename, atomnum, atom, bondnum, bond, cinfo, minfo);
	if (strcmp("charmm", cinfo.outtype) == 0
		|| strcmp("23", cinfo.outtype) == 0)
                wcharmm(ofilename, ifilename, atomnum, atom, bondnum, bond, cinfo, minfo);
	if (strcmp("mol2", cinfo.outtype) == 0
		|| strcmp("2", cinfo.outtype) == 0)  
		wmol2(ofilename, atomnum, atom, bondnum, bond, arom, cinfo, minfo);
	if (strcmp("mopint", cinfo.outtype) == 0
		|| strcmp("9", cinfo.outtype) == 0) {
		adjust_sequence_order(atomnum, atom, bondnum, bond);
		wmopint(ofilename, atomnum, atom, minfo);
	}
	if (strcmp("mopcrt", cinfo.outtype) == 0
		|| strcmp("10", cinfo.outtype) == 0)
		wmopcrt(ofilename, atomnum, atom, minfo);
	if (strcmp("mopout", cinfo.outtype) == 0
		|| strcmp("12", cinfo.outtype) == 0)
		wmopout(ofilename, atomnum, atom, cinfo, minfo);
	if (strcmp("gcrt", cinfo.outtype) == 0
		|| strcmp("8", cinfo.outtype) == 0)
		wgcrt(ofilename, atomnum, atom, minfo);
	if (strcmp("gzmat", cinfo.outtype) == 0
		|| strcmp("7", cinfo.outtype) == 0) {
		adjust_sequence_order(atomnum, atom, bondnum, bond);
		wgzmat(ofilename, atomnum, atom, minfo);
	}
	if (strcmp("gout", cinfo.outtype) == 0
		|| strcmp("11", cinfo.outtype) == 0)
		wgout();
	if (strcmp("jcrt", cinfo.outtype) == 0
		|| strcmp("18", cinfo.outtype) == 0)
		wjcrt(ofilename, atomnum, atom, minfo);
	if (strcmp("jzmat", cinfo.outtype) == 0
		|| strcmp("19", cinfo.outtype) == 0) {
		adjust_sequence_order(atomnum, atom, bondnum, bond);
		wjzmat(ofilename, atomnum, atom, minfo);
	}
	if (strcmp("jout", cinfo.outtype) == 0
		|| strcmp("20", cinfo.outtype) == 0)
		wjout();
	if (strcmp("pdb", cinfo.outtype) == 0
		|| strcmp("3", cinfo.outtype) == 0)
		wpdb(ofilename, atomnum, atom);
	if (strcmp("mpdb", cinfo.outtype) == 0
		|| strcmp("4", cinfo.outtype) == 0)
		wmpdb(ofilename, atomnum, atom);
	if (strcmp("csd", cinfo.outtype) == 0
		|| strcmp("14", cinfo.outtype) == 0)
		wcsd(ofilename, atomnum, atom);
	if (strcmp("mdl", cinfo.outtype) == 0
		|| strcmp("sdf", cinfo.outtype) == 0
		|| strcmp("sd", cinfo.outtype) == 0
		|| strcmp("15", cinfo.outtype) == 0)
		wmdl(ofilename, atomnum, atom, bondnum, bond, cinfo);
	if (strcmp("alc", cinfo.outtype) == 0
		|| strcmp("13", cinfo.outtype) == 0)
		walc(ofilename, atomnum, atom, bondnum, bond);
	if (strcmp("hin", cinfo.outtype) == 0
		|| strcmp("16", cinfo.outtype) == 0)
		whin(ofilename, atomnum, atom, bondnum, bond);
	if (strcmp("prepi", cinfo.outtype) == 0
		|| strcmp("5", cinfo.outtype) == 0)  
		wprep(ofilename, ifilename, atomnum, atom, bondnum, bond, cinfo,
			  &minfo, 1);
	if (strcmp("prepc", cinfo.outtype) == 0
		|| strcmp("6", cinfo.outtype) == 0)
		wprep(ofilename, ifilename, atomnum, atom, bondnum, bond, cinfo,
			  &minfo, 0);
	if (strcmp("rst", cinfo.outtype) == 0
		|| strcmp("17", cinfo.outtype) == 0)
		wrst(ofilename, atomnum, atom);
	if (strcmp("divcrt", cinfo.outtype) == 0
		|| strcmp("21", cinfo.outtype) == 0)
		wdivcrt(ofilename, atomnum, atom, minfo);
	if (strcmp("divout", cinfo.outtype) == 0
		|| strcmp("22", cinfo.outtype) == 0)
		wdivout(ofilename, atomnum, atom, cinfo, minfo);
/*	info(atomnum, atom, bondnum, bond, arom, cinfo, minfo); */
/*
        free(atom);
        free(bond);
        free(arom);
*/
	if(wa_flag==1)	
		write_at(at_filename);
	if(wb_flag==1)	
		write_bt(bt_filename);
	if (cinfo.pfindex == 1) {
		status = system("rm -f ANTECHAMBER* ATOMTYPE.INF BCCTYPE.INF NEWPDB.PDB PREP.INF");
	        if(status != 0) {
                	fprintf(stderr, "Error: cannot run \"%s\" in main() of antechamber.c properly, exit\n", "rm -f ANTECHAMBER* ATOMTYPE.INF BCCTYPE.INF NEWPDB.PDB PREP.INF");
                	exit(1);
        	}
	
	}
	printf("\n");
	return (0);
}
