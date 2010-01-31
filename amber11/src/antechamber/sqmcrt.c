void wsqmin(char *filename, int atomnum, ATOM atom[], MOLINFO minfo)
{
	FILE *fpout;
	int i, nelectrons; 

	if ((fpout = fopen(filename, "w")) == NULL) {
		fprintf(stderr, "Cannot open file %s to write in wsqmcrt(), exit\n", filename);
		exit(1);
	}

    /*  write initial keywords and control parameters:  */
	fprintf(fpout, "Run semi-empirical minimization\n");
	fprintf(fpout, " &qmmm\n" );
	fprintf(fpout, "  %s  qmcharge=%d,\n /\n",
       minfo.ekeyword, minfo.icharge );

	/* element(atomnum, atom); */

	nelectrons = 0;
	for (i = 0; i < atomnum; i++){
		fprintf(fpout, "%4d %5s  %12.4lf  %12.4lf  %12.4lf \n",
			atom[i].atomicnum, atom[i].name, atom[i].x, atom[i].y, atom[i].z);
		nelectrons += atom[i].atomicnum;
	}
	fprintf(fpout, "\n");
	fclose(fpout);
	nelectrons -= minfo.icharge;
	fprintf( stderr, "Total number of electrons: %d; net charge: %d\n",
		nelectrons,minfo.icharge );
    /*  check that the number of electrons is even:   */
	if( nelectrons%2 != 0 ){
		fprintf( stderr, "INFO: Number of electrons is odd: %d\n", nelectrons) ;
		fprintf( stderr, "      Please check the total charge (-nc flag) and spin multiplicity (-m flag)\n"); 
	}
}
