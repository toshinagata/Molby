%{
extern int mytaskid;

#ifdef flex

#undef YY_INPUT
#define YY_INPUT(b, r, ms) (r = mmo_yyinput(b, ms))

#else

#undef input
#undef unput
#define input() (*gopts++)
#define unput(c) (*--gopts = c)

#endif

#undef ECHO
#define ECHO yyout=stdout;if(mytaskid==0) fprintf(yyout, "\tmm_options:  %s\n",yytext)

%}
%%

cut[=\ ][0-9eE\.\+\-]+		{ ECHO; cut  = atof( &yytext[4] ); }
cutnp[=\ ][0-9eE\.\+\-]+	{ ECHO; cutnp  = atof( &yytext[6] ); }
scnb[=\ ][0-9eE\.\+\-]+		{ ECHO; scnb = atof( &yytext[5] ); }
scee[=\ ][0-9eE\.\+\-]+		{ ECHO; scee = atof( &yytext[5] ); }
ntpr[=\ ][1-9][0-9]*		{ ECHO; ntpr = atoi(&yytext[5]); 
                                    if( ntpr==0) ntpr = 9999999; }
hcp[=\ ][012]			{ ECHO; hcp = atoi( &yytext[4] ); }
hcp_h[=\ ][0-9eE\.\+\-]+	{ ECHO; hcp_h  = atof( &yytext[6] ); }
nsnb[=\ ][0-9]+			{ ECHO; nsnb = atoi( &yytext[5] ); }
nsnp[=\ ][0-9]+			{ ECHO; nsnp = atoi( &yytext[5] ); }
ntwx[=\ ][0-9]+			{ ECHO; ntwx = atoi( &yytext[5] ); }
nchk[=\ ][0-9]+			{ ECHO; nchk = atoi( &yytext[5] ); }
nchk2[=\ ][0-9]+		{ ECHO; nchk2 = atoi( &yytext[6] ); }

dt[=\ ][0-9eE\.\+\-]+		{ ECHO; dt  = atof( &yytext[3] ); }
t[=\ ][0-9eE\.\+\-]+		{ ECHO; t = atof( &yytext[2] ); }
tautp[=\ ][0-9eE\.\+\-]+	{ ECHO; tautp = atof( &yytext[6] ); }
temp0[=\ ][0-9eE\.\+\-]+	{ ECHO; temp0 = atof( &yytext[6] ); }
tempi[=\ ][0-9eE\.\+\-]+	{ ECHO; tempi = atof( &yytext[6] ); }
vlimit[=\ ][0-9eE\.\+\-]+	{ ECHO; vlimit = atof( &yytext[7] ); }
gamma_ln[=\ ][0-9eE\.\+\-]+	{ ECHO; gamma_ln = atof( &yytext[9] ); }
genmass[=\ ][0-9eE\.\+\-]+	{ ECHO; genmass = atof( &yytext[8] ); }
wcons[=\ ][0-9eE\.\+\-]+	{ ECHO; wcons = atof( &yytext[6] ); }
ntpr_md[=\ ][1-9][0-9]*		{ ECHO; ntpr_md = atoi(&yytext[8]); 
                                    if( ntpr_md==0) ntpr_md = 9999999; }
zerov[=\ ][0-9]+		{ ECHO; zerov = atoi( &yytext[6] ); }
k4d[=\ ][0-9eE\.\+\-]+		{ ECHO; k4d  = atof( &yytext[4] ); }

ips[=\ ][0123]			{ ECHO; ips = atoi( &yytext[4] ); }

gb[=\ ][01258]			{ ECHO; gb = atoi( &yytext[3] ); }
gb_debug[=\ ][01]		{ ECHO; gb_debug = atoi( &yytext[9] ); }
nr_debug[=\ ][01]		{ ECHO; nr_debug = atoi( &yytext[9] ); }
gb2_debug[=\ ][01]		{ ECHO; gb2_debug = atoi( &yytext[10] ); }
gbsa_debug[=\ ][01]		{ ECHO; gbsa_debug = atoi( &yytext[11] ); }
e_debug[=\ ][01]		{ ECHO; e_debug = atoi( &yytext[8] ); }
epsext[=\ ][0-9eE\.\+\-]+	{ ECHO; epsext = atof( &yytext[7] ); }
kappa[=\ ][0-9eE\.\+\-]+	{ ECHO; kappa = atof( &yytext[6] ); }
rgbmax[=\ ][0-9eE\.\+\-]+	{ ECHO; rgbmax = atof( &yytext[7] ); }

surften[=\ ][0-9eE\.\+\-]+	{ ECHO; surften = atof( &yytext[8] ); }
dradius[=\ ][0-9eE\.\+\-]+	{ ECHO; dradius = atof( &yytext[8] ); }
deltar[=\ ][0-9eE\.\+\-]+	{ ECHO; deltar = atof( &yytext[7] ); }
rwater[=\ ][0-9eE\.\+\-]+	{ ECHO; rwater = atof( &yytext[7] ); }
alphanp[=\ ][0-9eE\.\+\-]+	{ ECHO; alphanp= atof( &yytext[8] ); }
kappanp[=\ ][0-9eE\.\+\-]+	{ ECHO; kappanp = atof( &yytext[8] ); }
min_volume[=\ ][0-9eE\.\+\-]+	{ ECHO; min_volume  = atof( &yytext[11] ); }
max_set_size[=\ ][0-9]+		{ ECHO; max_set_size = atoi( &yytext[13] ); }
MPI_min_tasks[=\ ][0-9]+	{ ECHO; MPI_min_tasks = atoi( &yytext[14] ); }
cull_np_lists[=\ ][01]  	{ ECHO; cull_np_lists = atoi( &yytext[14] ); }
use_lower_tri[=\ ][01]  	{ ECHO; use_lower_tri = atoi( &yytext[14] ); }
dynamic_loops[=\ ][01]  	{ ECHO; dynamic_loops = atoi( &yytext[14] ); }
static_arrays[=\ ][01]  	{ ECHO; static_arrays = atoi( &yytext[14] ); }
blocksize[=\ ][0-9]+		{ ECHO; blocksize = atoi( &yytext[10] ); }
gbsa[=\ ][0-5]	                { ECHO; gbsa = atoi( &yytext[5] ); }
dim[=\ ][34]                    { ECHO; dim = atoi( &yytext[4] ); }

diel[=\ ]Yammp			{ ECHO; dield = -5; }
diel[=\ ]R94			{ ECHO; dield = -4; }
diel[=\ ]RL94			{ ECHO; dield = -3; }
diel[=\ ]RL			{ ECHO; dield = -2; }
diel[=\ ]R			{ ECHO; dield =  0; }
diel[=\ ]C			{ ECHO; dield =  1; }
dielc[=\ ][0-9eE\.\+\-]+	{ ECHO; dielc = atof( &yytext[6] ); }

rattle[=\ ][012]			{ ECHO; irattle = atoi( &yytext[7] ); }

[\ \n\t,]+						;

.						{	
#ifdef flex
                if(mytaskid==0)fprintf( stderr, "Error: mm_option: unknown keyword or invalid value:\n" );
		if(mytaskid==0)fprintf( stderr, "  first character of the option is '%s'.\n", yytext );
#else
		if(mytaskid==0)fprintf( stderr, "Error: mm_option: unknown keyword or invalid value: %s\n", --gopts );
#endif

		exit( 1 );
	}

%%

int yywrap(){
	return( 1 );
}

#ifdef flex


int  mmo_yyinput( buf, max_size )

char *buf;
int  max_size;
{
	int n;

	if( max_size > &gopts[mmoinputlim] - mmoinputptr ) 
		n =  &gopts[mmoinputlim] - mmoinputptr;
	else
		n = max_size;

	if( n>0 ){
		memcpy( buf, mmoinputptr, n );
		mmoinputptr += n;
	}
	return n;
}

#endif
