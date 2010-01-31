/* arpkdrv.f -- translated by f2c (version 20030306).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Common Block Declarations */

struct {
    integer logfil, ndigit, mgetv0, msaupd, msaup2, msaitr, mseigt, msapps, 
	    msgets, mseupd, mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, 
	    mneupd, mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd;
} debug_;

#define debug_1 debug_

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__3 = 3;
static integer c__6 = 6;
static integer c__2 = 2;
static integer c_n6 = -6;
static integer c__5 = 5;

/*    Original ARPACK/EXAMPLE/SIMPLE driver program rewritten by */
/*    Istvan Kolossvary to link with libLMOD. */

/* Subroutine */ int dsarpack_(integer *n_dim__, integer *n_eig_in__, integer 
	*n_eig_out__, integer *ncv_in__, integer *itr_in__, doublereal *
	eigval_tol__, doublereal *eigvals, doublereal *eigvecs, integer *
	spectrum, integer *need_eigvecs__, integer *ierr, integer *
	debug_arpack__, doublereal *v, doublereal *workl, doublereal *workd, 
	doublereal *d__, doublereal *resid, doublereal *ax, logical *select, 
	doublereal *xyz, doublereal *grad, integer *return_flag__, integer *
	label)
{
    /* Initialized data */

    static integer l12 = 1;
    static integer l18 = 2;
    static integer arpack_error__ = -2;

    /* System generated locals */
    integer v_dim1, v_offset, d_dim1, d_offset, i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);

    /* Local variables */
    static integer i__, j, n, ido, ldv, ncv, nev;
    static doublereal tol;
    static integer status_flag__;
    static char bmat[1];
    static integer info;
    static logical rvec;
    static integer maxn, mode1;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static doublereal sigma;
    static char which[2];
    static integer nconv;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *), dmout_(integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, char *, ftnlen);
    static integer ipntr[11], iparam[11];
    extern /* Subroutine */ int dsaupd_(integer *, char *, integer *, char *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, ftnlen, ftnlen), dseupd_(logical *, char *, 
	    logical *, doublereal *, doublereal *, integer *, doublereal *, 
	    char *, integer *, char *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, ftnlen, ftnlen, 
	    ftnlen);
    static integer maxncv, maxnev, ishfts, maxitr, lworkl;
    extern /* Subroutine */ int hessvec_(integer *, doublereal *, doublereal *
	    , doublereal *, doublereal *, integer *, integer *);

    /* Fortran I/O blocks */
    static cilist io___13 = { 0, 6, 0, 0, 0 };
    static cilist io___14 = { 0, 6, 0, 0, 0 };
    static cilist io___15 = { 0, 6, 0, 0, 0 };
    static cilist io___16 = { 0, 6, 0, 0, 0 };
    static cilist io___27 = { 0, 6, 0, 0, 0 };
    static cilist io___28 = { 0, 6, 0, 0, 0 };
    static cilist io___29 = { 0, 6, 0, 0, 0 };
    static cilist io___30 = { 0, 6, 0, 0, 0 };
    static cilist io___33 = { 0, 6, 0, 0, 0 };
    static cilist io___34 = { 0, 6, 0, 0, 0 };
    static cilist io___35 = { 0, 6, 0, 0, 0 };
    static cilist io___36 = { 0, 6, 0, 0, 0 };
    static cilist io___38 = { 0, 6, 0, 0, 0 };
    static cilist io___39 = { 0, 6, 0, 0, 0 };
    static cilist io___40 = { 0, 6, 0, 0, 0 };
    static cilist io___42 = { 0, 6, 0, 0, 0 };
    static cilist io___43 = { 0, 6, 0, 0, 0 };
    static cilist io___44 = { 0, 6, 0, 0, 0 };
    static cilist io___45 = { 0, 6, 0, 0, 0 };
    static cilist io___46 = { 0, 6, 0, 0, 0 };
    static cilist io___47 = { 0, 6, 0, 0, 0 };
    static cilist io___48 = { 0, 6, 0, 0, 0 };
    static cilist io___49 = { 0, 6, 0, 0, 0 };
    static cilist io___50 = { 0, 6, 0, 0, 0 };
    static cilist io___51 = { 0, 6, 0, 0, 0 };
    static cilist io___52 = { 0, 6, 0, 0, 0 };
    static cilist io___53 = { 0, 6, 0, 0, 0 };
    static cilist io___54 = { 0, 6, 0, 0, 0 };
    static cilist io___55 = { 0, 6, 0, 0, 0 };
    static cilist io___56 = { 0, 6, 0, 0, 0 };
    static cilist io___57 = { 0, 6, 0, 0, 0 };
    static cilist io___58 = { 0, 6, 0, 0, 0 };
    static cilist io___59 = { 0, 6, 0, 0, 0 };
    static cilist io___60 = { 0, 6, 0, 0, 0 };
    static cilist io___61 = { 0, 6, 0, 0, 0 };
    static cilist io___62 = { 0, 6, 0, 0, 0 };
    static cilist io___63 = { 0, 6, 0, 0, 0 };




/*     %-----------------% */
/*     | Dummy Arguments | */
/*     %-----------------% */



/*     %---------------% */
/*     | Include Files | */
/*     %---------------% */


/*     This code shows how to use ARPACK to find a few eigenvalues */
/*     (lambda) and corresponding eigenvectors (x) for the standard */
/*     eigenvalue problem: */

/*                        A*x = lambda*x */

/*     where A is an n by n real symmetric matrix. */

/*     The main points illustrated here are */

/*        1) How to declare sufficient memory to find NEV */
/*           eigenvalues of largest magnitude.  Other options */
/*           are available. */

/*        2) Illustration of the reverse communication interface */
/*           needed to utilize the top level ARPACK routine DSAUPD */
/*           that computes the quantities needed to construct */
/*           the desired eigenvalues and eigenvectors(if requested). */

/*        3) How to extract the desired eigenvalues and eigenvectors */
/*           using the ARPACK routine DSEUPD. */

/*     The only thing that must be supplied in order to use this */
/*     routine on your problem is to change the array dimensions */
/*     appropriately, to specify WHICH eigenvalues you want to compute */
/*     and to supply a matrix-vector product */

/*                         w <-  Av */

/*     in place of the call to AV( ) below. */

/*     Once usage of this routine is understood, you may wish to explore */
/*     the other available options to improve convergence, to solve generalized */
/*     problems, etc.  Look at the file ex-sym.doc in DOCUMENTS directory. */
/*     This codes implements */

/* \Example-1 */
/*     ... Suppose we want to solve A*x = lambda*x in regular mode, */
/*         where A is derived from the central difference discretization */
/*         of the 2-dimensional Laplacian on the unit square with */
/*         zero Dirichlet boundary condition. */
/*     ... OP = A  and  B = I. */
/*     ... Assume "call av (n,x,y)" computes y = A*x */
/*     ... Use mode 1 of DSAUPD. */

/* \BeginLib */

/* \Routines called: */
/*     dsaupd  ARPACK reverse communication interface routine. */
/*     dseupd  ARPACK routine that returns Ritz values and (optionally) */
/*             Ritz vectors. */
/*     dnrm2   Level 1 BLAS that computes the norm of a vector. */
/*     daxpy   Level 1 BLAS that computes y <- alpha*x+y. */

/* \Author */
/*     Richard Lehoucq */
/*     Danny Sorensen */
/*     Chao Yang */
/*     Dept. of Computational & */
/*     Applied Mathematics */
/*     Rice University */
/*     Houston, Texas */

/* \SCCS Information: %Z% */
/* FILE: %M%   SID: %I%   DATE OF SID: %G%   RELEASE: %R% */

/* \Remarks */
/*     1. None */

/* \EndLib */

/* ----------------------------------------------------------------------- */

/*     %-------------------------------------------------------% */
/*     | Storage Declarations:                                 | */
/*     |                                                       | */
/*     | The maximum dimensions for all arrays are             | */
/*     | set here to accommodate a problem size of             | */
/*     | N .le. MAXN                                           | */
/*     |                                                       | */
/*     | NEV is the number of eigenvalues requested.           | */
/*     |     See specifications for ARPACK usage below.        | */
/*     |                                                       | */
/*     | NCV is the largest number of basis vectors that will  | */
/*     |     be used in the Implicitly Restarted Arnoldi       | */
/*     |     Process.  Work per major iteration is             | */
/*     |     proportional to N*NCV*NCV.                        | */
/*     |                                                       | */
/*     | You must set:                                         | */
/*     |                                                       | */
/*     | MAXN:   Maximum dimension of the A allowed. (dynamic) | */
/*     | MAXNEV: Maximum NEV allowed. (dynamic)                | */
/*     | MAXNCV: Maximum NCV allowed. (dynamic)                | */
/*     %-------------------------------------------------------% */

/*     %--------------------------------------% */
/*     | F90 Allocatable Arrays (on the heap) | */
/*     %--------------------------------------% */

/*     Double precision,allocatable,save :: v(:,:) */
/*     integer,save :: v_row_allocated = 0, v_col_allocated = 0 */

/*     %----------------------------------------------% */
/*     | Originally, as F77 parameters, the following | */
/*     | integers were used to dimension work arrays. | */
/*     | They are replaced by dummy arguments used to | */
/*     | dimension the work arrays as F90 automatic   | */
/*     | arrays, but the integers are still used for  | */
/*     | passing the dimensions to lower level ARPACK | */
/*     | routines dsaupd, dseupd and dmout.           | */
/*     %----------------------------------------------% */


/* \SCCS Information: @(#) */
/* FILE: debug.h   SID: 2.3   DATE OF SID: 11/16/95   RELEASE: 2 */

/*     %---------------------------------% */
/*     | See debug.doc for documentation | */
/*     %---------------------------------% */

/*     %-------------------------------------------% */
/*     | Local F90 Automatic Arrays (on the stack) | */
/*     %-------------------------------------------% */

/*    &                 workl(ncv_in*(ncv_in+8)), */
/*    &                 workd(3*n_dim), d(ncv_in,2), resid(n_dim), */
/*    &                 ax(n_dim), */
/*     logical          select(ncv_in) */

/*     %---------------% */
/*     | Local Scalars | */
/*     %---------------% */

/*     integer          v_row_needed, v_col_needed */

/*     %------------% */
/*     | Parameters | */
/*     %------------% */


/*     %-----------------------------% */
/*     | BLAS & LAPACK routines used | */
/*     %-----------------------------% */


/*     %--------------------% */
/*     | Intrinsic function | */
/*     %--------------------% */

    /* Parameter adjustments */
    --grad;
    --xyz;
    --ax;
    --resid;
    --workd;
    --eigvecs;
    --eigvals;
    --select;
    d_dim1 = *ncv_in__;
    d_offset = 1 + d_dim1;
    d__ -= d_offset;
    --workl;
    v_dim1 = *n_dim__;
    v_offset = 1 + v_dim1;
    v -= v_offset;

    /* Function Body */

/*     %-----------------------% */
/*     | Executable Statements | */
/*     %-----------------------% */

    if (*label == 0) {
	goto L1;
    }
    switch (*label) {
	case 1:  goto L12;
	case 2:  goto L18;
    }
L1:

/*     %------------------------------------------------% */
/*     | Values used to calculate work array dimensions | */
/*     %------------------------------------------------% */

    maxn = *n_dim__;
    maxnev = *n_eig_in__;
    maxncv = *ncv_in__;
    ldv = maxn;

/*     %---------------------------------------------------% */
/*     | The include arpack_debug.h statement above and    | */
/*     | assignments here initiate trace output from the   | */
/*     | internal actions of ARPACK.  See debug.doc in the | */
/*     | DOCUMENTS directory for usage.  Initially, the    | */
/*     | most useful information will be a breakdown of    | */
/*     | time spent in the various stages of computation   | */
/*     | given by setting msaupd = 1.                      | */
/*     %---------------------------------------------------% */

    debug_1.ndigit = -5;
    debug_1.logfil = 6;
    debug_1.msgets = 0;
    debug_1.msaitr = 0;
    debug_1.msapps = 0;
    if (*debug_arpack__ == 1) {
	debug_1.msaupd = 1;
	debug_1.msaup2 = 1;
    } else {
	debug_1.msaupd = 0;
	debug_1.msaup2 = 0;
    }
/*      msaup2 = 0 */
    debug_1.mseigt = 0;
    debug_1.mseupd = 0;

/*   *** Allocatable array v will be allowed to grow to its largest size; */
/*   ***  it is never deallocated: */
/*     v_row_needed = n_dim        !!! ldv */
/*     v_col_needed = ncv_in       !!! maxncv */
/*     if( allocated(v) )then */
/*       if( (v_row_needed .gt. v_row_allocated) */
/*    & .or. (v_col_needed .gt. v_col_allocated) )then */
/*         deallocate(v,stat=ierr) */
/*         if( ierr .ne. 0 )then */
/*           write( logfil, '(a,i16,1x,i8)' ) */
/*    &       'ARPACK: could not deallocate v' */
/*           go to 9000 */
/*         endif */
/*       endif */
/*     endif */
/*     if( .not. allocated(v) )then */
/*       allocate( v(v_row_needed,v_col_needed), stat=ierr ) */
/*       if( ierr .ne. 0 )then */
/*         write( logfil, '(a,2i10)' ) */
/*    &     'ARPACK: could not allocate v' */
/*         go to 9000 */
/*       endif */
/*       v_row_allocated = v_row_needed */
/*       v_col_allocated = v_col_needed */
/*     endif */
/*     v = zero !!! zero out entire v array */

/*     %-------------------------------------------------% */
/*     | The following sets dimensions for this problem. | */
/*     %-------------------------------------------------% */

    n = *n_dim__;

/*     %----------------------------------------------% */
/*     |                                              | */
/*     | Specifications for ARPACK usage are set      | */
/*     | below:                                       | */
/*     |                                              | */
/*     |    1) NEV = N_EIG_IN  asks for N_EIG_IN      | */
/*     |       eigenvalues to be computed.            | */
/*     |                                              | */
/*     |    2) NCV = NCV_IN sets the length of the    | */
/*     |       Arnoldi factorization                  | */
/*     |                                              | */
/*     |    3) This is a standard problem             | */
/*     |         (indicated by bmat  = 'I')           | */
/*     |                                              | */
/*     |    4) Ask for the NEV eigenvalues of         | */
/*     |       smallest magnitude                     | */
/*     |         (indicated by which = 'SM')          | */
/*     |       See documentation in SSAUPD for the    | */
/*     |       other options SA, LA, LM, BE.          | */
/*     |                                              | */
/*     | Note: NEV and NCV must satisfy the following | */
/*     | conditions:                                  | */
/*     |              NEV <= MAXNEV                   | */
/*     |          NEV + 1 <= NCV <= MAXNCV            | */
/*     %----------------------------------------------% */

    nev = *n_eig_in__;
    ncv = *ncv_in__;
    *(unsigned char *)bmat = 'I';
    if (*spectrum == 1) {
	s_copy(which, "SM", (ftnlen)2, (ftnlen)2);
    } else if (*spectrum == 2) {
	s_copy(which, "SA", (ftnlen)2, (ftnlen)2);
    } else if (*spectrum == 3) {
	s_copy(which, "LM", (ftnlen)2, (ftnlen)2);
    } else if (*spectrum == 4) {
	s_copy(which, "LA", (ftnlen)2, (ftnlen)2);
    } else if (*spectrum == 5) {
	s_copy(which, "BE", (ftnlen)2, (ftnlen)2);
    } else {
	s_wsle(&io___13);
	do_lio(&c__9, &c__1, " ERROR with _SSIMP: Spectrum .NE. (SM|SA|LA|LM"
		"|BE)", (ftnlen)50);
	e_wsle();
	goto L9000;
    }

    if (n > maxn) {
	s_wsle(&io___14);
	do_lio(&c__9, &c__1, " ERROR with _SSIMP: N is greater than MAXN ", (
		ftnlen)43);
	e_wsle();
	goto L9000;
    } else if (nev > maxnev) {
	s_wsle(&io___15);
	do_lio(&c__9, &c__1, " ERROR with _SSIMP: NEV is greater than MAXNEV "
		, (ftnlen)47);
	e_wsle();
	goto L9000;
    } else if (ncv > maxncv) {
	s_wsle(&io___16);
	do_lio(&c__9, &c__1, " ERROR with _SSIMP: NCV is greater than MAXNCV "
		, (ftnlen)47);
	e_wsle();
	goto L9000;
    }

/*     %-----------------------------------------------------% */
/*     |                                                     | */
/*     | Specification of stopping rules and initial         | */
/*     | conditions before calling DSAUPD                    | */
/*     |                                                     | */
/*     | TOL  determines the stopping criterion.             | */
/*     |                                                     | */
/*     |      Expect                                         | */
/*     |           abs(lambdaC - lambdaT) < TOL*abs(lambdaC) | */
/*     |               computed   true                       | */
/*     |                                                     | */
/*     |      If TOL .le. 0,  then TOL <- macheps            | */
/*     |           (machine precision) is used.              | */
/*     |                                                     | */
/*     | IDO  is the REVERSE COMMUNICATION parameter         | */
/*     |      used to specify actions to be taken on return  | */
/*     |      from DSAUPD. (See usage below.)                | */
/*     |                                                     | */
/*     |      It MUST initially be set to 0 before the first | */
/*     |      call to DSAUPD.                                | */
/*     |                                                     | */
/*     | INFO on entry specifies starting vector information | */
/*     |      and on return indicates error codes            | */
/*     |                                                     | */
/*     |      Initially, setting INFO=0 indicates that a     | */
/*     |      random starting vector is requested to         | */
/*     |      start the ARNOLDI iteration.  Setting INFO to  | */
/*     |      a nonzero value on the initial call is used    | */
/*     |      if you want to specify your own starting       | */
/*     |      vector (This vector must be placed in RESID.)  | */
/*     |                                                     | */
/*     | The work array WORKL is used in DSAUPD as           | */
/*     | workspace.  Its dimension LWORKL is set as          | */
/*     | illustrated below.                                  | */
/*     |                                                     | */
/*     %-----------------------------------------------------% */

    lworkl = ncv * (ncv + 8);
    tol = *eigval_tol__;
    info = 0;
    ido = 0;

/*     %---------------------------------------------------% */
/*     | Specification of Algorithm Mode:                  | */
/*     |                                                   | */
/*     | This program uses the exact shift strategy        | */
/*     | (indicated by setting PARAM(1) = 1).              | */
/*     | IPARAM(3) specifies the maximum number of Arnoldi | */
/*     | iterations allowed.  Mode 1 of DSAUPD is used     | */
/*     | (IPARAM(7) = 1). All these options can be changed | */
/*     | by the user. For details see the documentation in | */
/*     | DSAUPD.                                           | */
/*     %---------------------------------------------------% */

    ishfts = 1;
    maxitr = *itr_in__;
    mode1 = 1;

    iparam[0] = ishfts;

    iparam[2] = maxitr;

    iparam[6] = mode1;

/*     %------------------------------------------------% */
/*     | M A I N   L O O P (Reverse communication loop) | */
/*     %------------------------------------------------% */

L10:

/*        %---------------------------------------------% */
/*        | Repeatedly call the routine DSAUPD and take | */
/*        | actions indicated by parameter IDO until    | */
/*        | either convergence is indicated or maxitr   | */
/*        | has been exceeded.                          | */
/*        %---------------------------------------------% */

    dsaupd_(&ido, bmat, &n, which, &nev, &tol, &resid[1], &ncv, &v[v_offset], 
	    &ldv, iparam, ipntr, &workd[1], &workl[1], &lworkl, &info, (
	    ftnlen)1, (ftnlen)2);

    if (ido == -1 || ido == 1) {

/*           %--------------------------------------% */
/*           | Perform matrix vector multiplication | */
/*           |              y <--- OP*x             | */
/*           | The user should supply his/her own   | */
/*           | matrix vector multiplication routine | */
/*           | here that takes workd(ipntr(1)) as   | */
/*           | the input, and return the result to  | */
/*           | workd(ipntr(2)).                     | */
/*           %--------------------------------------% */

	status_flag__ = 0;
L11:
	hessvec_(&n, &workd[ipntr[0]], &workd[ipntr[1]], &xyz[1], &grad[1], 
		return_flag__, &status_flag__);
	if (status_flag__ == 0) {
	    goto L13;
	}
	if (status_flag__ < 0) {
	    goto L9000;
	}
	*label = l12;
	return 0;
L12:
	goto L11;
L13:

/*           %-----------------------------------------% */
/*           | L O O P   B A C K to call DSAUPD again. | */
/*           %-----------------------------------------% */

	goto L10;

    }

/*     %----------------------------------------% */
/*     | Either we have convergence or there is | */
/*     | an error.                              | */
/*     %----------------------------------------% */

    if (info < 0) {

/*        %--------------------------% */
/*        | Error message. Check the | */
/*        | documentation in DSAUPD. | */
/*        %--------------------------% */

	s_wsle(&io___27);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	s_wsle(&io___28);
	do_lio(&c__9, &c__1, " Error with _saupd, info = ", (ftnlen)27);
	do_lio(&c__3, &c__1, (char *)&info, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___29);
	do_lio(&c__9, &c__1, " Check documentation in _saupd ", (ftnlen)31);
	e_wsle();
	s_wsle(&io___30);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	goto L9000;

    } else {

/*        %-------------------------------------------% */
/*        | No fatal errors occurred.                 | */
/*        | Post-Process using DSEUPD.                | */
/*        |                                           | */
/*        | Computed eigenvalues may be extracted.    | */
/*        |                                           | */
/*        | Eigenvectors may be also computed now if  | */
/*        | desired.  (indicated by rvec = .true.)    | */
/*        |                                           | */
/*        | The routine DSEUPD now called to do this  | */
/*        | post processing (Other modes may require  | */
/*        | more complicated post processing than     | */
/*        | mode1.)                                   | */
/*        |                                           | */
/*        %-------------------------------------------% */

	if (*need_eigvecs__ == 1) {
	    rvec = TRUE_;
	} else {
	    rvec = FALSE_;
	}

	dseupd_(&rvec, "All", &select[1], &d__[d_offset], &v[v_offset], &ldv, 
		&sigma, bmat, &n, which, &nev, &tol, &resid[1], &ncv, &v[
		v_offset], &ldv, iparam, ipntr, &workd[1], &workl[1], &lworkl,
		 ierr, (ftnlen)3, (ftnlen)1, (ftnlen)2);

/*        %----------------------------------------------% */
/*        | Eigenvalues are returned in the first column | */
/*        | of the two dimensional array D and the       | */
/*        | corresponding eigenvectors are returned in   | */
/*        | the first NCONV (=IPARAM(5)) columns of the  | */
/*        | two dimensional array V if requested.        | */
/*        | Otherwise, an orthogonal basis for the       | */
/*        | invariant subspace corresponding to the      | */
/*        | eigenvalues in D is returned in V.           | */
/*        %----------------------------------------------% */

	if (*ierr != 0) {

/*           %------------------------------------% */
/*           | Error condition:                   | */
/*           | Check the documentation of DSEUPD. | */
/*           %------------------------------------% */

	    s_wsle(&io___33);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    s_wsle(&io___34);
	    do_lio(&c__9, &c__1, " Error with _seupd, info = ", (ftnlen)27);
	    do_lio(&c__3, &c__1, (char *)&(*ierr), (ftnlen)sizeof(integer));
	    e_wsle();
	    s_wsle(&io___35);
	    do_lio(&c__9, &c__1, " Check the documentation of _seupd. ", (
		    ftnlen)36);
	    e_wsle();
	    s_wsle(&io___36);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    goto L9000;

	} else if (*debug_arpack__ == 1) {

	    nconv = iparam[4];
	    *n_eig_out__ = nconv;
	    if (nconv <= 0) {
		s_wsle(&io___38);
		do_lio(&c__9, &c__1, " ", (ftnlen)1);
		e_wsle();
		s_wsle(&io___39);
		do_lio(&c__9, &c__1, " ARPACK: Not a single mode converged.", 
			(ftnlen)37);
		e_wsle();
		s_wsle(&io___40);
		do_lio(&c__9, &c__1, " ", (ftnlen)1);
		e_wsle();
		goto L9000;
	    }

/*           %--------------------------------------------% */
/*           | "UnDO" DO 20 j=1,nconv loop, because it is | */
/*           | illegal to jump in and out from a DO loop. | */
/*           %--------------------------------------------% */

	    j = 1;
L16:

/*              %---------------------------% */
/*              | Compute the residual norm | */
/*              |                           | */
/*              |   ||  A*x - lambda*x ||   | */
/*              |                           | */
/*              | for the NCONV accurately  | */
/*              | computed eigenvalues and  | */
/*              | eigenvectors.  (iparam(5) | */
/*              | indicates how many are    | */
/*              | accurate to the requested | */
/*              | tolerance)                | */
/*              %---------------------------% */

	    status_flag__ = 0;
L17:
	    hessvec_(&n, &v[j * v_dim1 + 1], &ax[1], &xyz[1], &grad[1], 
		    return_flag__, &status_flag__);
	    if (status_flag__ == 0) {
		goto L19;
	    }
	    if (status_flag__ < 0) {
		goto L9000;
	    }
	    *label = l18;
	    return 0;
L18:
	    goto L17;
L19:

	    d__1 = -d__[j + d_dim1];
	    daxpy_(&n, &d__1, &v[j * v_dim1 + 1], &c__1, &ax[1], &c__1);
	    d__[j + (d_dim1 << 1)] = dnrm2_(&n, &ax[1], &c__1);
	    d__[j + (d_dim1 << 1)] /= (d__1 = d__[j + d_dim1], abs(d__1));

	    ++j;
	    if (j > nconv) {
		goto L20;
	    }

	    goto L16;

L20:

/*           %-----------------------------% */
/*           | Display computed residuals. | */
/*           %-----------------------------% */

#if 0
	    dmout_(&c__6, &nconv, &c__2, &d__[d_offset], &maxncv, &c_n6, 
		    "Ritz values and relative residuals", (ftnlen)34);
#endif

/*           %-------------------------------------------% */
/*           | Print additional convergence information. | */
/*           %-------------------------------------------% */

	    if (info == 1) {
		s_wsle(&io___42);
		do_lio(&c__9, &c__1, " ", (ftnlen)1);
		e_wsle();
		s_wsle(&io___43);
		do_lio(&c__9, &c__1, " Maximum number of iterations reached.",
			 (ftnlen)38);
		e_wsle();
		s_wsle(&io___44);
		do_lio(&c__9, &c__1, " ", (ftnlen)1);
		e_wsle();
	    } else if (info == 3) {
		s_wsle(&io___45);
		do_lio(&c__9, &c__1, " ", (ftnlen)1);
		e_wsle();
		s_wsle(&io___46);
		do_lio(&c__9, &c__1, " No shifts could be applied during imp"
			"licit", (ftnlen)43);
		do_lio(&c__9, &c__1, " Arnoldi update, try increasing NCV.", (
			ftnlen)36);
		e_wsle();
		s_wsle(&io___47);
		do_lio(&c__9, &c__1, " ", (ftnlen)1);
		e_wsle();
	    }

	    s_wsle(&io___48);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    s_wsle(&io___49);
	    do_lio(&c__9, &c__1, " _SSIMP ", (ftnlen)8);
	    e_wsle();
	    s_wsle(&io___50);
	    do_lio(&c__9, &c__1, " ====== ", (ftnlen)8);
	    e_wsle();
	    s_wsle(&io___51);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    s_wsle(&io___52);
	    do_lio(&c__9, &c__1, " Size of the matrix is ", (ftnlen)23);
	    do_lio(&c__3, &c__1, (char *)&n, (ftnlen)sizeof(integer));
	    e_wsle();
	    s_wsle(&io___53);
	    do_lio(&c__9, &c__1, " The number of Ritz values requested is ", (
		    ftnlen)40);
	    do_lio(&c__3, &c__1, (char *)&nev, (ftnlen)sizeof(integer));
	    e_wsle();
	    s_wsle(&io___54);
	    do_lio(&c__9, &c__1, " The number of Arnoldi vectors generated", (
		    ftnlen)40);
	    do_lio(&c__9, &c__1, " (NCV) is ", (ftnlen)10);
	    do_lio(&c__3, &c__1, (char *)&ncv, (ftnlen)sizeof(integer));
	    e_wsle();
	    s_wsle(&io___55);
	    do_lio(&c__9, &c__1, " What portion of the spectrum: ", (ftnlen)
		    31);
	    do_lio(&c__9, &c__1, which, (ftnlen)2);
	    e_wsle();
	    s_wsle(&io___56);
	    do_lio(&c__9, &c__1, " The number of converged Ritz values is ", (
		    ftnlen)40);
	    do_lio(&c__3, &c__1, (char *)&nconv, (ftnlen)sizeof(integer));
	    e_wsle();
	    s_wsle(&io___57);
	    do_lio(&c__9, &c__1, " The number of Implicit Arnoldi update", (
		    ftnlen)38);
	    do_lio(&c__9, &c__1, " iterations taken is ", (ftnlen)21);
	    do_lio(&c__3, &c__1, (char *)&iparam[2], (ftnlen)sizeof(integer));
	    e_wsle();
	    s_wsle(&io___58);
	    do_lio(&c__9, &c__1, " The number of OP*x is ", (ftnlen)23);
	    do_lio(&c__3, &c__1, (char *)&iparam[8], (ftnlen)sizeof(integer));
	    e_wsle();
	    s_wsle(&io___59);
	    do_lio(&c__9, &c__1, " The convergence criterion is ", (ftnlen)30)
		    ;
	    do_lio(&c__5, &c__1, (char *)&tol, (ftnlen)sizeof(doublereal));
	    e_wsle();
	    s_wsle(&io___60);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	}

/*        %----------------------------% */
/*        | Return eigvals and eigvecs | */
/*        %----------------------------% */

	nconv = iparam[4];
	*n_eig_out__ = nconv;
	if (nconv <= 0) {
	    s_wsle(&io___61);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    s_wsle(&io___62);
	    do_lio(&c__9, &c__1, " ARPACK: Not a single mode converged.", (
		    ftnlen)37);
	    e_wsle();
	    s_wsle(&io___63);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    goto L9000;
	}

	i__1 = nconv;
	for (j = 1; j <= i__1; ++j) {
	    eigvals[j] = d__[j + d_dim1];

	    i__2 = n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		eigvecs[(j - 1) * n + i__] = v[i__ + j * v_dim1];
/* L30: */
	    }
/* L40: */
	}

    }

/*     %--------------------------------% */
/*     | Done with subroutine dsarpack. | */
/*     %--------------------------------% */

    *label = 0;
    return 0;

L9000:

/* !! Error */
    if (status_flag__ == 0) {
	status_flag__ = arpack_error__;
    }

    *label = status_flag__;
    return 0;

} /* dsarpack_ */

