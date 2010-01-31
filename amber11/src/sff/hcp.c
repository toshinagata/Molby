
/***********************************************************************
 *                          calc_dist2()
 *
 * Calculate distance square
 *
 * Calling Parameters:
 * xi, yi, zi - point i
 * xj, yj, zj - point j
 * 
 * return: dist2
 ************************************************************************/
REAL_T calc_dist2(REAL_T xi, REAL_T yi, REAL_T zi, REAL_T xj, REAL_T yj, REAL_T zj)
{
   REAL_T dist2, xij, yij, zij;

   xij = xi - xj;
   yij = yi - yj;
   zij = zi - zj;
   dist2 = xij * xij + yij * yij + zij * zij;

   return(dist2);
}

/***********************************************************************
 *                          calc_approx_q()
 *
 * Calculate approximate charges and geometric center of components
 *
 * Calling Parameters:
 * x - atom coordinates
 * x_hcp1 - geometric centers of residues
 * x_hcp2 - geometric centers of strands
 * q_hcp1 - approximate charges for residues
 * q_hcp2 - approximate charges for strands
 ************************************************************************/
calc_approx_q(REAL_T * x, REAL_T * x_hcp1, REAL_T * x_hcp2, 
   REAL_T * q_hcp1, REAL_T * q_hcp2, INT_T hcp, PARMSTRUCT_T * prm)
{

   int s, r, a, m;                    /* index for strand, residue, atom, charge */
   int a_from, a_to;                  /* for atoms within a residue */
   int natoms_hcp1, natoms_hcp2;      /* total number of atoms */
   int nstrands = 1;                  /* UNTIL I CAN ADD STRUCTURE FOR STRANDS */
   REAL_T q;                          /* charge */
   REAL_T min_q = 0.001;              /* to prevent divide by zero */

   for (s = 0; s < nstrands; s++) /* for each strand */
   {
      /* initialize strands */
      natoms_hcp2 = 0;
      x_hcp2[3*s+0] = x_hcp2[3*s+1] = x_hcp2[3*s+2] = 0.0;

      for (m = 0; m < hcp; m++) 
      { 
          q_hcp2[4*(s*hcp+m)+0] = q_hcp2[4*(s*hcp+m)+1] = 
             q_hcp2[4*(s*hcp+m)+2] = q_hcp2[4*(s*hcp+m)+3] = 0.0;
      }

      for (r = 0; r < prm->Nres; r++)         /* for each residue */
      {
         /* initialize residues */
         natoms_hcp1 = 0;
         x_hcp1[3*r+0] = x_hcp1[3*r+1] = x_hcp1[3*r+2] = 0.0;

         for (m = 0; m < hcp; m++) 
         { 
             q_hcp1[4*(r*hcp+m)+0] = q_hcp1[4*(r*hcp+m)+1] = 
                q_hcp1[4*(r*hcp+m)+2] = q_hcp1[4*(r*hcp+m)+3] = 0.0;
         }

         a_from = prm->Ipres[r] - 1;
         if (r + 1 < prm->Nres) { a_to = prm->Ipres[r + 1] - 1; }
         else { a_to = prm->Natom; }
         for (a = a_from; a < a_to; a++)    /* for each atoms */
         {
            m = 0;
            q = prm->Charges[a];
            switch (hcp)
            {
               /* cutoff - should never get here */
               case 0: 
	          break;
               /* center of 1 charge (m=0) */
               case 1: 
	          break;
               /* centers of 2 charges (pos:m=1 & neg:m=0) */
               case 2: 
                  if (q > 0.0) { m = 1; }
                  break;
               /* centers of 3 charges (m=0/1/2) */
               case 3:
                  if (q > -0.333)   
                  { 
                     m = 1; 
                     if (q > 0.333) { m = 2; }
                  }
                  break;
            } /* end switch */
            ++natoms_hcp1;
            /* add to geometric center */
            x_hcp1[3*r+0] += x[3*a+0];
            x_hcp1[3*r+1] += x[3*a+1];
            x_hcp1[3*r+2] += x[3*a+2];
            /* add to center of charge */
            q_hcp1[4*(r*hcp+m)+0] += x[3*a+0] * q;
            q_hcp1[4*(r*hcp+m)+1] += x[3*a+1] * q;
            q_hcp1[4*(r*hcp+m)+2] += x[3*a+2] * q;
            q_hcp1[4*(r*hcp+m)+3] += q;
         }  /* end for each atoms */
         /* roll up to next level */
         natoms_hcp2 += natoms_hcp1;
         x_hcp2[3*s+0] += x_hcp1[3*r+0];
         x_hcp2[3*s+1] += x_hcp1[3*r+1];
         x_hcp2[3*s+2] += x_hcp1[3*r+2];
         for (m = 0; m < hcp; m++) 
         { 
             q_hcp2[4*(s*hcp+m)+0] += q_hcp1[4*(r*hcp+m)+0];
             q_hcp2[4*(s*hcp+m)+1] += q_hcp1[4*(r*hcp+m)+1];
             q_hcp2[4*(s*hcp+m)+2] += q_hcp1[4*(r*hcp+m)+2];
             q_hcp2[4*(s*hcp+m)+3] += q_hcp1[4*(r*hcp+m)+3];
         }
         /* compute geometric center and approx charges */
         if (natoms_hcp1 > 0)
         {
            x_hcp1[3*r+0] /= natoms_hcp1;
            x_hcp1[3*r+1] /= natoms_hcp1;
            x_hcp1[3*r+2] /= natoms_hcp1;
         }
         for (m = 0; m < hcp; m++) 
         {
             q = q_hcp1[4*(r*hcp+m)+3];
             if (fabs(q) > min_q)
             {  
                q_hcp1[4*(r*hcp+m)+0] /= q;
                q_hcp1[4*(r*hcp+m)+1] /= q;
                q_hcp1[4*(r*hcp+m)+2] /= q;
             }
         }
      }  /* end for each residue */
      /* compute geometric center and approx charges */
      if (natoms_hcp2 > 0)
      {
         x_hcp2[3*s+0] /= natoms_hcp2;
         x_hcp2[3*s+1] /= natoms_hcp2;
         x_hcp2[3*s+2] /= natoms_hcp2;
      }
      for (m = 0; m < hcp; m++) 
      { 
          q = q_hcp2[4*(s+m)+3];
          if (fabs(q) > min_q)
          {  
             q_hcp2[4*(s*hcp+m)+0] /= q;
             q_hcp2[4*(s*hcp+m)+1] /= q;
             q_hcp2[4*(s*hcp+m)+2] /= q;
          }
      }
   }  /* end for each strand */

}  /* end calc_approx_q() */




/***********************************************************************
 *                           nblist_hcp()
 *
 * Build HCP non-bonded pairlist for atoms and residues.
 * 1. Recalculate approximate charges for chains and residues
 * 2. For each atom, identify chains, residues and atoms within threshod distances
 *
 * Calling parameters:
 * npairs_hcpx - (x=0/1/2) number of pairs in each pairlist
 * pairs_hcp0 - list of atom, a list of atom numbers (level 0)
 * pairs_hcp1 - for each atom, a list of residue numbers (lelvel 1)
 * pairs_hcp2 - for each atom, a list of chain numbers (lelvel 2)
 * x          - atom coordinates
 ************************************************************************/

int nblist_hcp(int * npairs_hcp0, int * npairs_hcp1, int * npairs_hcp2,
   int ** pairlist_hcp0, int **pairlist_hcp1, int **pairlist_hcp2, 
   REAL_T * x, REAL_T * x_hcp1, REAL_T * x_hcp2, REAL_T * q_hcp1, REAL_T * q_hcp2,
   int hcp, REAL_T hcp_h, PARMSTRUCT_T * prm)
{
   int s, s1, r, r1, a, a1;                /* index for strand, residue and atom */
   int a_from, a_to, a1_from, a1_to;       /* for atoms within a residue */
   int nstrands = 1;                       /* UNTIL I CAN ADD STRUCTURE FOR STRANDS */
   int *temp_hcp0 = NULL;                  /* temp pairlists copied to parlist_hcpx */
   int *temp_hcp1 = NULL;
   int *temp_hcp2 = NULL;
   REAL_T dist2, dist2_hcp1, dist2_hcp2;   /* square of threshold distances */

   /* calculate approximate charges and geometric centers */
   calc_approx_q(x, x_hcp1, x_hcp2, q_hcp1, q_hcp2, hcp, prm);
   
   temp_hcp0 = ivector(0, prm->Natom);
   temp_hcp1 = ivector(0, prm->Nres);
   temp_hcp2 = ivector(0, nstrands);

   /* square of threshold distances */
   dist2_hcp1 = hcp_h * hcp_h;
   dist2_hcp2 = 75 * 75;

   /* free up old pairlist */
   if (pairlist_hcp0 != NULL) {
      for (a = 0; a < prm->Natom; a++) {
//         if (pairlist_hcp0[a] != NULL) {
            free_ivector(pairlist_hcp0[a], 0, 1); // npairs_hcp0[a]);
//         }
      }
   }
   if (pairlist_hcp1 != NULL) {
      for (a = 0; a < prm->Natom; a++) {
//         if (pairlist_hcp1[a] != NULL) {
            free_ivector(pairlist_hcp1[a], 0, 1); // npairs_hcp1[a]);
//         }
      }
   }
   if (pairlist_hcp2 != NULL) {
      for (a = 0; a < prm->Natom; a++) {
//         if (pairlist_hcp2[a] != NULL) {
            free_ivector(pairlist_hcp2[a], 0, 1); // npairs_hcp2[a]);
//         }
      }
   }


   /* calculate pairlists */
   for (s = 0; s < nstrands; s++) /* for each strand */
   {
      for (r = 0; r < prm->Nres; r++)         /* for each residue */
      {
         a_from = prm->Ipres[r] - 1;
         if (r < prm->Nres - 1) { a_to = prm->Ipres[r + 1] - 1; }
         else { a_to = prm->Natom; }
         for (a = a_from; a < a_to; a++)    /* for each atoms */
         {
//           printf("For atom %d in res %d: \n", a, r);

            npairs_hcp0[a] = npairs_hcp1[a] = npairs_hcp2[a] = 0;
            for (s1 = 0; s1 < nstrands; s1++) /* for each other strand */
            {
               dist2 = 0;
               if (s != s1)                  /* if not in same strand calc dist to it */
               {
                  dist2 = calc_dist2(x[3*a+0], x[3*a+1], x[3*a+2],
                     x_hcp2[3*s1+0], x_hcp2[3*s1+1], x_hcp2[3*s1+2]);
               }
               if (dist2 > dist2_hcp2)
               {
                  temp_hcp2[npairs_hcp2[a]] = s1;
                  ++npairs_hcp2[a];
               }
               else
               {
                  for (r1 = 0; r1 < prm->Nres; r1++)      /* for each other residue */
                  {
                     dist2 = 0;
                     if (r != r1)                         /* if not in the same residue */
                     {
                        dist2 = calc_dist2(x[3*a+0], x[3*a+1], x[3*a+2],
                           x_hcp1[3*r1+0], x_hcp1[3*r1+1], x_hcp1[3*r1+2]);
                     }
                     if (dist2 > dist2_hcp1)
                     {
//                        printf("Res: %d, d2: %f, x1 %f, %f, %f, x2 %f, %f, %f \n", r1, dist2,
//                           x[3*a+0], x[3*a+1], x[3*a+2], x_hcp1[3*r1+0], x_hcp1[3*r1+1], x_hcp1[3*r1+2]);
//                        fflush(stdout);
                        temp_hcp1[npairs_hcp1[a]] = r1;
                        ++npairs_hcp1[a];
                     }
                     else
                     {
                        a1_from = prm->Ipres[r1] - 1;
                        if (r1 < prm->Nres - 1) { a1_to = prm->Ipres[r1 + 1] - 1; }
                        else { a1_to = prm->Natom; }
                        for (a1 = a1_from; a1 < a1_to; a1++)    /* for each atoms */
                        {
                           if (a != a1)
                           {
//                              printf("Atom: %d, npairs %d, res dist %f\n", a1, npairs_hcp0[a], dist2);
//                              fflush(stdout);
                              temp_hcp0[npairs_hcp0[a]] = a1;
                              ++npairs_hcp0[a];
                           }  /* not same atom */
                        }  /* end for other atoms */
                     }  /* end else residue inside threshold dist */
                  }  /* end for other residues */
               }  /* end else strand inside threshold dist */
            }  /* end for other strand */

            /* copy temp pairlists to global array */
//            printf("Atom %d pairlist: %d atoms, %d res, %d strands \n",
//               a, npairs_hcp0[a], npairs_hcp1[a], npairs_hcp2[a]);
//            fflush(stdout);
            if (npairs_hcp0[a] > 0)
            {
               pairlist_hcp0[a] = ivector(0, npairs_hcp0[a]);
               for (a1 = 0; a1 < npairs_hcp0[a]; a1++)
               {
                  pairlist_hcp0[a][a1] = temp_hcp0[a1];
               }
            }
            if (npairs_hcp1[a] > 0)
            {
               pairlist_hcp1[a] = ivector(0, npairs_hcp1[a]);
               for (r1 = 0; r1 < npairs_hcp1[a]; r1++)
               {
                  pairlist_hcp1[a][r1] = temp_hcp1[r1];
               }
            }
            if (npairs_hcp2[a] > 0)
            {
               pairlist_hcp2[a] = ivector(0, npairs_hcp2[a]);
               for (s1 = 0; s1 < npairs_hcp2[a]; s1++)
               {
                  pairlist_hcp2[a][s1] = temp_hcp2[s1];
               }
            } 

         }  /* end for each atom */
      }  /* end for each residue */
   }  /* end for each strand  */

   free_ivector(temp_hcp0, 0, prm->Natom);
   free_ivector(temp_hcp1, 0, prm->Nres);
   free_ivector(temp_hcp2, 0, nstrands);


}  /* end nblist_hcp() */

/***********************************************************************
                            NBOND_HCP() 
************************************************************************/

/* 
 * Calculate the non-bonded energy and first derivatives for HCP.
 * The non-bonded pair list must be modified by the excluded atom list 
 * 
 * Calling parameters are as follows:
 *
 * npairs_hcp0   - number of entries in level 0 pairlist for each atom
 * pairlist_hcp0 - the non-bonded level 0 pair list
 * N14 - set to 0 for the non-bonded pair list, 1 for the 1-4 pair list
 * x - the atomic coordinate array
 * f - the gradient vector
 * enb - Van der Waals energy return value, passed by reference
 * eel - Coulombic energy return value, passed by reference
 * enbfac - scale factor for Van der Waals energy
 * eelfac - scale factor for Coulombic energy
 */

static int nbond_hcp(int * npairs_hcp0, int * npairs_hcp1, int * npairs_hcp2,
   int ** pairlist_hcp0, int ** pairlist_hcp1, int ** pairlist_hcp2,
   int * Iblo_hcp, int ** IexclAt_hcp, 
   REAL_T * x, REAL_T * q_hcp1, REAL_T * q_hcp2, REAL_T * f, REAL_T * enb, REAL_T * eel)
{
   int i, j, k, jj, jn, ic, npr, lpair, iaci, foff, threadnum, numthreads;	
   int *iexw;
   REAL_T dumx, dumy, dumz, dumw, cgi, cgj, r2inv, df2, r6, r10, f1, f2;
   REAL_T dedx, dedy, dedz, dedw, df, enbfaci, eelfaci, evdw, elec;
   REAL_T xi, yi, zi, wi, xij, yij, zij, wij, r, r2;
   REAL_T dis, kij, d0, diff, rinv, rs, rssq, eps1, epsi, cgijr, pow;
   int ibig, isml;
   REAL_T enbfac, eelfac; /* CAN REMOVE WITH MINOR CHANGES ! */

   REAL_T QMIN = 0.00001;

#define SIG 0.3
#define DIW 78.0
#define C1 38.5

   evdw = 0.0;
   elec = 0.0;
   enbfac = 1.0; /* SCALING FACTOR ONLY FOR N14 INTERACTIONS - REMOVE! */
   eelfac = 1.0; /* SCALING FACTOR ONLY FOR N14 INTERACTIONS - REMOVE! */
   enbfaci = 1.0 / enbfac;
   eelfaci = 1.0 / eelfac;
	

#pragma omp parallel reduction (+: evdw, elec) \
  private (i, j, iexw, npr, iaci, \
           xi, yi, zi, wi, xij, yij, zij, wij, dumx, dumy, dumz, dumw, \
           cgi, jn, r2, r2inv, r, rinv, rs, rssq, pow, \
           eps1, epsi, cgijr, df2, ic, r6, f2, f1, df, dis, d0, kij, \
           diff, ibig, isml, dedx, dedy, dedz, dedw, r10, \
           threadnum, numthreads, foff, lpair)
   	{
      /*
       * Get the thread number and the number of threads for multi-threaded
       * execution under OpenMP.  For all other cases, including ScaLAPACK,
       * MPI and single-threaded execution, use the values that have been
       * stored in mytaskid and numtasks, respectively.
       */

#if defined(OPENMP)
      threadnum = omp_get_thread_num();
      numthreads = omp_get_num_threads();
#else
      threadnum = mytaskid;
      numthreads = numtasks;
#endif

      /*
       * Compute an offset into the gradient array for this thread,
       * but only if OPENMP is defined.
       */

#ifdef OPENMP
      foff = dim * prm->Natom * threadnum;
#else
      foff = 0;
#endif

      /* iexw stores the list of excluded atoms of each atom */
      iexw = ivector(-1, prm->Natom+1);
      for (i = -1; i < prm->Natom+1; i++) {
         iexw[i] = -1;
      }

      /*
       * Loop over all atoms i.
       *
       * Explicitly assign threads to loop indices for the following loop,
       * in a manner equivalent to (static, N) scheduling with OpenMP, and
       * identical to the manner in which threads are assigned in nblist.
       *
       * Synchronization of OpenMP threads will occur following this loop
       * because the parallel region ends after this loop.  Following
       * synchronization, a reduction of the sumdeijda array will be
       * performed.
       *
       * Synchronization of MPI tasks will occur via the MPI_Allreduce
       * function that is called from within mme34.
       */
		 		
   for (i = 0; i < prm->Natom; i++) {
	#if defined(OPENMP) || defined(MPI) || defined(SCALAPACK)
         	if (!myroc(i, BLOCKSIZE, numthreads, threadnum))
            	continue;
	#endif

         iaci = prm->Ntypes * (prm->Iac[i] - 1);
         dumx = dumy = dumz = 0.0;
         xi = x[dim * i + 0];
         yi = x[dim * i + 1];
         zi = x[dim * i + 2];

         if (dim == 4) {
            dumw = 0.0;
            wi = x[dim * i + 3];
         }

         cgi = eelfaci * prm->Charges[i];
//            printf("nbond: i %d x %f, %f, %f q %f\n", 
//               i, xi, yi, zi, cgi);
//            fflush(stdout);

         /*
          * Expand the excluded list into the iexw array 
          * by storing i at array address j
          * THERE'S GOT TO BE A MORE INTUITIVE WAY TO DO THIS!
          */
         for (j = 0; j < Iblo_hcp[i]; j++) {
            iexw[IexclAt_hcp[i][j] - 1] = i;
         } 

         for (jj = 0; jj < npairs_hcp0[i]; jj++)
         {
            j = pairlist_hcp0[i][jj]; 	    
            /*
             * If the 'N14' calling parameter is clear, check whether
             * this i,j pair is exempted by the excluded atom list.
             */

//            if (j>=0 && j!=i && (N14 != 0 || iexw[j] != i) && mlca_ptr->type=='a'){
//               xij = xi - mlca_ptr->pos[0];
//               yij = yi - mlca_ptr->pos[1];
//               zij = zi - mlca_ptr->pos[2];
//               printf("atom j %d, npairs %d, jj %d \n", j, npairs_hcp0[i], jj);
//               fflush(stdout);
            if (iexw[j] != i)
            {
//                  printf("xj %f, %f, %f, qj %f \n", 
//                     x[dim*j], x[dim*j+1], x[dim*j+2], prm->Charges[j]);
//                  fflush(stdout);
               xij = xi - x[dim * j + 0];
               yij = yi - x[dim * j + 1];
               zij = zi - x[dim * j + 2];

/*             FROM HERE TO END OF IF MOVE TO A SUBROUTINE! */

               r2 = xij * xij + yij * yij + zij * zij;
               
               if (dim == 4) {
                  wij = wi - x[dim * j + 3];
                  r2 += wij * wij;
               }

               r2inv = 1.0 / r2;
               r = sqrt(r2);
               rinv = r * r2inv;

               /* Calculate the energy and derivatives according to dield. */

               if (dield == -3) {

                  /* special code Ramstein & Lavery dielectric, 94 force field */

                  rs = SIG * r;
                  rssq = rs * rs;
                  pow = exp(-rs);
                  eps1 = rssq + rs + rs + 2.0;
                  epsi = 1.0 / (DIW - C1 * pow * eps1);
                  cgijr = cgi * prm->Charges[j] * rinv * epsi;
                  elec += cgijr;
                  df2 = -cgijr * (1.0 + C1 * pow * rs * rssq * epsi);
                  ic = prm->Cno[iaci + prm->Iac[j] - 1] - 1;
                  if (ic >= 0) {
                     r6 = r2inv * r2inv * r2inv;
                     f2 = prm->Cn2[ic] * r6;
                     f1 = prm->Cn1[ic] * r6 * r6;
                     evdw += (f1 - f2) * enbfaci;
                     df = (df2 + (6.0 * f2 - 12.0 * f1) * enbfaci) * rinv;
                  } else {
                     df = df2 * rinv;
                  }

               } else if (dield == -4) {

                  /* distance-dependent dielectric code, 94 ff */
                  /* epsilon = r  */

                  rs = cgi * prm->Charges[j] * r2inv;
                  df2 = -2.0 * rs;
                  elec += rs;
                  ic = prm->Cno[iaci + prm->Iac[j] - 1] - 1;
                  if (ic >= 0) {
                     r6 = r2inv * r2inv * r2inv;
                     f2 = prm->Cn2[ic] * r6;
                     f1 = prm->Cn1[ic] * r6 * r6;
                     evdw += (f1 - f2) * enbfaci;
                     df = (df2 + (6.0 * f2 - 12.0 * f1) * enbfaci) * rinv;
                  } else {
                     df = df2 * rinv;
                  }

               } else if (dield == -5) {

                  /* non-bonded term from yammp  */

                  dis = r;
                  ic = prm->Cno[iaci + prm->Iac[j] - 1] - 1;
                  d0 = prm->Cn2[ic];
                  if (dis < d0) {
                     kij = prm->Cn1[ic];
                     diff = dis - d0;
                     evdw += kij * diff * diff;
                     df = 2.0 * kij * diff;
                  } else {
                     df = 0.0;
                  }
               } else {

                  /*
                   * Code for various dielectric models.
                   * The df2 variable should hold r(dV/dr).
                   */

                  if (dield == 0) {

                     /* epsilon = r  */

                     rs = cgi * prm->Charges[j] * r2inv;
                     df2 = -2.0 * rs;
                     elec += rs;

                  } else if (dield == 1) {

                     /* epsilon = 1  */
                     rs = (cgi * prm->Charges[j] * rinv);
                     df2 = -rs;
                     elec += rs;
			
                  } else if (dield == -2) {

                     /* Ramstein & Lavery dielectric, PNAS 85, 7231 (1988). */

                     rs = SIG * r;
                     rssq = rs * rs; 
                     pow = exp(-rs);
                     eps1 = rssq + rs + rs + 2.0;
                     epsi = 1.0 / (DIW - C1 * pow * eps1);
                     cgijr = cgi * prm->Charges[j] * rinv * epsi;
                     elec += cgijr;
                     df2 = -cgijr * (1.0 + C1 * pow * rs * rssq * epsi);
                  }

                  /* Calculate either Van der Waals or hydrogen bonded term. */

                  ic =prm->Cno[iaci + prm->Iac[j] - 1];
                  if (ic > 0 || enbfac != 1.0) 
		  {
			
                     if (ic > 0) 
		     {
                        ic--;
                     } 
		     else 
		     {
                        ibig = prm->Iac[i] > prm->Iac[j] ?
                            prm->Iac[i] : prm->Iac[j];
                        isml = prm->Iac[i] > prm->Iac[j] ?
                            prm->Iac[j] : prm->Iac[i];
                        ic = ibig * (ibig - 1) / 2 + isml - 1;
                     }
                     r6 = r2inv * r2inv * r2inv;
                     f2 = prm->Cn2[ic] * r6;
                     f1 = prm->Cn1[ic] * r6 * r6;
                     evdw += (f1 - f2) * enbfaci;
                     df = (df2 + (6.0 * f2 - 12.0 * f1) * enbfaci) * rinv;
#if 0
                     if (enbfac != 1.0)
                        nb14 += (f1 - f2) * enbfaci;
#endif
                  } 
		  else {
                     ic = -ic - 1;
                     r10 = r2inv * r2inv * r2inv * r2inv * r2inv;
                     f2 = prm->HB10[ic] * r10;
                     f1 = prm->HB12[ic] * r10 * r2inv;
                     evdw += (f1 - f2) * enbfaci;
                     df = (df2 + (10.0 * f2 - 12.0 * f1) * enbfaci) * rinv;
#if 0
                     hbener += (f1 - f2) * enbfaci;
#endif
                  }
//                  df = df2 * rinv; /* Ramu ignore vdw */
               }

               /*
                * The de term contains one more factor of Dij in the denominator
                * so that terms such as dedx do not need to include 1/Dij. 
                *
                * Update the gradient for atom j.
                */

               df *= rinv;

               dedx = df * xij;
               dedy = df * yij;
               dedz = df * zij;

               dumx += dedx;
               dumy += dedy;
               dumz += dedz;

/*   Ramu: not needed because force is calculated both ways 1-j and j-i */ 
/* 
               f[foff + dim * j + 0] -= dedx;
               f[foff + dim * j + 1] -= dedy;
               f[foff + dim * j + 2] -= dedz;
		
               if (dim == 4) {
                  dedw = df * wij;
                  dumw += dedw;
                  f[foff + dim * j + 3] -= dedw;
               }
*/
            }  /* end of if not on exlcuded atom list */
         }  /* end of for loop for pairlist_hcp0 */

//         printf("end of atom pairlist, starting residue pairlist\n");
//         fflush(stdout);

         /* calculate energy/force for level 1 approximations */
         for (jj = 0; jj < npairs_hcp1[i]; jj++)
         {
            j = pairlist_hcp1[i][jj];

            for (k = 0; k < hcp; k++)
            {
//               printf("res %d, charge %d, xj %f, %f, %f, qj %f\n", k, k
//                  q_hcp1[4*(j*hcp+k)+0], q_hcp1[4*(j*hcp+k)+1],
//                  q_hcp1[4*(j*hcp+k)+2], q_hcp1[4*(j*hcp+k)+3]);
//               fflush(stdout);

               cgj = q_hcp1[4*(j*hcp+k)+3];
               if (fabs(cgj) > QMIN)
               { 
                  xij = xi - q_hcp1[4*(j*hcp+k)+0];
                  yij = yi - q_hcp1[4*(j*hcp+k)+1];
                  zij = zi - q_hcp1[4*(j*hcp+k)+2];
                  r2 = xij * xij + yij * yij + zij * zij;

                  if (dim == 4) {
                     wij = wi - x[dim * j + 3];
                     r2 += wij * wij;
                  }

                  r2inv = 1.0 / r2;
                  r = sqrt(r2);
                  rinv = r * r2inv;

                  /* Calculate the energy and derivatives according to dield. */

                  if (dield == -3) {

                     /* special code Ramstein & Lavery dielectric, 94 force field */

                     rs = SIG * r;
                     rssq = rs * rs;
                     pow = exp(-rs);
                     eps1 = rssq + rs + rs + 2.0;
                     epsi = 1.0 / (DIW - C1 * pow * eps1);
                     cgijr = cgi * cgj * rinv * epsi;
                     elec += cgijr;
                     df2 = -cgijr * (1.0 + C1 * pow * rs * rssq * epsi);
                     df = df2 * rinv;
                  

                  } else if (dield == -4) {

                     /* distance-dependent dielectric code, 94 ff */
                     /* epsilon = r  */

                     rs = (cgi * cgj * r2inv);
                     df2 = -2.0 * rs;
                     elec += rs;
                     df = df2 * rinv;
                  

                  } else if (dield == -5) {

                     /* non-bonded term from yammp  */

                     df = 0.0;
         
                  } else {

                     /*
                      * Code for various dielectric models.
                      * The df2 variable should hold r(dV/dr).
                      */

                     if (dield == 0) {

                        /* epsilon = r  */
		   
                        rs = cgi * cgj * r2inv;
                        df2 = -2.0 * rs;
                        elec += rs;

                     } else if (dield == 1) {

                        /* epsilon = 1  */

                        rs = (cgi * cgj * rinv);
                        df2 = -rs;
                        elec += rs;

                     } else if (dield == -2) {

                        /* Ramstein & Lavery dielectric, PNAS 85, 7231 (1988). */

                        rs = SIG * r;
                        rssq = rs * rs;
                        pow = exp(-rs);
                        eps1 = rssq + rs + rs + 2.0;
                        epsi = 1.0 / (DIW - C1 * pow * eps1);
                        cgijr = cgi * cgj * rinv * epsi;
                        elec += cgijr;
                        df2 = -cgijr * (1.0 + C1 * pow * rs * rssq * epsi);
                     }

                     /* Ramu: no vdw force calculated for residue/chain appr */
                     df = df2 * rinv;
 	  
                  }

                  /*
                   * The de term contains one more factor of Dij in the denominator
                   * so that terms such as dedx do not need to include 1/Dij. 
                   *
                   * Update the gradient for atom j.
                   */

                  df *= rinv;

                  dedx = df * xij;
                  dedy = df * yij;
                  dedz = df * zij;

                  dumx += dedx;
                  dumy += dedy;
                  dumz += dedz;

               }  /* end-if q > QMIN */
 
            }  /* end of for each approximate charge */
	    
         }  /* end-for each level 1 pairlist */

         /* For atom i, the gradient is updated in the i-loop only. */

         f[foff + dim * i + 0] += dumx;
         f[foff + dim * i + 1] += dumy;
         f[foff + dim * i + 2] += dumz;
         if (dim == 4) {
            f[foff + dim * i + 3] += dumw;
         }

//         printf("hcp f on %d = %f, %f, %f\n", i, f[foff+dim*i+0], f[foff+dim*i+1], f[foff+dim*i+2]);
         
      }  /* end-for each atom */

      /* Deallocate the iexw array within this potentially parallel region. */
      free_ivector(iexw, -1, prm->Natom+1);
   }

   /* Return evdw and elec through by-reference calling parameters. */

   *enb = evdw/2.0; /* comment out if VDW calculated separately */
   *eel = elec/2.0;
   
//   printf("nbond hcp done\n");
//   fflush(stdout);

   return (0);
}

