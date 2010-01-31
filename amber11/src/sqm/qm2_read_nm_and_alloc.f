! <compile=optimized>
#include "copyright.h"
#include "../include/assert.fh"
#include "../include/dprec.fh"

!+++++++++++++++++++++++++++++++++++++++++++++
!This subroutine reads the QMMM namelist
!and also calls allocation routines
!for QMMM based on natom.
!
!Author:
!     Ross Walker
!+++++++++++++++++++++++++++++++++++++++++++++

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Reads the qmmm namelist and calls the qmmm memory allocation routines
#ifdef SQM
subroutine read_qmmm_nm_and_alloc( natom_in, igb, atnam, atnum, maxcyc, &
            grms_tol, ntpr )
#else
subroutine read_qmmm_nm_and_alloc( igb, ih, ix, x, cut, use_pme, ntb )
#endif

   use findmask
   use constants, only : RETIRED_INPUT_OPTION
   use qmmm_module, only : qmmm_struct, qm2_struct, qmmm_nml, &
                           AM1, PM3, MNDO, RM1, PDDGPM3_08, PDDGPM3, &
                           PDDGMNDO, DFTB, PM3CARB1, PM6, nelements, &
                           set_qmtheory, validate_qm_atoms, qmsort, &
                           allocate_qmmm, get_atomic_number, qmmm_div
   implicit none

!STATIC MEMORY
integer :: max_quantum_atoms  !Needed in read qmmm namelist since namelists cannot contain pointers
parameter ( max_quantum_atoms = 10000 )
!END STATIC MEMORY

!Passed in
   integer :: igb        !Value of igb from cntrl namelist
#ifdef SQM
   integer use_pme, ntb
   integer, intent(in) :: natom_in 
   character(len=8), intent(in) :: atnam(*)
   integer, intent(in)  :: atnum(*)
   integer, intent(out) :: maxcyc, ntpr
   _REAL_, intent(out)  :: grms_tol
#else
   character(len=4) ih(*)
   integer, intent(in) :: ix(*)
   _REAL_, intent(in) :: x(*)
   _REAL_, intent(in) :: cut !MM-MM cutoff in angstroms
   integer, intent(in) :: use_pme, ntb
#endif

!local  
   _REAL_ :: qmcut      ! local copied to qmmm_nml%qmcut - specified cutoff to use for QM-MM electrostatics.
                         ! Default = same as regular MM cutoff.
   _REAL_ :: lnk_dis     ! Distance from the QM atom to place link atom.
                         !A value of <0.0 means the link atom gets placed on the MM link pair atom's coordinates
                         !on every MD step.
   _REAL_ :: scfconv     ! local copied to qmmm_nml%scfconv - Convergence criteria for SCF routine. Default = 1.0D-8.
                         ! Minimum (tightest criteria) = 1.0D-16
   _REAL_ :: pseudo_diag_criteria !Criteria - maximum change in density matrix between successive SCF iterations
                                  !in which to allow pseudo diagonalisation. Default = 0.05.
   integer :: lnk_atomic_no !Atomic number of link atom
   integer :: lnk_method !controls how QM-MM valence terms are dealt with.
   integer :: qmgb       ! local copied to qmmm_nml%qmgb - flag for type of GB do with QM region
   integer :: qmtheory   ! local copied to qmmm_nml%qmtheory - flag for level of theory to use for QM region
   integer :: qmcharge   ! local copied to qmmm_nml%qmcharge - value of charge on QM system
   integer :: spin       ! local copied to qmmm_nml%spin - spin state of system
   integer :: i,j        ! temporary counter for local loops
   integer :: ifind
   integer :: qmqmdx     ! local copied to qmmm_nml%qmqm_analyt - 1 = analytical, 2 = numerical QM-QM derivatives in qm2
   integer :: verbosity  ! local copied to qmmm_nml%verbosity - Controls amount of info about QM part of calc that is printed (0=Def)
   integer :: tight_p_conv ! local copied to qmmm_nml%tight_p_conv - Controls convergence of density matrix. 0 (Def) = 0.05*sqrt(SCFCRT)
                           ! 1 = Converged both Energy and Density to SCFCONV
   integer :: printcharges !Local copied to qmmm_nml%printcharges as a logical. 1 = true - print mulliken and cm1a and cm2a charges
                           !on every step. 0 = false = don't print charges. Default = 0 (.false.)
   integer :: peptide_corr !Local copied to the logical qmmm_nml%peptide_corr
                           !Add MM correction to peptide linkages 0 = No (Default), 1 = Yes.
   integer :: itrmax       !Local copied to qmmm_nml%itrmax - Maximum number of scf cycles to run
                           !before assuming convergence has failed (default = 1000)
   integer :: qmshake      !Local copied to qmmm_nml%qmshake - shake QM atoms if ntc>1?
   integer :: qmmmrij_incore !Flag to store rij between qm-mm pairs and related equations in memory.
                             !1 (default) = store in memory. 0 = calc on fly.
   integer :: qmqm_erep_incore !Flag to store QM-QM 1 electron repulsion integrals in memory or to calculate
                               !them on the fly. Only available with QM-QM analytical derivatives.
                               !1 (default) = store in memory. 0 = calc on fly.
   integer :: pseudo_diag      !Whether to allow pseudo diagonalisations to be done when possible in SCF.
                               !0 (default) = Always do full diagonalisations.
                               !1 = do pseudo diagonalisations when possible.
   integer :: qm_ewald          !0 (default) do only regular QM-MM interaction in periodic calculations.
                                !1           do ewald based periodic QM-MM interactions.
                                !2           do ewald based periodic QM-MM but with the QM image charges
                                !            fixed at the previous steps mulliken charges during the SCF.
   integer :: qm_pme            !0 use regular Ewald for doing QM-MM interactions when qm_ewald>0.
                                !1 (default) use PME to do the reciprocal sum.
   integer :: kmaxqx, kmaxqy, kmaxqz !Maximum K space vectors
   integer :: ksqmaxq !Maximum K squared values for spherical cutoff in k space.
   integer :: writepdb
   integer :: qmmm_int !QM-MM interaction method
   integer :: adjust_q
   integer :: diag_routine !Controls diagonalization routine to use in SCF.
#ifdef OPENMP
   integer :: qmmm_omp_max_threads !Maximum number of openmp threads to use for parallel QMMM routines
#endif
   integer :: density_predict !Controls prediction of density matrix for next SCF step.
   integer :: fock_predict !Controls prediction of Fock matrix for next SCF step.
   integer :: nearest_qm_solvent !0 by default, set >0 to put this many solvent molecules in the QM region.
   integer :: nearest_qm_solvent_fq !frequency to recheck nearest qm solvents.
   character(len=4) :: nearest_qm_solvent_resname !Resname of the solvent residues - e.g. WAT

   _REAL_ :: fockp_d1 !prefactor for fock matrix prediction.
   _REAL_ :: fockp_d2 !prefactor for fock matrix prediction.
   _REAL_ :: fockp_d3 !prefactor for fock matrix prediction.
   _REAL_ :: fockp_d4 !prefactor for fock matrix prediction.
   logical :: mdin_qmmm=.false.

   integer :: idc
   integer :: divpb
!! (GMS)
   _REAL_  :: chg_lambda       ! Charge scaling factor for free energy calculation
!! DFTB options
   integer :: dftb_maxiter     ! Max # of iterations before resetting Broyden (default: 70 ) ==> qmmm_nml%dftb_maxiter
   integer :: dftb_disper      ! Use dispersion?  (default: 0 = false) ==> qmmm_nml%dftb_disper
   integer :: dftb_chg         ! DFTB CM3 charges (default: 0 = Mulliken, 1 = CM3) ==> qmmm_nml%dftb_chg
   _REAL_  :: dftb_telec       ! Electronic temperature, in Kelvins. (Default = 0.0K) ==> qmmm_nml%dftb_telec 
   _REAL_  :: dftb_telec_step  ! Telec step size for convergence accelerator (Default = 0.0K) ==> qmmm_nml%dftb_telec_step
   character(Len=256) :: dftb_3rd_order  ! 3rd order SCC-DFTB (default: 'NONE'== No third order)
                                       !     'PA' == Do 3rd order, Proton Affinities parameterization
                                       !     'PR' ==               Phosphate reactions parameterization
                                       !     'READ' == read the parameters from a user-specified file (TO IMPLEMENT)

#include "memory.h"

   !Apparently you can't use a pointer in a namelist :-( Therefore
   !we need a local scratch array that will be big enough that 
   !the iqmatoms list never exceeds it
   integer :: iqmatoms( max_quantum_atoms )

   character(len=1024) :: qmmask
   character(len=12) :: qm_theory 
        !Options=PM3,AM1,MNDO,PDDG-PM3,PM3PDDG,PDDG-MNDO,PDDGMNDO,
        !        PM3-CARB1,PM3CARB1,DFTB,SCC-DFTB,RM1, PM6
   integer, dimension(:), pointer :: isqm
   integer :: ier=0

   namelist /qmmm/ qmcut, iqmatoms,qmmask,qmgb,qm_theory, qmtheory, &
                   qmcharge, qmqmdx, verbosity, tight_p_conv, scfconv, &
                   printcharges, peptide_corr, itrmax, qmshake, &
                   qmqm_erep_incore, qmmmrij_incore, &
                   lnk_dis, lnk_atomic_no, lnk_method, spin, pseudo_diag,   &
                   pseudo_diag_criteria, &
                   qm_ewald, qm_pme, kmaxqx, kmaxqy, kmaxqz, ksqmaxq, &
                   writepdb, qmmm_int, adjust_q, diag_routine, &
                   density_predict, fock_predict, &
                   fockp_d1, fockp_d2, fockp_d3, fockp_d4, idc, divpb, &
                   dftb_maxiter, dftb_disper, dftb_3rd_order, dftb_chg, &
                   dftb_telec, dftb_telec_step, &
#ifdef OPENMP
                   qmmm_omp_max_threads, &
#endif
#ifdef SQM
                   maxcyc, ntpr, grms_tol,  &
#endif
                   chg_lambda, nearest_qm_solvent, nearest_qm_solvent_fq, &
                   nearest_qm_solvent_resname

!Setup defaults
#ifdef SQM
   qmcut = 9999.d0
   use_pme = 0
   ntb = 0
   maxcyc = 9999
   grms_tol = 0.02
   ntpr=10
#else
   qmcut = cut
#endif
   lnk_dis=1.09d0  !Methyl C-H distance
   lnk_atomic_no=1 !Hydrogen
   lnk_method=1 !treat MMLink as being MM atom.
   qmgb = 2 !Gets set to zero if igb==6 or igb==0.
   qm_theory = ''
   qmtheory = RETIRED_INPUT_OPTION 
   qmcharge = 0
   spin = 1
   qmqmdx = 1
   verbosity = 0
   tight_p_conv = 0
   scfconv = 1.0D-8
   printcharges = 0
   peptide_corr = 0
   itrmax = 1000
   qmshake = 1
   qmmask=''
   iqmatoms(1:max_quantum_atoms) = 0
   qmmmrij_incore = 1
   qmqm_erep_incore = 1
   pseudo_diag = 1
   pseudo_diag_criteria = 0.05d0
   qm_ewald=1 !Default is to do QMEwald, with varying charges, if ntb=0 or use_pme=0 then this will get turned off
   qm_pme = 1 !use pme for QM-MM
   kmaxqx=5; kmaxqy=5; kmaxqz=5    !Maximum K space vectors
   ksqmaxq=27 !Maximum K squared values for spherical cutoff in k space.
   writepdb = 0 !Set to 1 to write a pdb on the first step with just the QM region in it.
   qmmm_int = 1 !Default, do full interaction without extra Gaussian terms for PM3 / AM1 etc.
   adjust_q = 2 !Default adjust q over all atoms.
   diag_routine = 1 !Use default internal diagonalizer.
#ifdef OPENMP
   qmmm_omp_max_threads = 1 !Use just 1 openmp thread by default.
#endif
   density_predict = 0 !Use density matrix from previous MD step.
   fock_predict = 0 !Do not attempt to predict the Fock matrix.
   fockp_d1 = 2.4d0
   fockp_d2 = -1.2d0
   fockp_d3 = -0.8d0
   fockp_d4 = 0.6d0
   idc = 0
   divpb = 0
   nearest_qm_solvent = 0 !Do not change the QM region to keep the nearest solvents.
   nearest_qm_solvent_fq = 1 !Do update of QM solvent on every step if nearest_qm_solvent > 0
   nearest_qm_solvent_resname='WAT ' !by default assume solvent is resname WAT

   !DFTB
   dftb_maxiter     = 70   
   dftb_disper      = 0
   dftb_chg         = 0
   dftb_telec       = 0.0d0
   dftb_telec_step  = 0.0d0
   chg_lambda  = 1.0d0
   dftb_3rd_order   = 'NONE'


   !Read qmmm namelist
   rewind 5

   call nmlsrc('qmmm',5,ifind)
   if (ifind /= 0) mdin_qmmm=.true.

   !Read qmmm namelist
   rewind 5
   if ( mdin_qmmm ) then
     read(5,nml=qmmm)
   else
     write(6, '(1x,a,/)') 'Could not find qmmm namelist'
     call mexit(6,1)
   endif

   !Parse qm_theory
   call set_qmtheory(qmtheory,qm_theory)
   
#ifdef SQM
      natom = natom_in
      qmmm_struct%nquant = natom
      qmmm_struct%nquant_nlink = natom
      do i=1,natom
         iqmatoms(i) = i
      end do
#else

   if( qmmask /= '' ) then  !  get the quantum atoms from the mask
      write(6,'(a)') ''
      write(6,'(a)') 'LOADING THE QUANTUM ATOMS AS GROUPS'

      allocate(isqm( natom ), stat=ier)
      REQUIRE(ier==0)

      call atommask( natom, nres, 0, ih(m04), ih(m06), &
         ix(i02), ih(m02), x(lcrd), qmmask, isqm )
      qmmm_struct%nquant = sum(isqm(1:natom))
      write(6,'(a,a,a,i5,a)') '     Mask ', qmmask(1:len_trim(qmmask)), &
         ' matches ',qmmm_struct%nquant,' atoms'

      j = 0
      do i=1,natom
         if( isqm(i)>0 ) then
            j = j+1
            iqmatoms(j) = i
         end if
      end do

      deallocate(isqm, stat=ier)
      REQUIRE(ier==0)

   else  !  get the count from the input iqmatoms array

      do i = 1,max_quantum_atoms
         if( iqmatoms(i) == 0 ) exit
      end do
      qmmm_struct%nquant = i-1

   end if

!Initialize nlink to 0
   qmmm_struct%nlink = 0
   qmmm_struct%nquant_nlink = qmmm_struct%nquant

#endif
! Test to see if QM atom selection is legal.
   call validate_qm_atoms(iqmatoms,qmmm_struct%nquant,natom)

!  check we don't bust our statically allocated max_quantum_atoms
   call int_legal_range('QMMM: (number of quantum atoms) ', &
      qmmm_struct%nquant, 1, max_quantum_atoms )

   call qmsort(iqmatoms) !ensure the list of qm atoms is sorted numerically

! --- Variable QM solvent region - has to be very early here because we
!     will be changing nquant and iqmatoms.  ---
   call int_legal_range('QMMM: (QM-MM nearest_qm_solvent) ',nearest_qm_solvent,0,nres)
   if (nearest_qm_solvent > 0 ) then
     call int_legal_range('QMMM: (QM-MM nearest_qm_solvent_fq) ',nearest_qm_solvent_fq,1,99999999)
   end if
   qmmm_nml%nearest_qm_solvent = nearest_qm_solvent
   qmmm_nml%nearest_qm_solvent_fq = nearest_qm_solvent_fq
   qmmm_nml%nearest_qm_solvent_resname = nearest_qm_solvent_resname
   if (qmmm_nml%nearest_qm_solvent > 0) then
#ifdef SQM
     write(6,*) 'SQM does not support the nearest solvent stuff (yet)'
     call mexit(6,1)
#else
     !We need to work out how many atoms nearest_qm_solvent * natoms_per_solvent_residue equals.
     !We then need to find the nearest atoms and update nquant and iqmatoms respectively.
     call qm2_setup_vsolv(qmmm_struct%nquant, max_quantum_atoms, iqmatoms, &
                          nres, ih(m02), ix(i02),ix(i70), natom)
     ! check again that we don't bust our statically allocated max_quantum_atoms
     call int_legal_range('QMMM: (number of quantum atoms) ', &
      qmmm_struct%nquant, 1, max_quantum_atoms )
#endif
   end if
! --- End Variable QM water region ---

   call float_legal_range('QMMM: (QM-MM Cutoff) ', qmcut,0.0D0,1.0D30)
   call int_legal_range('QMMM: (QM GB Method) ', qmgb,0,3)
   call int_legal_range('QMMM: (QM Theory) ', qmtheory,1,10)
   call int_legal_range('QMMM: (QM-QM Derivatives) ', qmqmdx,1,2)
   call int_legal_range('QMMM: (Verbosity) ', verbosity,0,5)
   call int_legal_range('QMMM: (Max SCF Iterations) ', itrmax,1,10000000)
   call int_legal_range('QMMM: (Shake on QM atoms) ', qmshake,0,1)
   call int_legal_range('QMMM: (Density Matrix Convergence) ', tight_p_conv,0,1)
   call float_legal_range('QMMM: (SCF Convergence) ', scfconv,1.0D-16,1.0D0)
   call int_legal_range('QMMM: (PRINT CHARGES) ', printcharges,0,1)
   call int_legal_range('QMMM: (Spin State) ', spin,1,1)
!RCW: Currently limit spin state to singlets only since the code for spin>1 does not exist / work at present.
!     WARNING - IF WE LATER ALLOW SPIN>1 qm2_densit will need updating.
   call int_legal_range('QMMM: (Peptide Correction) ',peptide_corr,0,1)
   call int_legal_range('QMMM: (QM-MM RIJ in Core) ',qmmmrij_incore,0,1)
   call int_legal_range('QMMM: (QM-QM E-Rep in Core) ',qmqm_erep_incore,0,1)
   call int_legal_range('QMMM: (Link Atomic Number) ',lnk_atomic_no,1,nelements)
   call int_legal_range('QMMM: (QM-MM Link Method) ',lnk_method,1,2)
   call int_legal_range('QMMM: (Pseudo Diag) ',pseudo_diag,0,1)
   call int_legal_range('QMMM: (QM Ewald) ',qm_ewald,0,2)
   call int_legal_range('QMMM: (QM PME) ',qm_pme,0,1)
   call int_legal_range('QMMM: (QM Ewald kmaxqx) ',kmaxqx,1,99999999)
   call int_legal_range('QMMM: (QM Ewald kmaxqy) ',kmaxqy,1,99999999)
   call int_legal_range('QMMM: (QM Ewald kmaxqz) ',kmaxqz,1,99999999)
   call int_legal_range('QMMM: (QM Ewald ksqmaxq) ',ksqmaxq,1,kmaxqx*kmaxqy*kmaxqz)
   call int_legal_range('QMMM: (QM-MM qmmm_int) ',qmmm_int,0,2)
   call int_legal_range('QMMM: (QM-MM adjust_q) ',adjust_q,0,2)
   call int_legal_range('QMMM: (QM-MM diag_routine) ',diag_routine,0,7)
#ifdef OPENMP
   call int_legal_range('QMMM: (QM-MM qmmm_omp_max_threads) ',qmmm_omp_max_threads,1,32)
#endif
   call int_legal_range('QMMM: (QM-MM density_predict) ',density_predict,0,1)
   call int_legal_range('QMMM: (QM-MM fock_predict) ',fock_predict,0,1)
   call float_legal_range('QMMM: (Pseudo Diag Criteria) ',pseudo_diag_criteria,1.0D-12,1.0D0)
   if (lnk_dis>0.0d0) then
     !if lnk_dis is less than 0.0d0 then the link atom is just placed on top of
     !the MM link pair atom.
     call float_legal_range('QMMM: (Link Atom Distance) ',lnk_dis,0.7D0,4.0D0)
   endif

!! GMS
   call float_legal_range('QMMM: (QM-MM chg_lambda)'     , chg_lambda     , 0.0D0 , 1.0D0  )
   qmmm_nml%chg_lambda  = chg_lambda
!! DFTB
   call int_legal_range(  'QMMM: (QM-MM dftb_maxiter ) ' , dftb_maxiter   , 1     , 10000  )
   call int_legal_range(  'QMMM: (QM-MM dftb_disper) '   , dftb_disper    , 0     , 1      )
   call int_legal_range(  'QMMM: (QM-MM dftb_chg   ) '   , dftb_chg       , 0     , 1      )
   call float_legal_range('QMMM: (QM-MM dftb_telec)'     , dftb_telec     , 0.0D0 , 1.0D4  )

   if (dftb_3rd_order /= 'NONE') then
      call check_dftb_3rd_order(dftb_3rd_order)
      qmmm_nml%dftb_3rd_order = dftb_3rd_order
   endif

   qmmm_nml%dftb_maxiter   = dftb_maxiter
   qmmm_nml%dftb_disper      = dftb_disper
   qmmm_nml%dftb_chg         = dftb_chg
   qmmm_nml%dftb_telec       = dftb_telec
   qmmm_nml%dftb_telec_step  = dftb_telec_step

   if (dftb_chg > 0) printcharges=1

   qmmm_nml%qmcut = qmcut
   qmmm_nml%qmcut2 = qmcut*qmcut
   qmmm_nml%lnk_dis = lnk_dis
   qmmm_nml%lnk_atomic_no = lnk_atomic_no
   qmmm_nml%lnk_method = lnk_method

!DFTB Limitations - current things not supported in DFTB
!These are silent limitations that are non fatal - just to
!avoid problems with the default. Fatal errors are handled
!later in this routine.
   if (qmtheory == DFTB) then
     qmmmrij_incore=0
     qmqm_erep_incore=0
     pseudo_diag=0
     qmqmdx=1
   end if

!Divcon limitations - current things not supported in sander.DIVCON
!These are silent changes to remove defaults. Fatal errors are handled
!later in this routine.
   if (idc /= 0) then
     qmmmrij_incore=0
     qmqm_erep_incore=0
     pseudo_diag=0
   end if

   !qmgb values:
   ! 0 - do GB but leave QM charges as zero and add nothing to Fock matrix. This is like a vacuum
   !     QM molecule in a solvated MM system.
   ! 1 - do GB using the prmtop fixed resp charges for the GB calculation.
   ! 2 - do GB using Mulliken charges that are consistent with the GB field by modifying the fock
   !     matrix at every SCF step. (default)
   ! 3 - do GB using QM gas phase Mulliken charges - This is really a debugging option since the charges
   !     will not be consistent with the GB field since the fock matrix is not modified. This similarly
   !     means that the gradients will not be accurate. A warning will be printed at every QM call if this
   !     option is selected.

   !Make sure igb in &cntrl namelist is compatible with qmgb setting.
   if (igb==0 .or. igb==6) then
      !no qmgb available
      qmgb = 0
   end if
   !Print warning about qmgb being for debugging only.
   if (qmgb==3) then
     write(6,*) "QMMM: ------------------------------ WARNING --------------------------------"
     write(6,*) "QMMM: qmgb = 3 is designed for debugging purposes only. It gives GB"
     write(6,*) "QMMM:          energies based on gas phase QM Mulliken charges. These charges"
     write(6,*) "QMMM:          are NOT consistent with the GB field felt by the QM region and"
     write(6,*) "QMMM:          so any gradients calculated using this approach will"
     write(6,*) "QMMM:          NOT BE ACCURATE."
     write(6,*) "QMMM:          This option is really designed for:"
     write(6,*) "QMMM:                SINGLE POINT ENERGY EVALUATIONS ONLY"
     write(6,*) "QMMM: ------------------------------ WARNING --------------------------------"
   end if
   qmmm_nml%qmgb = qmgb

   qmmm_nml%qmtheory = qmtheory
   qmmm_struct%AM1_OR_PM3 = (qmmm_nml%qmtheory == AM1 .OR. qmmm_nml%qmtheory == PM3 &
                             .OR. qmmm_nml%qmtheory == PDDGPM3 .OR. qmmm_nml%qmtheory == PM3CARB1 .OR. &
                             qmmm_nml%qmtheory == RM1 .OR. qmmm_nml%qmtheory == PDDGPM3_08 .OR. qmmm_nml%qmtheory == PM6)
   qmmm_struct%PDDG_IN_USE = (qmmm_nml%qmtheory == PDDGPM3 .OR. qmmm_nml%qmtheory == PDDGMNDO  &
                             .OR. qmmm_nml%qmtheory == PDDGPM3_08 )
   qmmm_nml%qmcharge = qmcharge
   qmmm_nml%spin = spin
   qmmm_nml%verbosity = verbosity
   qmmm_nml%itrmax = itrmax
   qmmm_nml%qmshake = qmshake
   qmmm_nml%pseudo_diag_criteria = pseudo_diag_criteria
   if (qmqmdx /= 1) then
      qmmm_nml%qmqm_analyt = .false. !Do numerical QM-QM derivatives in qm2
   else
      qmmm_nml%qmqm_analyt = .true.  !Do analytical QM-QM dericatives in qm2
   end if
   if (tight_p_conv /= 1) then
      qmmm_nml%tight_p_conv = .false. !Loose density matrix convergence (0.05*sqrt(SCFCRT))
   else
      qmmm_nml%tight_p_conv = .true.  !Tight density matrix convergence (SCFCRT)
   end if
   !Write a warning about excessively tight convergence requests.
   if ( scfconv < 1.0D-12 ) then
     write(6,'(" QMMM: WARNING - SCF Conv = ",G8.2)') scfconv
     write(6,*) "QMMM:           There is a risk of convergence problems when the"
     write(6,*) "QMMM:           requested convergence is less that 1.0D-12 kcal/mol."
   end if
   qmmm_nml%scfconv = scfconv
   !How tight do we want the density convergence?
   if (qmmm_nml%tight_p_conv) then
      qmmm_nml%density_conv = qmmm_nml%scfconv
   else
      qmmm_nml%density_conv = 0.05D0 * sqrt(qmmm_nml%scfconv)
   end if

   if ( printcharges /= 1) then
      qmmm_nml%printcharges=.false.
   else
      qmmm_nml%printcharges=.true.
   end if
   if ( peptide_corr == 0) then
      qmmm_nml%peptide_corr = .false.
   else
      qmmm_nml%peptide_corr =  .true.
   end if 
   if ( qmmmrij_incore == 0 .or. qmmm_int==0 ) then
      qmmm_nml%qmmmrij_incore = .false.
   else
      qmmm_nml%qmmmrij_incore = .true. !Only available with qmmm_int=1 or qmmm_int=2
   end if
   if ( qmqm_erep_incore == 0 .or. qmqmdx == 2 ) then
      qmmm_nml%qmqm_erep_incore = .false.
   else
      !Only available with analytical derivatives.
      qmmm_nml%qmqm_erep_incore = .true.
   end if
   if ( pseudo_diag == 1 ) then
     qmmm_nml%allow_pseudo_diag = .true.
   else
     qmmm_nml%allow_pseudo_diag = .false.
   end if
 
   qmmm_nml%qm_ewald = qm_ewald
   qmmm_nml%ksqmaxq = ksqmaxq
   qmmm_nml%kmaxqx = kmaxqx
   qmmm_nml%kmaxqy = kmaxqy
   qmmm_nml%kmaxqz = kmaxqz
   !If ntb=0 or use_pme =0 then we can't do qm_ewald so overide what the user may
   !have put in the namelist and set the value to false.
   if (ntb==0 .or. use_pme==0) then
     qmmm_nml%qm_ewald = 0
     qmmm_nml%qm_pme = .false.
   end if
   if (qmmm_nml%qm_ewald>0 .and. qm_pme>0) then
     qmmm_nml%qm_pme=.true.
   else
     qmmm_nml%qm_pme=.false.
   end if

   if ( writepdb == 0 ) then
     qmmm_nml%writepdb=.false.
   else
     qmmm_nml%writepdb=.true.
   end if
   qmmm_nml%qmmm_int = qmmm_int

   qmmm_nml%idc = idc
   qmmm_nml%divpb = divpb

!Setup some specific calculation flags that depend on namelist variables.
!Need to make sure these get copied to other threads in an MPI run.

   !Will we be calculating the Mulliken charges on every SCF iteration?
   !Default is no. Will be set to true in a bit if certain options, such as qm_ewald
   !require it.
#ifdef SQM
   qm2_struct%calc_mchg_scf = .true.
#else
   if (qmmm_nml%qm_ewald==1 .or. qmmm_nml%qmgb>1) then
     !We will be needing the mulliken charges on every SCF iteration
     qm2_struct%calc_mchg_scf = .true.
   else
     qm2_struct%calc_mchg_scf = .false.
   end if
#endif

   !DFTB Calculates Mulliken charges anyway so we might as well store them in the correct place.
   if (qmmm_nml%qmtheory == DFTB) qm2_struct%calc_mchg_scf = .true.

!At this point we know nquant and natom so we can allocate our arrays that depend on nquant or natom
!Note if this is a LES run qmmm_struct%nquant is nqaunt
!Note non master mpi threads need to call this allocation routine manually themselves.
   call allocate_qmmm( natom )

#ifdef SQM
   qmmm_struct%iqm_atomic_numbers(1:natom) = atnum(1:natom)
   qmmm_nml%iqmatoms(1:natom) = iqmatoms(1:natom)
#else
#   ifndef NO_SANDER_DIVCON
   !DIVCON SPECIFIC STUFF
   if(qmmm_nml%divpb == 1)then
      !the +100 is for link atoms, needs some way to figure out actual number of link atoms
      allocate(qmmm_div%all_atom_numbers(natom+100), stat=ier)
      REQUIRE(ier == 0)

      do i=1,natom
         call get_atomic_number(ih(m04+i-1),x(lmass+i-1),qmmm_div%all_atom_numbers(i))
      enddo
   endif
   !END DIVCON SPECIFIC STUFF
#   endif

   do i = 1, qmmm_struct%nquant
      qmmm_nml%iqmatoms(i) = iqmatoms(i)
      !Get the atomic numbers (used to be done in rdparm2...)
      j = iqmatoms(i)
      call get_atomic_number( ih(m04+j-1),x(lmass+j-1),qmmm_struct%iqm_atomic_numbers(i) )
   end do
#endif

!Now we have a list of atom numbers for QM atoms we can build a true false (natom long) list
!specifying what the quantum atoms are. Useful for doing quick .OR. operations against other
!lists.
   qmmm_struct%atom_mask = .false. !Note, sets entire natom long array to false

   do i = 1, qmmm_struct%nquant
     qmmm_struct%atom_mask(qmmm_nml%iqmatoms(i)) = .true.
   end do

   qmmm_nml%adjust_q = adjust_q

   qmmm_nml%diag_routine = diag_routine
#ifdef OPENMP
   qmmm_nml%qmmm_omp_max_threads = qmmm_omp_max_threads

   !For the time being the number of threads to use for diag and pdiag
   !routines is set to max_threads - later this will be optimized if
   !diag_routine=0.
   qmmm_omp%diag_threads = qmmm_omp_max_threads
   qmmm_omp%pdiag_threads = qmmm_omp_max_threads
#endif
   qmmm_nml%density_predict = density_predict
   qmmm_nml%fock_predict = fock_predict
   qmmm_nml%fockp_d1 = fockp_d1
   qmmm_nml%fockp_d2 = fockp_d2
   qmmm_nml%fockp_d3 = fockp_d3
   qmmm_nml%fockp_d4 = fockp_d4

! --- CHECK FOR LIMITATIONS ---

  !--- You cannot mix Fock prediction with density prediction. ---
  if (qmmm_nml%fock_predict > 0 .and. qmmm_nml%density_predict > 0) then
      call sander_bomb('read_qmmm_nm_and_alloc','Fock matrix and Density matrix prediction are mutually exclusive.', &
                       'Cannot have fock_predict > 0 and density_predict > 0')
  end if

  !--- For Fock prediction the 4 pre-factors must sum to 1.0d0 ---
  if (qmmm_nml%fock_predict == 1) then
    if (abs(1.0d0-(qmmm_nml%fockp_d1 + qmmm_nml%fockp_d2 + qmmm_nml%fockp_d3 + qmmm_nml%fockp_d4)) > 1.0d-6) then
       write(6,*) 'QMMM: Failure, fockp_d1 to d4 must sum to 1.0d0 - current sum is', &
                  (qmmm_nml%fockp_d1 + qmmm_nml%fockp_d2 + qmmm_nml%fockp_d3 + qmmm_nml%fockp_d4)
       call sander_bomb('read_qmmm_nm_and_alloc','Fock matrix prediction coefficients do not sum to 1.0.', &
                         'adjust fockp_d1 to fockp_d4 so that they sum to 1.0.')
    end if
  end if

  !--- You cannot use variable solvent with GB calculations.
  if (qmmm_nml%qmgb > 0 .and. qmmm_nml%nearest_qm_solvent > 0) then
      call sander_bomb('read_qmmm_nm_and_alloc','Nearest QM solvent and qmgb are mutually exclusive.', &
                       'Cannot have nearest_qm_solvent > 0 and qmgb > 0')
  end if

  !--- DIVCON LIMITATIONS ---
  !Divcon currently only works with gas phase simulations.
  if (qmmm_nml%idc>0) then
    if (ntb /= 0) then
      !This covers qmewald, qm_pme as well.
      call sander_bomb('read_qmmm_nm_and_alloc','idc /= 0 (Divcon on) but periodic boundaries are in use.', &
                       'Periodic boundaries are currently only supported when idc == 0')
    end if
    if (qmmm_nml%qmgb/=0) then
      call sander_bomb('read_qmmm_nm_and_alloc','idc /= 0 (Divcon on) but qmgb/=0. QMMM GB is currently', &
                       'only supported when idc == 0')
    end if
    if (qmmm_nml%peptide_corr) then
      call sander_bomb('read_qmmm_nm_and_alloc','idc /= 0 (Divcon on) but peptide_corr /= 0.', &
                       'Peptide correction for divcon is handled by the divcon.in file. Set peptide_corr = 0 to proceed.')
    end if
    if (qmmm_nml%printcharges) then
      call sander_bomb('read_qmmm_nm_and_alloc','idc /= 0 (Divcon on) but printcharges /= 0.', &
                       'Printcharges is not available with idc > 0.')
    end if
    if (qmmm_nml%nearest_qm_solvent > 0) then
        call sander_bomb('read_qmmm_nm_and_alloc','Nearest QM solvent is not available with idc>0.', &
                         'Cannot have nearest_qm_solvent > 0 and idc > 0')
    end if
#ifdef MPI
    !No support for parallel Divcon
    write(6,*) 'Divcon capability (idc>0) can only run in serial mode for now'
    call mexit(6,1)
#endif
    !Write a warning about qmtheory being ignored with Divcon.
    write(6,'("|QMMM: WARNING DIVCON IN USE")')
    write(6,'("|QMMM: qm_theory IS IGNORED WHEN USING DIVCON - QM HAMILTONIAN MUST BE SELECTED")')
    write(6,'("|QMMM: IN DIVCON.IN FILE.")')
  end if
  !--- END DIVCON LIMITATIONS ---

  !--- DFTB LIMITATIONS ---
  if (qmmm_nml%qmtheory == DFTB ) then
    if (qmmm_nml%peptide_corr) then
      call sander_bomb('read_qmmm_nm_and_alloc','qm_theory=DFTB but peptide_corr /= 0.', &
                       'Peptide correction is not available, or required for  DFTB. Set peptide_corr = 0 to proceed.')
    end if
  end if
  !--- END DFTB LIMITATIONS ---

! --- END CHECK FOR LIMITATIONS ---

  return

end subroutine read_qmmm_nm_and_alloc
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

