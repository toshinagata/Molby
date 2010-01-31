#include "copyright.h"
#include "../include/dprec.fh"
#include "../include/assert.fh"
!
!        SQM  stand-alone quantum program
!

program sqm

   use qmmm_module, only : qmmm_nml,qmmm_struct, qmmm_mpi, qm2_struct, qm_gb

   implicit none

   _REAL_ x(3000), f(3000), escf, reff(1000), onereff(1000), work(18000), scf_mchg(1000)
   character(len=8) atnam(1000)
   _REAL_ born_radii(1000), one_born_radii(1000)
   _REAL_ intdiel, extdiel, Arad
   integer natom, ier, atnum(1000), xmin_iter
   
   character(len=80) arg ! temp for each of the command line arguments
   integer iarg !         index of the current argument
   integer last_arg_index !   index of the last argument
   integer ntpr
   character owrite
   character(len=256) mdin, mdout 

   integer igb, maxcyc
   _REAL_ grms_tol
   logical :: master=.true.

   ! ==== Initialise first_call flags for QMMM ====
   qmmm_struct%qm_mm_first_call = .true.
   qmmm_struct%fock_first_call = .true.
   qmmm_struct%fock2_2atm_first_call = .true.
   ! qmmm_struct%qm2_deriv_qm_analyt_first_call = .true.
   qmmm_struct%qm2_allocate_e_repul_first_call = .true.
   ! qmmm_struct%qm2_rotate_qmqm_first_call = .true.
   qmmm_struct%qm2_calc_rij_eqns_first_call = .true.
   qmmm_struct%qm2_scf_first_call = .true.
   qmmm_struct%zero_link_charges_first_call = .true.
   qmmm_struct%adj_mm_link_pair_crd_first_call = .true.

   !     --- default file names ---
   
   mdin   = 'mdin'
   mdout  = 'mdout'
   iarg = 0
   owrite = 'N'  ! output status: New
   last_arg_index = command_argument_count()
   do while (iarg < last_arg_index)

      iarg = iarg + 1
      call getarg(iarg,arg)

      if (arg == '-i') then
         iarg = iarg + 1
         call getarg(iarg,mdin)
      else if (arg == '-o') then
         iarg = iarg + 1
         call getarg(iarg,mdout)
      else if (arg == '-O') then
         owrite = 'R'   ! output status: Replace
      else if (arg == ' ') then
         continue
      else
         write(0,'(/,5x,a,a)') 'Error unknown flag: ',arg
         call mexit(6, 1)
      end if 
   end do  !  while (iarg < last_arg_index)

   igb = 0
   call amopen(5,mdin,'O','F','R')
   call amopen(6,mdout,owrite,'F','W')

   write(6,*) '           --------------------------------------------------------'
   write(6,*) '                             AMBER SQM VERSION 1.0                 '
   write(6,*) ''
   write(6,*) '                                     By'
   write(6,*) '                             Ross C. Walker (SDSC)'
   write(6,*) '                                     and'
   write(6,*) '                            David A. Case (Rutgers)'
   write(6,*) ''              
   write(6,*) '           --------------------------------------------------------'
   write(6,*) ''                  

   call getsqmx( natom, x, atnam, atnum )
   call read_qmmm_nm_and_alloc(natom,igb,atnam,atnum,maxcyc,grms_tol,ntpr )
   call qm_assign_atom_types

   ! Set default QMMM MPI parameters - for single cpu operation.
   ! These will get overwritten by qmmm_mpi_setup if MPI is on.
   ! qmmm_mpi%master = master
   qmmm_mpi%commqmmm_master = master
   qmmm_mpi%numthreads = 1
   qmmm_mpi%mytaskid = 0
   qmmm_mpi%natom_start = 1
   qmmm_mpi%natom_end = natom
   qmmm_mpi%nquant_nlink_start = 1
   qmmm_mpi%nquant_nlink_end = qmmm_struct%nquant_nlink
   call allocate_qmgb(qmmm_struct%nquant_nlink)

   allocate( qmmm_struct%dxyzqm(3, qmmm_struct%nquant_nlink), stat = ier )
   REQUIRE(ier == 0) !Deallocated in deallocate qmmm

   call xmin(natom, x,f,escf, xmin_iter, maxcyc, born_radii, &
       one_born_radii, intdiel, extdiel, Arad, scf_mchg, grms_tol, ntpr )

   write(6,*) 'Final SCF energy is ', escf
   write(6,*) ''
   call qm2_print_charges(1,qmmm_nml%dftb_chg,qmmm_struct%nquant_nlink, &
                            scf_mchg,qmmm_struct%iqm_atomic_numbers)

   write(6,*) ''
   write(6,*) 'Final Structure'
   call qm_print_coords(0,.true.)
   write(6,*)
   
   write(6,*) '          --------- Calculation Completed ----------'
   write(6,*)

   call mexit(6,0)

end program sqm

subroutine sqm_energy(natom,coords,f,escf, &
                 born_radii,one_born_radii, &
                 intdiel, extdiel, Arad, scf_mchg ) 
!
!     Argument list variables:
!
!     coords(natom*3)                 - Cartesian coordinates for all atoms.
!                                       Amber array
!     natom                           - Total number of REAL atoms.
!     qmmm_struct%nquant              - Number of REAL quantum atoms as specified in mdin.
!     iqmatoms(qmmm_struct%nquant)
!                                     - Atom numbers for quantum atoms link atoms given values of -1
!     qmmm_struct%iqm_atomic_numbers(qmmm_struct%nquant) - Atomic numbers for qm atoms.
!     qmmm_struct%nlink               - Number of link atoms.
!     f((natom)*3)                    - Atomic forces.
!     escf                            - Heat of formation from QM.
!     qmmm_struct%qm_coords(3,qmmm_struct%nquant+qmmm_struct%nlink)  
!                                     - Cartesian coordinates of quantum atoms.
!                                       (Extracted from coords by qm_extract_coords)

!     Locally defined arrays:
!     dxyzqm(3,qmmm_struct%nquant+qmmm_struct%nlink)     
!                                - Quantum mechanical derivatives from qm-mm
!                                                        interactions.
!     dxyzcl(3,natom)          - Classical derivatives from qm-mm interaction.
!    born_radii(1->natom)      - Effective GB radii - only used when doing qm with gb (and qm_gb==2)
!                                Calculated via an initial call to egb.
!    one_born_radii(1->natom)  - 1.0d0/born_radii(i)
!    scf_mchg                  - nquant long, gets filled with the mulliken charges during scf.

   use qmmm_module, only : qmmm_nml,qmmm_struct, qm2_struct, qm2_rij_eqns, &
                           element_sym, qm_gb, DFTB, qmmm_mpi, qmmm_scratch
   use constants, only : EV_TO_KCAL, KCAL_TO_EV, zero, one, alpb_alpha
  
   implicit none

#include "../include/assert.fh"

#ifdef MPI
   include 'mpif.h'
#endif

! Passed in
   integer, intent(in) :: natom
   _REAL_ , intent(inout)  :: coords(natom*3) !Amber array - adjusted for link atoms
   _REAL_ , intent(out) :: f(natom*3)
   _REAL_ , intent(out) :: escf
   _REAL_ , intent(in) :: born_radii(natom), one_born_radii(natom)
   _REAL_ , intent(in) :: intdiel, extdiel, Arad
   _REAL_ , intent(inout) :: scf_mchg(qmmm_struct%nquant_nlink)

!Locals
   _REAL_ , dimension(2,3) :: bxbnd
   _REAL_ :: mulliken_charge, total_mulliken_charge, total_energy
   _REAL_ :: alpb_beta
   _REAL_ :: scaled_mm_charges(2)

   integer :: ier=0
   integer i, j, m, offset, qm_no, i3

!Locals for link atoms
   _REAL_ :: forcemod(3)
   integer :: lnk_no, mm_no

!=============================================================================
!                   START OF QMMM SETUP: allocate list memory
!=============================================================================

!  If this is the first call to the routine, do some initial allocation
!  that has not been done elsewhere.
   if (qmmm_struct%qm_mm_first_call) then

     allocate ( qmmm_struct%qm_coords(3,qmmm_struct%nquant_nlink), stat=ier )
                !Stores the REAL and link atom qm coordinates
     REQUIRE(ier == 0)

     !Allocation for QM_GB (qmgb==2)
     if (qmmm_nml%qmgb == 2) then
       !Calculate dielectric factor
       if (qm_gb%alpb_on) then
         alpb_beta=alpb_alpha*(intdiel/extdiel)
         qm_gb%intdieli = one/(intdiel*(one + alpb_beta))
         qm_gb%extdieli = one/(extdiel*(one + alpb_beta))
         qm_gb%one_Arad_beta = alpb_beta/Arad
       else
         qm_gb%intdieli = 1.0d0/intdiel
         qm_gb%extdieli = 1.0d0/extdiel
       end if
       qm_gb%mmcut2 = 999.d0
     end if
   end if ! ---- first call endif ----------

   ! call qm_extract_coords(coords)
   i3 = 0
   do i=1,natom
      qmmm_struct%qm_coords(1,i) = coords(i3+1)
      qmmm_struct%qm_coords(2,i) = coords(i3+2)
      qmmm_struct%qm_coords(3,i) = coords(i3+3)
      i3 = i3 + 3
   end do

!=============================================================================
!                   START OF REST OF QMMM SETUP
!=============================================================================
   if(qmmm_struct%qm_mm_first_call) then 
       if (qmmm_mpi%commqmmm_master) write(6,'(/80(1H-)/''  3.1 QM CALCULATION INFO'',/80(1H-)/)')
       call qm2_load_params_and_allocate() !Load the parameters
             !Also does a lot of memory allocation and pre-calculates all
             !the STO-6G orbital expansions.

       if (qmmm_mpi%commqmmm_master) then
          ! call qm_print_dyn_mem(natom,qmmm_struct%qm_mm_pairs)
          call qm_print_coords(0,.true.)
          !Finally print the result header that was skipped in sander.
           write(6,'(/80(1H-)/''   4.  RESULTS'',/80(1H-)/)')
        end if
   end if !if (qmmm_struct%qm_mm_first_call)

!======================END OF QMMM SETUP ======================================

   !Calculate RIJ and many related equations here. Necessary memory allocation
   !is done inside the routine.
!Parallel
   call qm2_calc_rij_and_eqns(qmmm_struct%qm_coords, qmmm_struct%nquant_nlink, &
          qmmm_struct%qm_xcrd, natom, qmmm_struct%qm_mm_pairs)
                                !and store them in memory to save time later.

   !============================
   ! Calculate SCF Energy
   !============================
   call qm2_energy(escf, scf_mchg, natom, born_radii, one_born_radii, &
                   coords, scaled_mm_charges)

   !=============================
   ! Calculation of Forces
   !=============================

   qmmm_struct%dxyzqm=zero
   if (qmmm_nml%qmtheory==DFTB) then
     call qm2_dftb_get_qm_forces(qmmm_struct%dxyzqm)
   else
     !standard semi-empirical
     call qm2_get_qm_forces(qmmm_struct%dxyzqm)
   end if

   !NOW PUT THE CALCULATED gradient (not force!) INTO THE SANDER FORCE ARRAY
   do i=1,qmmm_struct%nquant
     m = qmmm_nml%iqmatoms(i)
     m = (m-1)*3
     f(m+1) = qmmm_struct%dxyzqm(1,i)
     f(m+2) = qmmm_struct%dxyzqm(2,i)
     f(m+3) = qmmm_struct%dxyzqm(3,i)
   enddo

   !=============================
   !   Print Mulliken Charges
   !=============================

   if (qmmm_nml%printcharges .and. qmmm_mpi%commqmmm_master) then
     call qm2_print_charges(1,qmmm_nml%dftb_chg,qmmm_struct%nquant_nlink, &
                            scf_mchg,qmmm_struct%iqm_atomic_numbers)
   end if

   !=============================
   ! End Print Mulliken Charges
   !=============================

   !Print some extra information if verbosity level is > 0

   if (qmmm_mpi%commqmmm_master .AND. qmmm_nml%verbosity > 0) then
      !Verbosity level of 1 or more = print more accurate SCF energy
      write (6,'("QMMM:")')
      write (6,'("QMMM: SCF Energy =",f22.14," KCal/mol, ",f22.14," KJ/mol")') escf, escf*4.184d0
      !If verbosity level is greater than 1 we also print the nuclear and electronic energies.
      if (qmmm_nml%verbosity > 1) then
         write (6,'("QMMM:")')
         write (6,'("QMMM:        Electronic energy = ",f18.8," eV (",f18.8," KCal/mol)")') &
              qmmm_struct%elec_eng, qmmm_struct%elec_eng*EV_TO_KCAL
         if  ( qmmm_nml%qmtheory == DFTB ) then
            write (6,'("QMMM:         Repulsive energy = ",f18.8," eV (",f18.8," KCal/mol)")') &
                 qmmm_struct%enuclr_qmqm,qmmm_struct%enuclr_qmqm*EV_TO_KCAL
!!            write (6,'("QMMM:        Careful: Dispersion Energy Already Included.")')
!!            write (6,'("QMMM:        Dispersion Energy = ",f18.8," eV (",f18.8," KCal/mol)")') &
!!                 dftb_edisp*AU_TO_KCAL
            total_energy = qmmm_struct%elec_eng + qmmm_struct%enuclr_qmqm
         else
           write (6,'("QMMM: QM core - QM core energy = ",f18.8," eV (",f18.8," KCal/mol)")') &
                  qmmm_struct%enuclr_qmqm,qmmm_struct%enuclr_qmqm*EV_TO_KCAL
           write (6,'("QMMM: QM core - MM atom energy = ",f18.8," eV (",f18.8," KCal/mol)")') &
                  qmmm_struct%enuclr_qmmm,qmmm_struct%enuclr_qmmm*EV_TO_KCAL
           write (6,'("QMMM: Total core - core energy = ",f18.8," eV (",f18.8," KCal/mol)")') &
                  qmmm_struct%enuclr_qmmm+qmmm_struct%enuclr_qmqm, &
                  (qmmm_struct%enuclr_qmmm+qmmm_struct%enuclr_qmqm)*EV_TO_KCAL
           total_energy = qmmm_struct%elec_eng + qmmm_struct%enuclr_qmmm + qmmm_struct%enuclr_qmqm
         end if
         write (6,'("QMMM:             Total energy = ",f18.8," eV (",f18.8," KCal/mol)")') &
                total_energy, total_energy*EV_TO_KCAL
         !If verbosity level is greater than 3 we also print the force array on the QM atoms
         if (qmmm_nml%verbosity > 3) then
            write (6,'("QMMM:")')
            write (6,'("QMMM: Forces on QM atoms from SCF calculation")')
            write (6,'("QMMM: Atm ",i6,": ",3f20.14)') (j,qmmm_struct%dxyzqm(1,j), qmmm_struct%dxyzqm(2,j), &
                   qmmm_struct%dxyzqm(3,j), j=1,qmmm_struct%nquant_nlink)
            if (qmmm_nml%verbosity > 4) then
               !Also print info in KJ/mol
               write (6,'("QMMM:")')
               write (6,'("QMMM: Forces on QM atoms from SCF calculation (KJ/mol)")')
               write (6,'("QMMM: Atm ",i6,": ",3f20.14)') (j,qmmm_struct%dxyzqm(1,j)*4.184d0, &
                      qmmm_struct%dxyzqm(2,j)*4.184d0, qmmm_struct%dxyzqm(3,j)*4.184d0, &
                      j=1,qmmm_struct%nquant_nlink)
            end if   
         end if
      end if
   end if

   qmmm_struct%qm_mm_first_call = .false.

   return

end subroutine sqm_energy

!======================END OF QM_MM ======================================

!-------------------------------------------------
!     --- FLOAT_LEGAL_RANGE ---
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Check the range of a float; abort on illegal values.
subroutine float_legal_range(string,param,lo,hi)
   implicit none
   _REAL_ param,lo,hi
   character(len=*)string

   if ( param < lo .or. param > hi )then
      write(6,59)
      write(6,60)string,param
      write(6,61)
      write(6,62)lo,hi
      write(6,63)
      call mexit(6,1)
   end if
   59 format(/,1x,'Ewald PARAMETER RANGE CHECKING: ')
   60 format(1x,'parameter ',a,' has value ',e12.5)
   61 format(1x,'This is outside the legal range')
   62 format(1x,'Lower limit: ',e12.5,' Upper limit: ',e12.5)
   63 format(1x,'Check ew_legal.h')
   return
end subroutine float_legal_range 

!-------------------------------------------------
!     --- INT_LEGAL_RANGE ---
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Check the range of an integer; abort on illegal values.
subroutine int_legal_range(string,param,lo,hi)
   implicit none
   integer param,lo,hi
   character(len=*)string

   if ( param < lo .or. param > hi )then
      write(6,59)
      write(6,60)string,param
      write(6,61)
      write(6,62)lo,hi
      write(6,63)
      call mexit(6,1)
   end if
   59 format(/,1x,'PARAMETER RANGE CHECKING: ')
   60 format(1x,'parameter ',a,' has value ',i8)
   61 format(1x,'This is outside the legal range')
   62 format(1x,'Lower limit: ',i8,' Upper limit: ',i8)
   63 format(1x,'The limits may be adjustable; search in the .h files ')
   return
end subroutine int_legal_range 

!-------------------------------------------------
!     --- SANDER_BOMB ---
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Print an error message and quit
subroutine sander_bomb(routine,string1,string2)
   implicit none
   character(len=*) routine,string1,string2

   write(6, '(1x,2a)') &
         'SANDER BOMB in subroutine ', routine
   write(6, '(1x,a)') string1
   write(6, '(1x,a)') string2
   call mexit(6,1)
end subroutine sander_bomb
!-------------------------------------------------

!     --- CROSS ---


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine cross here]
subroutine cross(v1,v2,v12)
   
   !    v12 is cross product of v1 and v2
   
   _REAL_ v1(3),v2(3),v12(3)
   v12(1) = v1(2)*v2(3)-v1(3)*v2(2)
   v12(2) = v1(3)*v2(1)-v1(1)*v2(3)
   v12(3) = v1(1)*v2(2)-v1(2)*v2(1)
   return
end subroutine cross 
subroutine getsqmx(natom,x,atnam,atnum)
   
   !     --- reads initial coords,

   implicit none
   _REAL_ x(*)
   integer i,j,i3,lun
   integer natom,ier,atnum(*)
   character(len=8) atnam(*)
   character(len=80) line

   lun = 5 
   !  skip over the &qmmm namelist at the beginning:
   do i=1,20
      read(5,'(a)') line
      if( line(1:2) == " /" ) go to 9
   end do
   write(0,*) 'Error in finding end of qmmm namelist'
   call mexit(6,1)
   
   9 i3=0
   do i=1,999
      read(lun,*, end=10, err=11) atnum(i),atnam(i),x(i3+1),x(i3+2),x(i3+3)
      i3 = i3+3
   end do
   10 natom = i-1
   return

   11 write(0,*) 'error in reading coordinates '
   call mexit(6,1)

end subroutine getsqmx 

!  following stub routine to avoid changing qmmm_modele.f for sqm:
#ifdef MPI
subroutine qm2_variable_solv_mpi_setup
   return
end subroutine qm2_variable_solv_mpi_setup
#endif
