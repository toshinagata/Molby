! <compile=optimized>
#include "copyright.h"
#include "../include/dprec.fh"

!Subroutines for the printing of charges on QM atoms.
!Author Ross Walker (SDSC, 2009)

subroutine qm2_print_charges(nstep,dftb_chg,nquant_nlink,scf_mchg,iqm_atomic_numbers)

  use qmmm_module, only : element_sym

  implicit none

!Passed in
  integer, intent(in) :: nstep, dftb_chg, nquant_nlink
  integer, intent(in) :: iqm_atomic_numbers(nquant_nlink)
  _REAL_, intent(in) :: scf_mchg(nquant_nlink)

!Local
  integer :: i
  _REAL_ :: total_mulliken_charge, total_cm3_chg
  _REAL_ :: scf_cm3(nquant_nlink)
  _REAL_ :: mulliken_charge


  total_mulliken_charge=0.0d0
  total_cm3_chg=0.d0
  write(6,'("    Atomic Charges for Step",i8," :")') nstep
  if (dftb_chg == 1) then
     write(6,'("  Atom    Element       Mulliken Charge       CM3 Charge")')
     call qm2_dftb_cm3(scf_mchg, scf_cm3)
  else
     write(6,'("  Atom    Element       Mulliken Charge")')
  end if

  do i=1,nquant_nlink
     !Mulliken charges have already been calculated and stored.
     mulliken_charge = scf_mchg(i)

     total_mulliken_charge=total_mulliken_charge+mulliken_charge
     if (dftb_chg == 1) then
        total_cm3_chg =  total_cm3_chg + scf_cm3(i)
        write(6,'(" ",i5,"      ",A2,"        ",F14.3,"    ",F14.3)') i, &
        element_sym(iqm_atomic_numbers(i)), &
        mulliken_charge, scf_cm3(i)
     else
        write(6,'(" ",i5,"      ",A2,"        ",F14.3)') i, &
        element_sym(iqm_atomic_numbers(i)), &
        mulliken_charge
     endif
  end do

  if (dftb_chg == 1) then
     write(6,'(" Total Charges: ",8X,F12.3,6X,F12.3)')  &
          total_mulliken_charge, total_cm3_chg
  else
     write(6,'(" Total Mulliken Charge =",F12.3)') &
          total_mulliken_charge
  endif

  return

end subroutine qm2_print_charges

