! <compile=optimized>
#include "copyright.h"
#include "../include/dprec.fh"
subroutine  qm2_hcore_qmmm(H,ENUCLR,qm_xcrd)

! This routine is repsonsible for calculating the QM-MM interactions between
! QM-MM pairs that are in the QMMM pair list. This routine calls qm2_rotate_qmmm.
! It is responsible for updating the one-electron matix H with the QM electron
! MM core interaction. It also calculates the QM-MM core-core interactions that
! are returned as an addition to ENUCLR.
!
! QM core charge interacts with MM resp charge. The MM resp charge, in electron
! units instead of AMBER's internal units is stored as the 4th index of the qm_COORD array.
! Hence the qm_coord array runs: x1, y1, z1, q1, x2, y2, z2, q2...

! Current code: Ross Walker (TSRI, 2005)

      use qmmm_module, only : qmmm_struct, qm2_struct, qm2_params, qmmm_nml, qmmm_mpi
      implicit none

!Passed in
      _REAL_, intent(inout) :: H(qm2_struct%matsize)
      _REAL_, intent(out) :: enuclr
      _REAL_ , intent(in) :: qm_xcrd(4,qmmm_struct%qm_mm_pairs)

!Local
! Keeps track of the nquant * ni_mm loop iterations
      _REAL_ qm_atom_coord(3)
      integer loop_count, i, j, ia, ib, i2, i1, j1, ii
      logical heavy_atom
      _REAL_ enucij
      _REAL_ E1B(10), E1B_light

!     Loop over REAL QM atoms. - no need to worry about link atoms here
      loop_count=0
!      do i=1,qmmm_struct%nquant_nlink
      do i=qmmm_mpi%nquant_nlink_start,qmmm_mpi%nquant_nlink_end
!
!       Definitions:
!
!       orb_loc(1,i) = first atomic orbital on atom i
!       orb_loc(2,i) =  last atomic orbital on atom i
! iqm_atomic_numbers(i) =    atomic number for atom i
!
        heavy_atom = (qm2_params%natomic_orbs(i) > 1)
        ia = qm2_params%orb_loc(1,i)
        qm_atom_coord(1) = qmmm_struct%qm_coords(1,i)
        qm_atom_coord(2) = qmmm_struct%qm_coords(2,i)
        qm_atom_coord(3) = qmmm_struct%qm_coords(3,i)

!       Loop over MM atoms that interact with QM atom i, update the
!       one-electron matrix H, and accumulate core-core repulsions
!       ENUCLR.

! Loop in steps of since qm_xcrd array is x,y,z,chg,x,y,z,chg...
! Duplicated code here but it avoids the if(heavy_atom) inside the loop
        if (heavy_atom) then
          ib = qm2_params%orb_loc(2,i) !Will be same as ia for light atom
          do j=1,qmmm_struct%qm_mm_pairs
            loop_count = loop_count+1
!          Get the QM-MM interactions:
!          QM electrons - MM core ---> E1B
!          QM core - MM core ----> enucij
            call qm2_rotate_qmmm_heavy(loop_count,i,qm_atom_coord,qm_xcrd(1,j),E1B,enucij)
!          Add E1B to the corresponding diagonal elements of the one-electron
!          matrix H.
            i2 = 0
            do i1=ia,ib
              ii = qm2_params%pascal_tri1(i1) + ia - 1
              do j1=ia,i1
                ii = ii + 1
                i2 = i2 + 1
                H(ii) = H(ii) + E1B(i2)
              end do 
            end do 
!         Add on QM core - MM core interactions.
            ENUCLR = ENUCLR + enucij
          end do
        else !if (heavy_atom)
          !It is a light atom (Hydrogen) - same notes as above for heavy atom
          do j=1,qmmm_struct%qm_mm_pairs
            loop_count = loop_count+1
            call qm2_rotate_qmmm_light(loop_count,i,qm_atom_coord,qm_xcrd(1,j),E1B_light,enucij)
            ii = qm2_params%pascal_tri1(ia) + ia 
            H(ii) = H(ii) + E1B_light
            ENUCLR = ENUCLR + enucij
          end do
        end if !if (heavy_atom)
      end do
      RETURN
end subroutine qm2_hcore_qmmm

SUBROUTINE qm2_rotate_qmmm_light(loop_count,IQM,xyz_qm,xyz_mm,E1B,ENUC)

!For light atoms
!-----------------------------------
!Written by Ross Walker (TSRI, 2005)
!-----------------------------------
!See heavy routine below for notes.

      use qmmm_module, only : qmmm_nml, qm2_struct, qmmm_struct, qm2_params, qm2_rij_eqns,alph_mm, &
                              EXPONENTIAL_CUTOFF, PM6
      use constants, only : A2_TO_BOHRS2, AU_TO_EV
      implicit none
! Passed in
      integer, intent(in) :: loop_count, iqm
      _REAL_, intent(in) :: xyz_qm(3), xyz_mm(4)
      _REAL_, intent(out) :: e1b, enuc

! Local
      _REAL_ :: RI
      _REAL_ :: oneBDDi1
      _REAL_ :: X1, X2, X3
      _REAL_ :: RR2, oneRIJ, RIJ, RIJ2
      _REAL_ :: scale, anam1, temp_real, temp_real2
      integer :: qmitype, i
    
#include "qm2_array_locations.h"
      if (qmmm_nml%qmtheory==PM6) then
        ! Pairwise core core
        call sander_bomb('qm2_hcore_qmmm','ERROR, PM6 Not currently supported.','')
      else
        if (qmmm_nml%qmmmrij_incore) then
          scale=qm2_rij_eqns%qmmmrijdata(QMMMEXP1,loop_count)+qm2_rij_eqns%qmmmrijdata(QMMMEXP2,loop_count)
          RI = AU_TO_EV*qm2_rij_eqns%qmmmrijdata(QMMMSQRTAEE,loop_count)
        else
          X1 = xyz_qm(1) - xyz_mm(1)
          X2 = xyz_qm(2) - xyz_mm(2)
          X3 = xyz_qm(3) - xyz_mm(3)
          RIJ2=X1*X1+X2*X2+X3*X3
          RR2=RIJ2*A2_TO_BOHRS2
          oneRIJ=1.0d0/sqrt(RIJ2)
          RIJ =RIJ2*oneRIJ
          scale=EXP(-qm2_params%cc_exp_params(iqm)*RIJ)+EXP(-ALPH_MM*RIJ)
          oneBDDi1=qm2_params%multip_2c_elec_params(3,iqm)**2
          RI = AU_TO_EV/sqrt(RR2+oneBDDi1)
        end if

        E1B=-xyz_mm(4)*RI
        ENUC = -qm2_params%core_chg(IQM)*E1B
        scale = ABS(scale*ENUC)

      end if

!Add in the extra core core repulsion terms if qmmm_int==2. This is equivalent of the
!way CHARMM and DYNAMO do QMMM interactions.
      IF(qmmm_nml%qmmm_int == 2) then
         if (qmmm_struct%AM1_OR_PM3) THEN
!          Add gaussians.
           if (qmmm_nml%qmmmrij_incore) then
             rij = qm2_rij_eqns%qmmmrijdata(QMMMRIJ,loop_count)
             onerij = qm2_rij_eqns%qmmmrijdata(QMMMONERIJ,loop_count)
           end if
           qmitype = qmmm_struct%qm_atom_type(iqm)
           anam1=0.0d0
           do I=1,qm2_params%num_fn(qmitype)
              temp_real=RIJ-qm2_params%FN3(i,qmitype)
              temp_real2=qm2_params%FN2(i,qmitype)*temp_real*temp_real
              if (temp_real2 < EXPONENTIAL_CUTOFF) then!Skip doing the exponential if it is essentially zero
                anam1=anam1+qm2_params%FN1(i,qmitype)*EXP(-temp_real2)
              end if
           end do
           anam1=anam1*qm2_params%core_chg(IQM)*xyz_mm(4)*oneRIJ
           scale = scale + anam1
         end if
      ENDIF

      ENUC = ENUC+scale !E(AB) = ZaZb(SASA,SBSB)(1+exp(-alphaA*RAB)+exp(-alphaB*RAB))
    return
end subroutine qm2_rotate_qmmm_light

SUBROUTINE qm2_rotate_qmmm_heavy(loop_count,IQM,xyz_qm,xyz_mm,E1B,ENUC)

!For heavy atoms
!-----------------------------------
!Written by Ross Walker (TSRI, 2005)
!-----------------------------------

!-----------------------------------
!This routine calculates the QMelectron - MMcore interactions and the
!QMcore-MMcore interaction for a give QM-MM pair. Initially calculated in
!a local frame by qm2_repp_qmmm they are subsequently rotated into the
!molecular frame by this routine. The core charge for the MM atom is the
!RESP charge from the prmtop file in electron units.

!On Input:
! loop_count - Offset into certain arrays
!        iqm - Current qm atom in our loop from 1 to nquant
!     xyz_qm - Cartesian coordinates of the QM atom.
!     xyz_mm - Cartesian coordinates and charge (4) of MM atom.
!        E1B - QM-MM electron core interaction in KCal/mol for each
!              of the 10 multipoles of the QM atom. (s,s), (s,px), (s,py)
!              (s,pz), (px,px), (px,py), (px,pz), (py,py), (py,pz), (pz,pz)
!       ENUC - QM-MM core core interaction in KCal/mol.
!-----------------------------------

      use qmmm_module, only : qmmm_nml, qmmm_struct, qm2_struct, qm2_params, qm2_rij_eqns, &
                              alph_mm, AXIS_TOL, EXPONENTIAL_CUTOFF
      use constants, only : A_TO_BOHRS, A2_TO_BOHRS2, AU_TO_EV, HALF_AU_TO_EV, FOURTH_AU_TO_EV
                            
      implicit none
! Passed in
      integer, intent(in) :: loop_count, iqm
      _REAL_, intent(in) :: xyz_qm(3), xyz_mm(4)
      _REAL_, intent(out) :: e1b(10), enuc

! Local
      _REAL_ RI(4)
      _REAL_ :: oneBDDi1, oneBDDi2, oneBDDi3
      _REAL_ :: sqrtrr2aqe, qmi_DD, qmi_QQ
      _REAL_ :: X1, X2, X3,Y1, Y2, Y3,Z1, Z2,Z3, oneZ3, X3_2
      _REAL_ :: chrgmm, RR2, oneRIJ, RIJ, RR, RIJ2, CHGMM_RI2, CHGMM_RI3, CHGMM_RI4
      _REAL_ :: scale, anam1, temp_real, temp_real2
      integer :: qmitype, i
      
      X1 = xyz_qm(1) - xyz_mm(1)
      X2 = xyz_qm(2) - xyz_mm(2)
      X3 = xyz_qm(3) - xyz_mm(3)
      chrgmm = xyz_mm(4)
#include "qm2_array_locations.h"

!---RIJ and related equation storage specifics---
      if (qmmm_nml%qmmmrij_incore) then
        scale=qm2_rij_eqns%qmmmrijdata(QMMMEXP1,loop_count)+qm2_rij_eqns%qmmmrijdata(QMMMEXP2,loop_count)
        RI(1) = AU_TO_EV*qm2_rij_eqns%qmmmrijdata(QMMMSQRTAEE,loop_count)
        oneRIJ=qm2_rij_eqns%qmmmrijdata(QMMMONERIJ,loop_count)
        RI(2) = HALF_AU_TO_EV*(qm2_rij_eqns%qmmmrijdata(QMMMSQRTRRADDADE,loop_count)- &
                               qm2_rij_eqns%qmmmrijdata(QMMMSQRTRRMDDADE,loop_count))

        RI(3) = RI(1) + FOURTH_AU_TO_EV*(qm2_rij_eqns%qmmmrijdata(QMMMSQRTRRAQQAQE,loop_count)+ &
                qm2_rij_eqns%qmmmrijdata(QMMMSQRTRRMQQAQE,loop_count))- &
                HALF_AU_TO_EV*qm2_rij_eqns%qmmmrijdata(QMMMSQRTRR2AQE,loop_count)

        RI(4) = RI(1) + HALF_AU_TO_EV*(qm2_rij_eqns%qmmmrijdata(QMMMSQRTRRAQQ2AQE,loop_count)- &
                                       qm2_rij_eqns%qmmmrijdata(QMMMSQRTRR2AQE,loop_count))
      else
        RIJ2=X1*X1+X2*X2+X3*X3
        RR2=RIJ2*A2_TO_BOHRS2
        oneRIJ=1.0d0/sqrt(RIJ2)
        RIJ =RIJ2*oneRIJ
        RR = RIJ*A_TO_BOHRS
        scale=EXP(-qm2_params%cc_exp_params(iqm)*RIJ)+EXP(-ALPH_MM*RIJ)
        oneBDDi1=qm2_params%multip_2c_elec_params(3,iqm)**2
        RI(1) = AU_TO_EV/sqrt(RR2+oneBDDi1)
        qmi_DD=qm2_params%multip_2c_elec_params(1,iqm)
        qmi_QQ=2.0d0*qm2_params%multip_2c_elec_params(2,iqm)
        oneBDDi2=qm2_params%multip_2c_elec_params(4,iqm)**2
        oneBDDi3=qm2_params%multip_2c_elec_params(5,iqm)**2
        RI(2)= HALF_AU_TO_EV*(1.0d0/SQRT((RR+qmi_DD)**2+oneBDDi2) - 1.0d0/SQRT((RR-qmi_DD)**2+oneBDDi2))

        SQRTRR2AQE=1.0d0/SQRT(RR2+oneBDDi3)
        RI(3)= RI(1) + FOURTH_AU_TO_EV*(1.0d0/SQRT((RR+qmi_QQ)**2+oneBDDi3) + 1.0d0/SQRT((RR-qmi_QQ)**2+oneBDDi3)) &
               - HALF_AU_TO_EV*SQRTRR2AQE

        RI(4) = RI(1) + HALF_AU_TO_EV*(1.0d0/SQRT(RR2+(qmi_QQ*qmi_QQ)+oneBDDi3) - SQRTRR2AQE)
      end if
!---END RIJ and related equation storage specifics---

! Step 1 is to calculate the interaction between the QM electrons and the MM
! atomic charge for this QM-MM pair in the local frame. MM charge is the RESP
! charge which is currently stored in the 4th index of the MM cartesian coordinate
! array. The reason the charge is packed into the 4th index of the coordinate array
! is so that we move linearly through memory = more cache hits. (RCW). - DONE ABOVE FOR SPEED

!All atoms have S orbitals
      E1B(1)=-chrgmm*RI(1)

!Calculate QM-QM nuclear term
!This is calculates as the QM core charge * chrgmm * RI(1)
!Reuse E1B term to save doing the multiplication of chrgmm*RI(1) again.
      ENUC = -qm2_params%core_chg(IQM)*E1B(1)
      scale = ABS(scale*ENUC)

!Add in the extra core core repulsion terms if qmmm_int==2. This is equivalent of the
!way CHARMM and DYNAMO do QMMM interactions.
      if(qmmm_nml%qmmm_int == 2 .and. qmmm_struct%AM1_OR_PM3) THEN
!        Add gaussians.
         if (qmmm_nml%qmmmrij_incore) then
           rij = qm2_rij_eqns%qmmmrijdata(QMMMRIJ,loop_count)
         end if
         qmitype = qmmm_struct%qm_atom_type(iqm)
         anam1=0.0d0
         do I=1,qm2_params%num_fn(qmitype)
            temp_real=RIJ-qm2_params%FN3(i,qmitype)
            temp_real2=qm2_params%FN2(i,qmitype)*temp_real*temp_real
            if (temp_real2 < EXPONENTIAL_CUTOFF) then!Skip doing the exponential if it is essentially zero
              anam1=anam1+qm2_params%FN1(i,qmitype)*EXP(-temp_real2)
            end if
         end do
         anam1=anam1*qm2_params%core_chg(IQM)*xyz_mm(4)*oneRIJ
         scale = scale + anam1
      end if
      ENUC = ENUC + scale

      X1 = X1*oneRIJ
      X2 = X2*oneRIJ
      X3 = X3*oneRIJ
      X3_2 = X3*X3

      if (abs(X3) > (1.0D0-AXIS_TOL)) then
         X3 = SIGN(1.0D0,X3)
         Y1 = 0.0D0
         Y2 = 1.0D0
         Y3 = 0.0D0
         Z1 = 1.0D0
         Z2 = 0.0D0
         Z3 = 0.0D0
      else
         Z3=SQRT(1.0D0-X3_2)
         oneZ3=1.0D0/Z3
         Y1=-oneZ3*X2*SIGN(1.D0,X1)
         Y2=ABS(oneZ3*X1)
         Y3=0.0D0
         Z1=-oneZ3*X1*X3
         Z2=-oneZ3*X2*X3
      endif

      CHGMM_RI2=-chrgmm*RI(2)
      CHGMM_RI3=-chrgmm*RI(3)
      CHGMM_RI4=-chrgmm*RI(4)

      E1B(2) = CHGMM_RI2*X1
      E1B(3) = CHGMM_RI3*X1*X1+CHGMM_RI4*((Y1*Y1)+(Z1*Z1))
      E1B(4) = CHGMM_RI2*X2
      E1B(5) = CHGMM_RI3*X2*X1+CHGMM_RI4*((Y2*Y1)+(Z2*Z1))
      E1B(6) = CHGMM_RI3*X2*X2+CHGMM_RI4*((Y2*Y2)+(Z2*Z2))
      E1B(7) = CHGMM_RI2*X3
      E1B(8) = CHGMM_RI3*X3*X1+CHGMM_RI4*Z3*Z1
      E1B(9) = CHGMM_RI3*X3*X2+CHGMM_RI4*Z3*Z2
      E1B(10)= CHGMM_RI3*X3_2+CHGMM_RI4*Z3*Z3

    return

end subroutine qm2_rotate_qmmm_heavy

