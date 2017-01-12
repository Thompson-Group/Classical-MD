module kinds

 !******************************************************************************
 !
 ! Fortran module for setting data types.
 ! Author - Jacob A. Harvey
 ! Copyright - CRUNCH Lab Feb 2012 
 !
 !*******************************************************************************
 
   implicit none
 
   integer, parameter :: dp = selected_real_kind(14,300)
   integer, parameter :: ip = selected_int_kind(12)
 
End module kinds


! Subroutine to calculate the electrostatic interactions between atoms on different molecules
!   within a given cuttoff radius.  This is the direct Coulomb's law
!
!         V_ij(r_ij) = C*q(i)*q(j)/r_ij
!
! This routine assumes an array of the charge products qq(i,j) = q(i)*q(j) has been calculated
!
! The units are kcal/mol (energy), electron charge, e (charge), and Angstroms (distance)
!
! Jan. 12, 2017 - Paul Burris & Ward Thompson
!

  Subroutine coulomb(natoms,rx,ry,rz,rcut,qq, fx_c,fy_c,fz_c,v_c)

  use kinds
  implicit none

  integer :: natoms
  real(kind=dp) :: rcut, v_c
  real(kind=dp), dimension(natoms,natoms) :: rx, ry, rz
  real(kind=dp), dimension(natoms) :: fx_c, fy_c, fz_c
  real(kind=dp), dimension(natoms,natoms) :: qq

!Working variables
  integer i, j
  real(kind=dp), parameter :: C_coul = 1.0_dp
  real(kind=dp) :: rij

! Loop over pairs of atoms

  v_c = 0.0_dp; fx_c = 0.0_dp; fy_c = 0.0_dp; fz_c = 0.0_dp

  do j = 1, natoms - 1
     do i = j + 1, natoms

        rij = sqrt( rx(i,j)**2 + ry(i,j)**2 + rz(i,j)**2 )  ! Calculate the atom-atom distance

        if (rij.le.rcut) then
           
           v_c = v_c + C_coul*qq(i,j)/rij                   ! Calculate the contribution to the potential

           fx_c(i) = fx_c(i) + C_coul*qq(i,j)*rx(i,j)/rij**3
           fy_c(i) = fy_c(i) + C_coul*qq(i,j)*ry(i,j)/rij**3
           fz_c(i) = fz_c(i) + C_coul*qq(i,j)*rz(i,j)/rij**3

        endif

     enddo
  enddo

 End Subroutine coulomb  
