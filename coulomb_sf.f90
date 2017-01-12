! Subroutine to calculate the electrostatic interactions between atoms on different molecules
!   within a given cuttoff radius.  This is the *SHIFTED FORCE* Coulomb's law
!
!         V_c(r_ij) = C*q(i)*q(j)/r_ij
!
!         V_ij(r_ij) = V_c(r_ij) - V_c(rcut) - (dV/dr)|_rcut * (r_ij - rcut)
!
! This routine assumes an array of the charge products qq(i,j) = q(i)*q(j) has been calculated
!
! The units are kcal/mol (energy), electron charge, e (charge), and Angstroms (distance)
!
! Jan. 12, 2017 - Paul Burris & Ward Thompson
!

  Subroutine coulomb_sf(natoms,rx,ry,rz,rcut,qq, fx_c,fy_c,fz_c,v_c)

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
  real(kind=dp) :: rij, fxtmp, fytmp, fztmp

! Loop over pairs of atoms

  v_c = 0.0_dp; fx_c = 0.0_dp; fy_c = 0.0_dp; fz_c = 0.0_dp

  do j = 1, natoms - 1
     do i = j + 1, natoms

        rij = sqrt( rx(i,j)**2 + ry(i,j)**2 + rz(i,j)**2 )  ! Calculate the atom-atom distance

        if (rij.le.rcut) then
           
           v_c = v_c + C_coul*qq(i,j)*(1.0_dp/rij - 1.0_dp/rcut + (rij - rcut)/rcut**2)  ! Calculate the contribution to the potential

           fxtmp = C_coul*qq(i,j)*rx(i,j)*(1.0_dp/rij**2 - 1.0_dp/rcut**2)/rij
           fytmp = C_coul*qq(i,j)*ry(i,j)*(1.0_dp/rij**2 - 1.0_dp/rcut**2)/rij
           fztmp = C_coul*qq(i,j)*rz(i,j)*(1.0_dp/rij**2 - 1.0_dp/rcut**2)/rij

           fx_c(i) = fx_c(i) + fxtmp 
           fy_c(i) = fy_c(i) + fytmp
           fz_c(i) = fz_c(i) + fztmp

           fx_c(j) = fx_c(j) - fxtmp 
           fy_c(j) = fy_c(j) - fytmp
           fz_c(j) = fz_c(j) - fztmp

        endif

     enddo
  enddo

 End Subroutine coulomb_sf

