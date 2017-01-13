! Subroutine to calculate the electrostatic interactions between atoms on different molecules
!   within a given cuttoff radius.  This is the *DAMPED SHIFTED FORCE* Coulomb's law
!
!         V_c(r_ij) = C*q(i)*q(j)/r_ij
!
!         V_ij(r_ij) = V_c(r_ij) - V_c(r_cut) - (dV/dr)|_r_cut * (r_ij - r_cut)
!
! This routine assumes an array of the charge products qq(i,j) = q(i)*q(j) has been calculated
!
! The units are kcal/mol (energy), electron charge, e (charge), and Angstroms (distance)
!
! Jan. 12, 2017 - Paul Burris & Ward Thompson
!

  Subroutine coulomb_dsf

  use common_variables
  implicit none

!Working variables
  integer i, j
  real(kind=dp), parameter :: C_coul = 332.0637301_dp
  real(kind=dp) :: rij, fxtmp, fytmp, fztmp
  real(kind=dp) :: Cpref, etmp, ecut, exptmp, expcut
  real(kind=dp) :: pi = 4.0_dp*atan(1.0_dp)

!External functions
  real(kind=dp), external :: erfc

! Loop over pairs of atoms

  v_c = 0.0_dp; fx_c = 0.0_dp; fy_c = 0.0_dp; fz_c = 0.0_dp

  do j = 1, n_atoms - 1
     do i = j + 1, n_atoms

        rij = sqrt( rx(i,j)**2 + ry(i,j)**2 + rz(i,j)**2 )  ! Calculate the atom-atom distance

        if (rij.le.r_cut) then
           
           !Define some quantities used multiple times
           Cpref = C_coul*qq(i,j)
           etmp = erfc(alpha*rij)/rij**2
           ecut = erfc(alpha*r_cut)/r_cut**2
           exptmp = 2.0_dp*alpha/sqrt(pi)*exp(-(alpha*rij)**2)/rij
           expcut = 2.0_dp*alpha/sqrt(pi)*exp(-(alpha*r_cut)**2)/r_cut

           v_c = v_c + Cpref*( etmp*rij - ecut*r_cut + ( ecut + expcut*(rij-r_cut) ) )  ! potential


           fxtmp = ( Cpref*rx(i,j)/rij )*( (etmp + exptmp) - (ecut + expcut) )  !force on i in x due to j
           fytmp = ( Cpref*ry(i,j)/rij )*( (etmp + exptmp) - (ecut + expcut) )  
           fztmp = ( Cpref*rz(i,j)/rij )*( (etmp + exptmp) - (ecut + expcut) )

           fx_c(i) = fx_c(i) + fxtmp 
           fy_c(i) = fy_c(i) + fytmp
           fz_c(i) = fz_c(i) + fztmp 

           fx_c(j) = fx_c(j) - fxtmp 
           fy_c(j) = fy_c(j) - fytmp
           fz_c(j) = fz_c(j) - fztmp 

        endif

     enddo
  enddo

 End Subroutine coulomb_dsf

