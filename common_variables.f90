Module common_variables
!**************************************************************************************************************
! This module contains arrays and scalars that are common to several subroutines in this md code.
! These variables do not need to be declared in other subroutines. 
!Zeke will allocate the arrays and read in values for the variables in the read_data code
! use this module as needed in other subroutines
!Written by Mesele and Pubudu
! Hackathon Thursday January 17, 2017
!**************************************************************************************************************

  use kinds
  Implicit none

  integer(kind=ip) :: n_atoms, n_bonds, n_angles, n_dih, n_imp
  integer(kind=ip) :: n_a_type, n_b_type, n_angle_type, n_dih_type, n_imp_type
  real(kind=dp) :: xlo, xhi, ylo, yhi, zlo, zhi, Lx, Ly, Lz, r_cut, alpha, v_a, v_b, v_c, v_v, v_tot
  integer(kind=ip), allocatable, dimension(:)   :: a_id, mol_id, a_type
  real(kind=dp), allocatable, dimension(:,:) :: rx, ry, rz, bond_table, angle_table, ee, ss, qq, dih_table
  real(kind=dp), allocatable, dimension(:) :: fx_v, fy_v, fz_v, fx_c, fy_c, fz_c, fx_tot, fy_tot, fz_tot
  real(kind=dp), allocatable, dimension(:) :: M, q, vx, vy, vz, ep, sig
  real(kind=dp), allocatable, dimension(:) :: x, y, z, fx_b, fy_b, fz_b, fx_a, fy_a, fz_a, k_r, req, k_ang, theta_eq
  real(kind=dp), allocatable, dimension(:) :: x_o, y_o, z_o, vx_o, vy_o, vz_o, fx_o, fy_o, fz_o
end module
