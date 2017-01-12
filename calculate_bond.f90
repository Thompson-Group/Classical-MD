subroutine calculate_bond(rx, ry, rz, BT_1, BT_2)
    use kinds
    implicit none
    ! Variable declarations
    !Loop Integers, Indices, etc.
    integer(kind=ip) :: n_atom, n_bond, n_b_type, i, b_type, b_atom_1, b_atom_2
    !Bonded Parameter Arrays
    real(kind=dp), allocatable,dimension(n_b_type) :: k, req
    !Arrays for storing information, all of n_bond length.
    real(kind=dp), allocatable,dimension(n_bond):: rx,ry,rz,rinst
    !Arrays for storing information, all of n_atom length.
    real(kind=dp), allocatable,dimension(n_atom) :: fx_b, fy_b, fz_b
    !Array read in from read_data.f90, stores bonding information.
    real(kind=dp), allocatable,dimension(n_bond,4)::bond_table
    real(kind=dp), bond_pe
    allocate(F(Nbond))
    

    !Bond Energy and Force Calculations
    IF (bond_style.eq."harmonic") THEN
        ! Do loop over number of bonds
        DO i = 1, n_bond
            !Index calculations
            b_type  =bond_table(i+n_bond)
            b_atom_1=bond_table(i+2*n_bond)
            b_atom_2=bond_table(i+3*n_bond)
            !Bond Length Calculations
            !rinst is the instantaneous bond length
            !r_sep is rinst-req
            rinst(i)=sqrt((rx(b_atom_1)-rx(b_atom_2))**2 &
                        + (ry(b_atom_1)-ry(b_atom_2))**2 &
                        + (rz(b_atom_1)-rz(b_atom_2))**2 )
            r_sep=rinst(i)-req(b_type)
            !Potential Energy Calculation
            u_bond(i)=0.5*k(b_type)*r_sep**2
            bond_pe = bond_pe + u_bond(i)

            !Force calculation 
            !Calculates the force on atom 1 as:
            !F_1 = -dU/dr * dr/dx
            !Note: du/dr = k*r_sep
            !and dr/dx = x_1/rinst(i)
            !Calculates the force on atom 2 as:
            !F_2 = -F_1

            fx_b(b_atom_1) = - k(b_type) * r_sep * rx(b_atom_1) / rinst(i)
            fx_b(b_atom_2) = - fx_b(b_atom_1)
            fy_b(b_atom_1) = - k(b_type) * r_sep * ry(b_atom_1) / rinst(i)
            fy_b(b_atom_2) = - fy_b(b_atom_1)
            fz_b(b_atom_1) = - k(b_type) * r_sep * rz(b_atom_1) / rinst(i)
            fz_b(b_atom_2) = - fz_b(b_atom_1)
        END DO
    END IF


end subroutine
