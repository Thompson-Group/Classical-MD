subroutine calculate_bond(rx, ry, rz, BT_1, BT_2)
    use kinds
    implicit none


    IF (bond_style.eq."harmonic") THEN
        call harmonic_bond_energy()
        call harmonic_bond_force()
    END IF


end subroutine
