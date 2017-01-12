subroutine read_data(data_filename, n_atoms, n_bonds, n_angles, n_dih, n_a_type, n_b_type, n_dih_type, &
        n_imp_type, n_imp, bond_table, angle_table, dih_table,restart, vx, vy ,vz, xlo, xhi,&
        ylo, yhi, zlo, zhi, a_id, mol_id, a_type, q, x, y, z) 
    use kinds
    implicit none
    character(len=50) :: ignore, data_filename
    integer(kind=ip) :: n_atoms, n_bonds, n_angles, n_dih, n_imp, n_a_type, n_b_type, n_ang_type, n_dih_type, n_imp_type
    integer(kind=ip) :: i
    real(kind=dp) :: xlo, xhi, ylo, yhi, zlo, zhi
    real(kind=dp), allocatable,dimension(:):: a_id, mol_id, a_type, q, x, y, z, vx, vy, vz
    real(kind=dp), allocatable,dimension(:,:) :: bond_table, angle_table, dih_table
    logical :: restart
    restart = .FALSE.
    !open(unit=30, file=data_filename)
    open(unit=30, file="test.in")
    read(30,*)
    read(30,*)
    read(30,*) n_atoms, ignore
    read(30,*) n_bonds, ignore
    read(30,*) n_angles, ignore
    read(30,*) n_dih, ignore
    read(30,*) n_imp, ignore
    allocate(a_id(n_atoms))
    allocate(mol_id(n_atoms))
    allocate(a_type(n_atoms))
    allocate(q(n_atoms))
    allocate(x(n_atoms))
    allocate(y(n_atoms))
    allocate(z(n_atoms))
    allocate(vx(n_atoms))
    allocate(vy(n_atoms))
    allocate(vz(n_atoms))
    allocate(bond_table(n_bonds,3))
    allocate(angle_table(n_angles,4))
    allocate(dih_table(n_dih, 5))
    
    read(30,*)
    read(30,*) n_a_type, ignore, ignore
    IF (n_bonds.ne.0) THEN
        read(30,*) n_b_type, ignore, ignore
    END IF
    IF (n_angles.ne.0) THEN
        read(30,*) n_ang_type, ignore, ignore
    END IF
    IF (n_dih.ne.0) THEN
        read(30,*) n_dih_type, ignore, ignore
    END IF
    IF (n_imp.ne.0) THEN
        read(30,*) n_imp_type, ignore, ignore
    END IF
    read(30,*)
    read(30,*) xlo, xhi, ignore, ignore
    read(30,*) ylo, yhi, ignore, ignore
    read(30,*) zlo, zhi, ignore, ignore
    read(30,*)
    !Reads in the Atoms Section
    read(30,*)
    read(30,*)
    DO i = 1, n_atoms
        read(30,*) a_id(i), mol_id(i), a_type(i), q(i), x(i), y(i), z(i)
    END DO
    IF (n_bonds.ne.0) THEN
        read(30,*)
        read(30,*)
        read(30,*)
        DO i = 1, n_bonds
            read(30,*) ignore, bond_table(i,1), bond_table(i,2), bond_table(i,3)
        END DO
    END IF
    IF (n_angles.ne.0) THEN
        read(30,*)
        read(30,*)
        read(30,*)
        DO i = 1, n_angles
            read(30,*) ignore, angle_table(i,1), angle_table(i,2), angle_table(i,3), angle_table(i,4)
        END DO
    END IF
    IF (n_dih.ne.0) THEN
        read(30,*)
        read(30,*)
        read(30,*)
        DO i = 1, n_dih
            read(30,*) ignore, dih_table(i,1), dih_table(i,2), dih_table(i,3), dih_table(i,4), dih_table(i,5)
        END DO
    END IF
    IF (n_imp.ne.0) THEN
        write(*,*) "Error: Identified impropers, when code does not support them."
    END IF
    IF (restart) THEN
        DO i = 1, n_atoms
            read(30,*) ignore, vx(i), vy(i), vz(i)
        END DO

    END IF
    close(unit=30)
end subroutine
