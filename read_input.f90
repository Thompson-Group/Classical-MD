subroutine read_input(data_filename,df_xyz,df_thermo,df_rest)
    use kinds
    use common_variables
    integer(kind=ip) :: i, j, ign_int,df_xyz, df_thermo, df_rest
    character(len=50) :: data_filename,word,input_filename, run_style, nvt_type
    logical :: restart
    open(unit=ninput, file=input_filename)
     

    read(ninput,*) word 
    if(word .eq. 'read_data') then
       read(ninput,*) data_filename
       restart = .false.   
    else if(word .eq. 'read_restart') then
       read(ninput,*) data_filename
       restart = .true. 
    else
       write(*,*) "Error : Data file not found at top of input file"
    end if

    call read_data(data_filename,restart)

        read(ninput,*)
    do
        read(ninput,*) word

        if (word .eq. 'bond_coeffs') then
            do i = 1, n_b_type
            ! Provide the bond coeffs in bond type order            
                read(ninput,*)  ign_int , k_r(i), req(i)
            enddo

        else if (word .eq. 'angle_coeffs') then 
            do i = 1, n_angle_type
            ! Provide the angle coeffs in angle type order  
                read(ninput,*) ign_int, k_ang(i), theta_eq(i)
            enddo
        else if (word .eq. 'cutoff_distance') then
            read(ninput,*) r_cut
        !coul_flag = 0: off, =1: shifted force, =2: damped shifted force
        else if (word .eq. 'coulomb_style') then 
            read(ninput,*) alpha, coul_flag 
        else if (word .eq. 'force_check') then
            read(ninput,*) fc_flag,del
        else if (word .eq. 'run_style') then 
            read(ninput,*)  run_style

            if (run_style .eq. 'nve') then
                nvt_type='none'
            else if (run_style .eq. 'nvt') then
                !If run_style is nvt, need to choose options: rescale or anderson 
                read(ninput,*) nvt_type
            else
                write(*,*) "Error: Couldn't find the run_style in input file"
            endif 
        else if (word .eq. 'timestep') then
            read(ninput,*) dt

            
        else if (word .eq. 'run') then 
            read(ninput,*) nstep

        else if (word .eq. 'dump_freq') then 
            read(ninput,*) df_xyz, df_thermo, df_rest

    
        else if (word .eq. 'temperature') then 
            read(ninput,*) temp
        else
            exit
        end if
        read(ninput,*)  
    end do

       

end subroutine
