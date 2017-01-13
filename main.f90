program classical_md
!************************************************************************
!
! Fortran program for running a classical molecular dynamics simulation.
! Use velocity verlet integrator.
!
! Built by the most awesomest group of scientists.
!
! Jan. 2017
!
!************************************************************************

       use common_variables

       implicit none

       real(kind = dp) :: temp,dt

       integer(kind = ip) :: nmd,nlog,df_xyz,df_thermo,df_rest

       character(len=50) :: data_filename

       logical :: restart

! read simulation input parameters
! Note: this subroutine also calls read_data.f90 to read in initial
! configurations

       call read_input(data_filename,df_xyz,df_thermo,df_rest,nvt_type)

! set initial velocites if not a restart

       if (.not. restart) call initvel

! calculate force of initial configuration

       call forces

! do a force check

       call force_check

! open output files

       call output_setup('start')

!*****************************START MD SIMULATION***********************!

       do i = 1 , nmd

! update positions in first stage of VV integrator

            call integrator(1)

! calculate forces of new positions

            call forces

! update velocities in second stage of VV integrator

            call integrator(2)

! calculate thermodynamic properties

            if (mod(i,nlog) .eq. 0) call thermo_stuff()

! apply thermostat

            if (thermostat .gt. 1) then

                 call thermostat(kb,KE,temp,temp_c,lambda)

            endif

! write traj

            if (mod(i,nhis) .eq. 0) call thermo_dump(i)

       enddo

! close output files

       call output_setup('end')

end program classical_md
