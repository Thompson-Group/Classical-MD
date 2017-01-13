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

       integer(kind = ip) :: df_xyz,df_thermo,df_rest,nstep,fc_flag
       integer(kind = ip) :: i,d

       character(len=50) :: data_filename,nvt_type

       logical :: restart,change

! read simulation input parameters
! Note: this subroutine also calls read_data.f90 to read in initial
! configurations

       call read_input(data_filename,df_xyz,df_thermo,df_rest,nvt_type,fc_flag,&
                       nstep)

! set initial velocites if not a restart

       if (.not. restart) call initvel

! calculate force of initial configuration

       call forces

! do a force check

       call force_check

! open output files

       call output_setup('start')

!*****************************START MD SIMULATION***********************!

       do i = 1 , nstep

! update positions in first stage of VV integrator

            call integrator(1)

! calculate forces of new positions

            call forces

! update velocities in second stage of VV integrator

            call integrator(2)

! calculate thermodynamic properties

            if (mod(i,df_thermo) .eq. 0) call thermo_stuff()

! apply thermostat

            
            if ((trim(nvt_type) .eq. "rescale" .or. &
                trim(nvt_type) .eq. "anderson")) then

                 call thermostat(nvt_type,temp_inst)

            endif

! write traj

            if (mod(i,df_xyz) .eq. 0) call thermo_dump(i)

       enddo

! close output files

       call output_setup('end')

end program classical_md
