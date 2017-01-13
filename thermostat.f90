subroutine termostat(nvt_type,KE,temp_inst)

     use common_variables
     use constants
     implicit none

     integer (kind = ip) :: i,j,type,dt
     integer :: n
     real (kind = dp) :: temp,temp_inst,lamda,KE,vxb,vyb,vzb,numx,nu
     integer,dimension(:),allocatable :: seed
     character :: nvt_type
 
     call random_seed(size=n)
     allocate(seed(n))
     call date_and_time(values=dt)
     seed(1:) = dt*(/(j,j = 1,n)/)
     call random_seed(put=seed)
 
     if(nvt_type .eq. "rescale") then 

          temp_inst = 2*KE/((3*n_atoms-6)*kb)
          do i = 1,n_atoms

         
              if (temp_inst .gt. 0) then
 
                 lamda = sqrt(temp/temp_inst)
              
                 vx(i) = lamda*vx_o(i)
                 vy(i) = lamda*vy_o(i)
                 vz(i) = lamda*vz_o(i)
 
              endif

           enddo

     elseif(nvt_type .eq.  "andersen") then

           
           do j = 1,n_atoms
 
               call random_number(numx)
               if(numx .lt. nu*dt) then
                 
                 call boltz_vel(vxb,vyb,vzb)                  
                 vx(j)= vxb
                 vy(j)= vyb
                 vz(j)= vzb  
               
               endif
           enddo
     endif

end subroutine

       





 
