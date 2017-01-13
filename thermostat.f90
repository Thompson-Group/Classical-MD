subroutine termostat(nvt_type)

     use common_variables
     use constants
     implicit none

     integer (kind = ip) :: i,j,step(8)
     integer :: n
     real (kind = dp) :: temp_inst,KE,vxb,vyb,vzb,numx,lamda,nu
     integer,dimension(:),allocatable :: seed
     character :: nvt_type
 
     call random_seed(size=n)
     allocate(seed(n))
     call date_and_time(values=step)
     seed(1:) = step*(/(j,j = 1,n)/)
     call random_seed(put=seed)

     
 
     if(nvt_type .eq. "rescale") then 
          do i = 1, n_atoms

            ke = ke + 0.5_dp*M(i)*(vx(i)**2 + vy(i)**2 + vz(i)**2)
          
          enddo
           
          temp_inst = 2*ke/((3*n_atoms-6)*kb)
         
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

       





 
