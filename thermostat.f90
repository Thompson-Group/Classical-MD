subroutine termostat(type,dt,n,seed,numx,nu,kb,KE,Temp,Temp_c,lamda,vxo,vyo,vzo)

     use common_variables
     implicit none

     integer (kind = ip) :: i,j,type,dt(4)
     integer :: n
     real (kind = dp) :: kb,Temp,Temp_c,lamda,KE,nu,vxb,vyb,vzb,numx
     real (kind = dp),dimension(n_atoms) :: vxo,vyo,vzo
     integer,dimension(:),allocatable :: seed
 
     call random_seed(size=n)
     allocate(seed(n))
     call date_and_time(values=dt)
     seed(1:) = dt(4)*(/(j,j = 1,n)/)
     call random_seed(put=seed)
 
     if(type .eq. 1) then 

  
          do i = 1,n_atoms

              Temp_c = 2*KE/n_atoms*kb
         
              if (Temp_c .gt. 0) then
 
                 lamda = sqrt(Temp/Temp_c)
              
                 vx(i) = lamda*vxo(i)
                 vy(i) = lamda*vyo(i)
                 vz(i) = lamda*vzo(i)
 
              endif

           enddo

     elseif(type .eq.  2) then

           
           do j = 1,n_atoms
 
               call random_number(numx)
               if(numx .lt. nu*dt(4)) then
                 
                 call boltz_vel(vxb,vyb,vzb)                  
                 vx(j)= vxb
                 vy(j)= vyb
                 vz(j)= vzb  
               
               endif
           enddo
     endif

end subroutine

       





 
