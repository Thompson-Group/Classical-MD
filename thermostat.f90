subroutine termostat(kb,KE,Temp,Temp_c,lamda,vxo,vyo,vzo)

     use common_variables
     implicit none

     integer (kind = ip) :: i
     real (kind = dp) :: kb,Temp,Temp_c,lamda,KE
     real (kind = dp),dimension(n_atoms) :: vxo,vyo,vzo
 
     do i = 1,n_atoms

         Temp_c = 2*KE/n_atoms*kb
        
           if (Temp_c .gt. 0) then
 
               lamda = sqrt(Temp/Temp_c)
              
               vx = lamda*vxo
               vy = lamda*vyo
               vz = lamda*vzo
 
           endif

      enddo

end subroutine

       





 
