subroutine integrator(Stage,mass,dt, xo,yo,zo,vxo,vyo,vzo,fxo,fyo,fzo)
!************************************************************************
!
! Subroutine for updating positions and velocities using the velocity
! verlet algorithm
!
!************************************************************************

      use common_variables

      implicit none
  
      integer (kind = ip) :: Stage,i,j
      real (kind = dp) :: dt
      real (kind = dp),dimension(n_atoms) :: mass,xo,yo,zo,vxo,vyo,vzo,fxo,fyo,fzo
    
! first stage of VV updates positions
   
      if(Stage .eq. 1) then

           do i = 1, n_atoms
     
! updates positions to t + dt

               x(i) = xo(i) + vxo(i)* dt +(Fxo(i)*dt**2)/2*mass(i)
               y(i) = yo(i) + vyo(i)* dt +(Fyo(i)*dt**2)/2*mass(i)
               z(i) = zo(i) + vzo(i)* dt +(Fzo(i)*dt**2)/2*mass(i)
 
! wrap positions into box

               x(i) = x(i) - (Lx * anint(x(i)/Lx))
               y(i) = y(i) - (Ly * anint(y(i)/Ly))
               z(i) = z(i) - (Lz * anint(z(i)/Lz))
              
! save new positions to old positions

               xo(i) = x(i)
               yo(i) = y(i)
               zo(i) = z(i)
          
           enddo

! second stage of VV updates velocities 

      else if(Stage .eq. 2) then 

            do j = 1, n_atoms

! calculate t + dt velocities
               
               vx(j) = vxo(j) + ((fx_tot(j) +Fxo(j))/2)*dt*mass(j)
               vy(j) = vyo(j) + ((fy_tot(j) +Fyo(j))/2)*dt*mass(j)
               vz(j) = vzo(j) + ((fz_tot(j) +Fzo(j))/2)*dt*mass(j)

! save new velocities to old velocities

               vxo(j) = vx(j)
               vyo(j) = vy(j)
               vzo(j) = vz(j)

! save old forces to new forces
           
               fxo(j) = fx_tot(j)
               fyo(j) = fy_tot(j)
               fzo(j) = fz_tot(j)
            
            enddo
     
       
      endif

end subroutine
