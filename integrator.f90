subroutine integrator(Stage,natoms, mass, timestep, xo,yo,zo,x,y,z,vxo,vyo,vzo,&
                      vx, vy,vz,Fxo,Fyo,Fzo, Fx,Fy,Fz)

      use kinds
      implicit none
  
      integer (kind = ip) :: Stage, natoms, i,j
      real (kind = dp) :: timestep
      real (kind = dp),dimension(natoms) :: mass, xo, yo, zo, x, y, z,vxo,vyo,vzo,Fxo,Fyo,Fzo
      real (kind = dp),dimension(natoms) :: Fx,Fy,Fz,vx,vy,vz
       
      if(Stage .eq. 1) then

           do i = 1, natoms
     
               x(i) = xo(i) + vxo(i)* timestep +(Fxo(i)*timestep**2)/2*mass(i)
               y(i) = yo(i) + vyo(i)* timestep +(Fyo(i)*timestep**2)/2*mass(i)
               z(i) = zo(i) + vzo(i)* timestep +(Fzo(i)*timestep**2)/2*mass(i)
               
               xo(i) = x(i)
               yo(i) = y(i)
               zo(i) = z(i)
          
           enddo

      else if(Stage .eq. 2) then 

            do j = 1, natoms
               
               vx(j) = vxo(j) + ((Fx(j) +Fxo(j))/2)*timestep*mass(j)
               vy(j) = vyo(j) + ((Fy(j) +Fyo(j))/2)*timestep*mass(j)
               vz(j) = vzo(j) + ((Fz(j) +Fzo(j))/2)*timestep*mass(j)

               vxo(j) = vx(j)
               vyo(j) = vy(j)
               vzo(j) = vz(j)
           
               Fxo(j) = Fx(j)
               Fyo(j) = Fy(j)
               Fzo(j) = Fz(j)
            
            enddo
     
       
      endif

end subroutine
