subroutine calc_vdw(natoms,rx,ry,rz,x,y,z,ee,ss,r_cut,Fx_v,Fy_v,Fz_v,E_v)

     use kinds
 
     implicit none

     integer(kind = ip) :: natoms,i,j
     real(kind = dp),dimension(natoms,natoms) :: rx,ry,rz,ee,ss
     real(kind = dp) :: E_v,r_cut,dist
     real (kind = dp),dimension(natoms) :: x,y,z,Fx_v,Fy_v,Fz_v
    
     Fx_v = 0.0_dp
     Fy_v = 0.0_dp
     Fz_v = 0.0_dp
     E_v = 0.0_dp
     
     do i = 1, natoms-1

         do j = i+1, natoms 
             
            dist=sqrt(rx(i,j)**2+ry(i,j)**2+rz(i,j)**2)
     
            if (dist .lt. r_cut) then
          
               E_v = E_v + 4.0_dp*ee(i,j)*((ss(i,j)/dist**12)-(ss(i,j)/dist**6))

               Fx_v(i) = Fx_v(i) + ((-48*x(i)/dist**2)*(0.5*(ss(i,j)/dist**6)-(ss(i,j)/dist**12)))
               Fy_v(i) = Fy_v(i) + ((-48*y(i)/dist**2)*(0.5*(ss(i,j)/dist**6)-(ss(i,j)/dist**12)))
               Fz_v(i) = Fz_v(i) + ((-48*z(i)/dist**2)*(0.5*(ss(i,j)/dist**6)-(ss(i,j)/dist**12)))

               Fx_v(j) = Fx_v(j) - Fx_v(i)
               Fy_v(j) = Fy_v(j) - Fy_v(i)
               Fz_v(j) = Fz_v(j) - Fz_v(i)

            endif
         
          enddo

    enddo

end subroutine 
