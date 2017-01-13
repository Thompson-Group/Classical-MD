subroutine calc_vdw()
!************************************************************************
! 
! Subroutine for calculating LJ energy and force.
!
!************************************************************************

     use common_variables
 
     implicit none

     integer(kind = ip) :: i,j
     real(kind = dp) :: dist

! initialize
    
     fx_v = 0.0_dp
     fy_v = 0.0_dp
     fz_v = 0.0_dp
     v_v = 0.0_dp

! loop over all atoms
     
     do i = 1, n_atoms-1

! loop over all later atoms

         do j = i+1, n_atoms 

! distance
             
            dist=sqrt(rx(i,j)**2+ry(i,j)**2+rz(i,j)**2)
     
            if (dist .lt. r_cut) then

! energy
          
               v_v = v_v + 4.0_dp*ee(i,j)*((ss(i,j)/dist**12)-(ss(i,j)/dist**6))

! force (-du/dx)

               fx_v(i) = fx_v(i) + ((-48*x(i)/dist**2)*(0.5*(ss(i,j)/dist**6)-(ss(i,j)/dist**12)))
               fy_v(i) = fy_v(i) + ((-48*y(i)/dist**2)*(0.5*(ss(i,j)/dist**6)-(ss(i,j)/dist**12)))
               fz_v(i) = fz_v(i) + ((-48*z(i)/dist**2)*(0.5*(ss(i,j)/dist**6)-(ss(i,j)/dist**12)))

! equal and opposite force

               fx_v(j) = fx_v(j) - fx_v(i)
               fy_v(j) = fy_v(j) - fy_v(i)
               fz_v(j) = fz_v(j) - fz_v(i)

            endif
         
          enddo

    enddo

end subroutine 
