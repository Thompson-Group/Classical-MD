subroutine boltz_vel(vxb,vyb,vzb,mass)
!************************************************************************
!
! Subroutine for picking new atomic velocity from Boltzmann distribution
!
!************************************************************************

       use common_variables
       use constants

       implicit none

       integer(kind = ip) :: clock(8),i
       integer :: n
       integer,allocatable,dimension(:) :: seed
       real(kind=dp) :: RT, mtot, sx, sy, sz, num1, num2, r1, r2
       real(kind=dp) :: rsq, fac, gset, gauss
       real(kind=dp) :: vxb, vyb, vzb, mass

       call random_seed(size = n)
       allocate(seed(n))
       call date_and_time(values=clock)
       seed(1:) = clock(6) * (/ (i, i = 1 , n) /)
       call random_seed(put=seed)

       RT = kb * temp

! get random number to assign sign of each velocity component

       call random_number(sx)
       call random_number(sy)
       call random_number(sz)

! call uniform numbers and transform to gaussian 
 
       transform:do

            call random_number(num1)
            call random_number(num2)

            r1 = 2.0_dp * num1 - 1.0_dp
            r2 = 2.0_dp * num2 - 1.0_dp
            rsq = r1*r1 + r2*r2
 
            if((rsq .lt. 1.0_dp) .and. (rsq .ne. 0.0_dp)) exit transform

       enddo transform

       fac = sqrt(-2.0_dp * log(rsq)/rsq)
       gset = r1*fac
       gauss = r2*fac

       vxb = gauss*sign(1.0_dp,2.0_dp*sx-1.0_dp)*sqrt(RT/mass)
       vyb = gauss*sign(1.0_dp,2.0_dp*sy-1.0_dp)*sqrt(RT/mass)
       vzb = gauss*sign(1.0_dp,2.0_dp*sz-1.0_dp)*sqrt(RT/mass)

     end subroutine boltz_vel

       
 
