subroutine dists(x,y,z,l,natoms,rx,ry,rz)
implicit none

integer :: natoms
real(kind=dp), dimension(natoms,natoms) :: rx, ry, rz, x, y, z
real(kind=dp), dimension(natoms) :: xtmp, ytmp, ztmp
real(kind=dp), dimension(3) :: l
! working variables
integer :: j, i


do j = 1, natoms-1
   do i = j+1, natoms
      xtmp(j) = x(j)-anint((x(j)-x(i))/l(1))*l(1)
      ytmp(j) = y(j)-anint((y(j)-y(i))/l(2))*l(2)
      ztmp(j) = z(j)-anint((z(j)-z(i))/l(3))*l(3)
      rx(i,j) = x(i)-xtmp(j)
      ry(i,j) = y(i)-ytmp(j)      
      rz(i,j) = z(i)-ztmp(j)
   enddo
enddo
rx(j,i) = rx(i,j)
ry(j,i) = ry(i,j)
rz(j,i) = rz(i,j)

end subroutine dists
