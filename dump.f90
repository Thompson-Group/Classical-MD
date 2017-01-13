subroutine dump(i)

use common_variables
implicit none

integer(kind=ip),intent(in) :: i
open(20,file=dumpfile)

	write(30,*) n_atoms
	write(30,*) 'Atoms. Timestep: ', i
	do j=1,n_atoms
		write(30,*) n_a_type(j), x(j), y(j), z(i)
	end do

end subroutine dump
