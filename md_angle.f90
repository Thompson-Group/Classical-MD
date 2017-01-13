module md_angle
! to calculate the angle bend energy for water

contains
	subroutine calc_angle(rx,ry,rz,nangles,natms,k,theta0,natyp,atp,At1,at2, at3,bond_energy, Fx, Fy, Fz)
	use kinds
	implicit none

	integer(kind=ip),intent(in) :: natms, nangles, natyp
	real(kind=dp),intent(in)    :: rx(natms,natms),ry(natms,natms),rz(natms,natms)
        real(kind=dp),intent(in)    :: k(natyp),theta0(natyp)
        real(kind=dp),intent(out)    ::bond_energy
        integer(kind=ip), intent(in) :: at1(nangles), at2(nangles), at3(nangles),atp(nangles)
	real(kind=dp),intent(out)   :: fx(natms),fy(natms),fz(natms)
	integer(kind=ip) :: i,a,b,c,d
	real(kind=dp) :: dot_prod, catmp,cfac,lenR1, lenR2, rad_ang

        fx = 0; fy = 0; fz= 0
	bond_energy = 0.0
	do i=1,nangles
		a = AT1(i)
		b = AT2(i)
		c = AT3(i)
                d = Atp(i)
		
		dot_prod = rx(a,b)*rx(b,c) + ry(a,b)*ry(b,c) + rz(a,b)*rz(b,c)
		lenR1 = sqrt(rx(a,b)**2 + ry(a,b)**2 + rz(a,b)**2)
		lenR2 = sqrt(rx(b,c)**2 + ry(b,c)**2 + rz(b,c)**2)
                
                catmp = dot_prod/(lenR1*lenR2)
		rad_ang = acos(dot_prod/(lenR1*lenR2))
		
                cfac = -k(d)*(rad_ang - theta0(d))/sqrt(1 - catmp **2)
		bond_energy = bond_energy + 0.5*k(d)*(rad_ang-theta0(d))**2 
                fx(a) = fx(a) - cfac*(rx(b,c)/(lenR1*lenR2) - catmp*rx(a,b)/lenR1**2)
                fy(a) = fy(a) - cfac*(ry(b,c)/(lenR1*lenR2) - catmp*ry(a,b)/lenR1**2)
                fz(a) = fz(a) - cfac*(rz(b,c)/(lenR1*lenR2) - catmp*rz(a,b)/lenR1**2)

		fx(b) = fx(b) + cfac*(rx(a,b)+rx(b,c))/(lenR1*lenR2) - catmp*(rx(a,b)/lenR1**2 + rx(b,c)/lenR2**2)
		fy(b) = fy(b) + cfac*(ry(a,b)+ry(b,c))/(lenR1*lenR2) - catmp*(ry(a,b)/lenR1**2 + ry(b,c)/lenR2**2)
		fz(b) = fz(b) + cfac*(rz(a,b)+rz(b,c))/(lenR1*lenR2) - catmp*(rz(a,b)/lenR1**2 + rz(b,c)/lenR2**2)
               
		fx(c) = fx(c) - cfac*(rx(a,b)/(lenR1*lenR2) - catmp*rx(b,c)/lenR2**2)
		fy(c) = fy(c) - cfac*(ry(a,b)/(lenR1*lenR2) - catmp*ry(b,c)/lenR2**2)
		fz(c) = fz(c) - cfac*(rz(a,b)/(lenR1*lenR2) - catmp*rz(b,c)/lenR2**2)

         enddo
   end subroutine calc_angle
 end module   md_angle

         
