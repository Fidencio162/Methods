	module Mathematical
	implicit double precision (a-h,o-z)
	double precision, parameter::wtrap=0.5d0
	double precision, parameter::wsimp13=1.d0/6.d0
	double precision, parameter::wsimp38=1.d0/24.d0
	contains
!----------------------------------------------
	subroutine Riemann(fint,fk0,dr)
! Riemann sum
	fint=fk0*dr
	return
	end
	subroutine trapezoidal(fint,fk0,fk1,dr)
! trapezoidal rule
	fint=wtrap*(fk0+fk1)*dr
	return
	end
! Simpson's rule 1/3
	subroutine Simpson13(fint,fk0,fk1,fk2,dr)
	fint=wtrap*wsimp13*(fk0+4.d0*fk1+fk2)*dr
	return
	end
! Simpson's rule 3/8
	subroutine Simpson38(fint,fk0,fk1,fk2,fk3,dr)
	fint=wsimp38*(fk0+3.d0*fk1+3.d0*fk2+fk3)*dr
	return
	end
	
	
	!	Spherical Bessel Functions
	function Besselj0(qm,rm)
	implicit double precision(a-h,o-z)
	Besselj0=dsin(qm*rm)/(qm*rm)
	return
	end
	!--------------------------------------------------------------
	function Besselj1(qm,rm)
	implicit double precision(a-h,o-z)
	Besselj1=dsin(qm*rm)/(qm*rm)**2-dcos(qm*rm)/(qm*rm)
	return
	end
	!--------------------------------------------------------------
	function Besselj2(qm,rm)
	implicit double precision(a-h,o-z)
	Besselj2=(3.d0/(qm*rm)**2-1.d0)*dsin(qm*rm)/(qm*rm)
     & 	-3.d0*dcos(qm*rm)/(qm*rm)**2
	return
	end
	! Determinant 2x2, 3x3 and 4x4
	function Det22(f11,f12,f21,f22)
	implicit double precision(a-h,o-z)
	Det22=f11*f22-f12*f21
	return
	end
	!--------------------------------------------------------------
	function Det33(f11,f12,f13,f21,f22,f23,f31,f32,f33)
	implicit double precision(a-h,o-z)
	Det33=f11*Det22(f22,f23,f32,f33)-f12*Det22(f21,f23,f31,f33)
     &	      +f13*Det22(f21,f22,f31,f32)
	return
	end
	!--------------------------------------------------------------
	function Det44(f11,f12,f13,f14,f21,f22,f23,f24,f31,f32,f33,f34,
     &		       f41,f42,f43,f44)
	implicit double precision(a-h,o-z)
	Det44=f11*Det33(f22,f23,f24,f32,f33,f34,f42,f43,f44)
     &	     -f12*Det33(f21,f23,f24,f31,f33,f34,f41,f43,f44)
     &	     +f13*Det33(f21,f22,f24,f31,f32,f34,f41,f42,f44)
     &	     -f14*Det33(f21,f22,f23,f31,f32,f33,f41,f42,f43)
	return
	end
	end module
