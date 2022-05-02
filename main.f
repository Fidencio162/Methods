	program main
	use Mathematical
	implicit double precision(a-h,o-z)
	parameter (ni=1000000)
	dimension fp(ni),r(ni)
	pi=dacos(-1.d0)
	dr=12.d0/ni
	print*, wtrap
	do k=1,ni
	  r(k)=dr*k
	  fp(k)=dexp(-r(k))
	enddo
	! Integrando
	sriem=0.d0
	strap=0.d0
	ssimp13=0.d0
	ssimp38=0.d0
	do i=1,ni
	!dx=r(i+1)-r(i)
	fp0=fp(i)
	call Riemann(fintr,fp0,dr)
	sriem=sriem+fintr
	enddo
	do i=1,ni-1
	dx=r(i+1)-r(i)
	fp1=fp(i+1)
	fp0=fp(i)
	call trapezoidal(fintt,fp0,fp1,dx)
	strap=strap+fintt
	enddo
	!--------------------------------------------------------------
	do i=1,ni-2
	dx=r(i+2)-r(i)
	fps2=fp(i+2)
	fps1=fp(i+1)
	fps0=fp(i)
	call Simpson13(fints,fps0,fps1,fps2,dx)
	ssimp13=ssimp13+fints
	enddo
	!--------------------------------------------------------------
	do i=1,ni-3
	dx=r(i+3)-r(i)
	fpss3=fp(i+3)
	fpss2=fp(i+2)
	fpss1=fp(i+1)
	fpss0=fp(i)
	call Simpson38(fintss,fpss0,fpss1,fpss2,fpss3,dx)
	ssimp38=ssimp38+fintss
	enddo
	print*, 'La integral es (Riemman_):  ',sriem
	print*, 'La integral es (Trapacio):  ',strap
	print*, 'La integral es (Simpson1/3):',ssimp13
	print*, 'La integral es (Simpson3/8):',ssimp38
	Detp3=Det33(2.d0,3.d0,1.d0,0.d0,-2.d0,5.d0,3.d0,7.d0,-4.d0)
	Detp=Det44(1.d0,4.d0,2.d0,1.d0,-1.d0,-1.d0,3.d0,2.d0,0.d0,5.d0
     &	    ,7.d0,-4.d0,2.d0,1.d0,-3.d0,2.d0)
	print*, 'El valor del determinante es:',Detp
	end program
