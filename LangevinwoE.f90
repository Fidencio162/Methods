program LangevinooE
	implicit double precision(a-h,o-z)
	parameter(maxt=100000)
	dimension x(maxt),y(maxt),t(maxt),vx(maxt),vy(maxt),et(maxt)
	iseed=-6539729
	open(unit=1,file='phasespaceq1.dat')
	!initial condition
	x1=0.d0
	y1=0.d0
	v1=0.d0
	v2=0.d0
	e1=0.d0
	!steps time
	dt=0.001d0
	!parameters
	S=0.01d0
	a=2.d0
	q=1.d0
	gam0=0.2d0
	d2=1.d0
	c=0.1d0
	do i=1,maxt
	t(i)=float(i)*dt
	x(i)=x1+v1*dt
	y(i)=y1+v2*dt
	vx(i)=(1.d0-gam0*dt+d2*e1*dt)*v1-a*x1*dt+dsqrt(2.d0*S*dt)*gasdev(iseed)
	vy(i)=(1.d0-gam0*dt+d2*e1*dt)*v2-a*y1*dt+dsqrt(2.d0*S*dt)*gasdev(iseed)
	v2=v1**2+v2**2
	et(i)=(1.d0-c*dt-d2*v2*dt)*e1+q*dt
	x1=x(i)
	y1=y(i)
	v1=vx(i)
	v2=vy(i)
	e1=et(i)
	write(1,*) t(i),x1,y1
	enddo
	
end program

! Random generator algorithm
      	FUNCTION ranf(Idum)
	implicit double precision(a-h,o-z)

      	PARAMETER (IM1=2147483563, IM2=2147483399, &
       AM=1./IM1, IMM1=IM1-1,IA1=40014, IA2=40692, &
     	IQ1=53668, IQ2=52774, IR1=12211,IR2=3791, &
     	NTAb=32, NDIv=1+IMM1/NTAb, EPS=1.2E-7, RNMx=1.-EPS)

	dimension iv(ntab)
      	SAVE iv, iy, idum2
      	DATA idum2/123456789/, iv/NTAb*0/, iy/0/
      	IF (Idum.LE.0) THEN
         Idum = MAX(-Idum, 1)
         idum2 = Idum
         DO j = NTAb + 8, 1, -1
            k = Idum/IQ1
            Idum = IA1*(Idum-k*IQ1) - k*IR1
            IF (Idum.LT.0) Idum = Idum + IM1
            IF (j.LE.NTAb) iv(j) = Idum
         END DO
         iy = iv(1)
      	END IF
     	k = Idum/IQ1
      	Idum = IA1*(Idum-k*IQ1) - k*IR1
      	IF (Idum.LT.0) Idum = Idum + IM1
      	k = idum2/IQ2
      	idum2 = IA2*(idum2-k*IQ2) - k*IR2
      	IF (idum2.LT.0) idum2 = idum2 + IM2
      	j = 1 + iy/NDIv
      	iy = iv(j) - idum2
      	iv(j) = Idum
      	IF (iy.LT.1) iy = iy + IMM1
      	ranf = MIN(AM*iy, RNMx)
      	RETURN
      	END
      
      	FUNCTION gasdev(idum)
	implicit double precision(a-h,o-z)
	SAVE iset,gset
	DATA iset/0/
	  if (idum.lt.0) iset=0
	  if (iset.eq.0) then
   1    v1=2.d0*ranf(idum)-1.d0
	  v2=2.d0*ranf(idum)-1.d0
	  rsq=v1**2+v2**2
	  if((rsq.ge.(1.d0)).or.(rsq.eq.(0.d0)))goto 1 
	  fac=dsqrt(-2.d0*dlog(rsq)/rsq)
	  gset=v1*fac
	  gasdev=v2*fac
	  iset=1
	  else
	  gasdev=gset
	  iset=0
	endif
	return
	END      
