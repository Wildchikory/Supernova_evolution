!! SN PROJECT

MODULE DATA
integer :: N, i
real*8 :: pi,cmpc,cmkpc,yr,kbol,mu,mp, e0
parameter (N=500) 
parameter(pi=3.141592)
parameter(cmpc=3.085d18)
parameter(cmkpc=1000.*cmpc)
parameter(yr=3.156d7)
parameter(kbol=1.38d-16)
parameter(mu=0.61)
parameter(mp=1.67d-24)
parameter(e0=1.d51)
parameter(yrsec=3.15d7)
parameter(KKeV=1.16d7)
END MODULE DATA

PROGRAM ZEUS
USE DATA
IMPLICIT NONE
real*8 :: xa(N), xb(N), xmax, xmin, deltax, dxa(N), dxb(N), time(N), rsed(N)
real*8 :: d(N), e(N), v(N), P(N), s(N), Temp(N) !DENSITA', ENERGIAINTERNA, VELOCITA', PRESSIONE, MOMENTO
real*8 :: q(N) !VISCOSITA' ARTIFICIALE
real*8 :: g2a(N), g2b(N), g31a(N), g31b(N), dvl1a(N), dvl1b(N) 
real*8 :: F1(N), F2(N), F3(N), M(N),  dstar(N),  e_dstar(N), vstar(N), e_d(N)
real*8 :: divV(N), dt0, a, b, lambda(N),lambda2(N*100), lambda_wo(N), dtemp
real*8 :: Ecin, Eter,  EterIN, lumx, tempt(N*100), Etot, sed(N)
real*8 :: dtmin, tmax, t, C2, gam, cv, k, t1, t2, t3, cfl
real*8 :: SedLaw !Costante per la legge di Sedov
real*8, dimension(9) :: rho_s, r_s, tempo, r_s2, tempo2
integer :: sdr, Num, stamp, ncicli, t0, o, step
real*8, EXTERNAL :: Cool
!temp2(N*100),

!CREAZIONE DOPPIA GRIGLIA (xa e xb)

xmin=0.
xmax=70.*cmpc

!GRIGLIA "a"
do i=1,N
	xa(i)= xmin+(xmax-xmin)*(i-1.)/(N-1.)
end do

deltax=xa(3)-xa(2)

! set grid
xa(1) = -deltax
do i=2,N
	xa(i)=xa(i-1)+deltax
end do
print*, 'xa1=', xa(1)/cmpc
print*, 'xaN=', xa(N)/cmpc
print*, 'deltaxa=', deltax/cmpc
!GRIGLIA "b"
do i=1, N-1
	xb(i)=0.5*(xa(i)+xa(i+1))
end do
xb(N)=xb(N-1)+(xb(N-1)-xb(N-2))   !! add the last calculated Delta_xb to xb(N-1)
print*, 'xb1=', xb(1)/cmpc
print*, 'xbN=', xb(N)/cmpc
print*, 'deltaxb=', deltax/cmpc
do i=2, N-1
	dxa(i)=xa(i+1)-xa(i)
	dxb(i)=xb(i)-xb(i-1)
end do

dxa(1)=xa(2)-xa(1)
dxa(N)=dxa(N-1)
dxb(1)=dxb(2)
dxb(N)=xb(N)-xb(N-1)

open(20,file='grid2.dat')
do i=1,N
   write(20,1001)xa(i),xb(i),dxa(i),dxb(i)
enddo
close(20)
1001 format(4(1pe12.4))

!DEFINIZIONE FATTORI DI SCALA METRICI 

do i=1, N
	g2a(i)=xa(i)
	g31a(i)=xa(i)
	g2b(i)=xb(i)
	g31b(i)=xb(i)
end do

do i=1, N-1
	dvl1a(i)=(xa(i+1)**3-xa(i)**3)/3.
end do
dvl1a(N)=dvl1a(N-1)
do i=2, N
	dvl1b(i)=(xb(i)**3-xb(i-1)**3)/3.
end do
	dvl1b(1)=dvl1b(2)




!IMPLEMENTAZIONE CONDIZIONI INIZIALI

 gam=5./3.
 cv=1.99d8    !! warning: this is right for gam = 5/3 !!
 t=0.
 c2=3.
 cfl=0.01

do i=1, N
  d(i)=2.d-24
  TEMP(i)=1.d4
  e(i)=cv*d(i)*TEMP(i)
  P(i)=(gam-1.)*e(i)       
  v(i)=0.
end do	

EterIN=2.d-24*cv*1.d4*4./3.*pi*(70.*cmpc)**3  !initial energy

!SN INJECTION

e(2)= e0/(4./3.*pi*xa(4)**3)
e(3)=e(2)
TEMP(2)=e(2)/(cv*d(2))
TEMP(3)=e(3)/(cv*d(3))
P(2)=(gam-1)*e(2)
P(3)=(gam-1)*e(3)

CALL Bcb(d)
CALL Bcb(p)
CALL Bcb(temp)
CALL Bcb(e)


ncicli=0
tmax =1.d4*yrsec
step=1

print*,'How do you want to calculate profiles? 1. Without cooling, 2. with cooling'
read(*,*) o !variable that tells apart the different possibilities for profiles

if (o==1) then
	open(30, file='lum_wo.dat')
	print*, 'wo'
else if (o==2) then
	open(30, file='luminosity.dat')
	print*, 'cool'
end if
1010 continue

do while (t<tmax)      !!!! HERE STARTS THE TIME INTEGRATION !!!!!
        ncicli=ncicli+1       
        do i=1,N
			P(i)=(gam-1.)*e(i)
        end do

!CALCOLO DTMIN
        dtmin=1.d30   !! any very large value !!
	do i=2, N-1
		 dtmin=min(dtmin,(xb(i)-xb(i-1))/(abs(v(i))+sqrt(gam*p(i)/d(i))))
	end do
	
	if (ncicli==0) then
		cfl=0.01
	 else if ((ncicli>0).and.(cfl<0.5)) then
		cfl=cfl+0.1*cfl
	 else if (cfl>=0.5) then
		cfl=0.5
	end if
        dtmin=cfl*dtmin
        t=t+dtmin
        !print*,'ncicli, dtmin = ',ncicli, real(dtmin),real(t)

	if (o==2) then
		do i=1, N  
			lambda(i)=cool(temp(i)) !cooling
			e(i)= e(i)-lambda(i)*(d(i)/(2.17*1.d-24))**2*dtmin !cooling
			P(i)=(gam-1.)*e(i)
			temp(i)=e(i)/(cv*d(i))!cooling
			if (temp(i)<1.d4) then
				temp(i)=1.d4  !cooling
				e(i)=cv*d(i)*temp(i)  !cooling
				P(i)=(gam-1.)*e(i)  !cooling
			end if
		end do
	end if
	
CALL Bcb(d)
CALL Bcb(p)
CALL Bcb(temp)
CALL Bcb(e)

!SOURCE STEP
!SUBSTEP I: AGGIORNAMENTO DELLA VELOCITÀ PER GRADIENTE DI P

	do i=2, N-1
		v(i)=v(i)-dtmin*2.*(P(i)-P(i-1))/((d(i)+d(i-1))*dxb(i))	
	end do
	CALL BCa(v)


!CALCOLO Q
	do i=2, N-1
		if ((v(i+1)-v(i))<0.) then
			q(i)=C2*d(i)*(v(i+1)-v(i))**2
		else 
			q(i)=0.
		end if
	end do
	CALL BCb(q)

!SUBSTEP II: AGGIORNAMENTO PER VISCOSITÀ ARTIFICIALE

	do i=2, N-1
		v(i)=v(i)-dtmin*2.*(q(i)-q(i-1))/((d(i)+d(i-1))*dxb(i))
	end do
	CALL BCa(v)

	do i=2, N-1
		e(i)=e(i)-dtmin*q(i)*(v(i+1)-v(i))/dxa(i)
	end do
	CALL BCb(e)

!SUBSTEP III: AGGIORNAMENTO PER RISCALDAMENTO DA COMPRESSIONE
	do i=2,N-1
		divV(i)=(g2a(i+1)*g31a(i+1)*v(i+1)-g2a(i)*g31a(i)*v(i))/dvl1a(i)
	end do
	CALL BCa(divV)

	do i=2, N-1
		e(i)=e(i)*(1.-0.5*dtmin*(gam-1.)*divV(i))/(1.+0.5*dtmin*(gam-1.)*divV(i))
	end do
	CALL BCb(e)

!!  Here update T when needed (not needed for the shock tube)



!!!!!!TRANSPORT STEP (use Upwind first order only)

	do i=2, N-1       !! here define the momentum density
		s(i)=0.5*(d(i)+d(i-1))*v(i)  !! this is at "i" !!
	end do	

	CALL BCa(s)

!AGGIORNAMENTO DENSITÀ

	do i=2, N-1       !! here select the value of the density at the interface "i"
		if (v(i)>0.) then
			dstar(i)=d(i-1)     !! at i !!
		else
			dstar(i)=d(i)
		end if
	end do
	dstar(N)=dstar(N-1)
	dstar(1)=dstar(3)

	do i=2, N
		F1(i)=dstar(i)*v(i)*g2a(i)*g31a(i)    !! at i !!	
	end do

!AGGIORNAMENTO ENERGIA

	do i=2, N-1
		M(i)=dstar(i)*v(i)
	end do
	CALL BCa(M)
	
	
	do i=2, N-1
		if (v(i)>0.) then
			e_dstar(i)=e(i-1)/d(i-1)   !! at i !!
		else
			e_dstar(i)=e(i)/d(i)
		end if
	end do
	e_dstar(N)=e_dstar(N-1)
	e_dstar(1)=e_dstar(3)

	!ORA AGGIORNO LA DENSITÀ
	do i=2, N-1
		d(i)=d(i)-dtmin*(F1(i+1)-F1(i))/dvl1a(i)
	end do 
	CALL BCb(d)
	

	do i=2, N
		F2(i)=e_dstar(i)*M(i)*g2a(i)*g31a(i)				
	end do
	CALL BCa(F2)

	do i=2, N-1
		e(i)=e(i)-dtmin*(F2(i+1)-F2(i))/dvl1a(i)
		TEMP(i)=e(i)/(cv*d(i)) !cambio temperatura
	end do

	CALL BCb(e)


!AGGIORNAMENTO MOMENTO 

	do i=2, N-1
		if ((v(i-1)+v(i))*0.5>0) then
			vstar(i)=v(i-1)       !! at i-1/2  !!
		else
			vstar(i)=v(i)
		end if
	end do

	CALL BCb (vstar)

	do i=1, N-1
		F3(i)=vstar(i+1)*0.5*(M(i)+M(i+1))*g2b(i)*g31b(i)   !! questo e' a i+1/2, occhio !!  
	end do
	
	do i=2, N-1
		s(i)=s(i)-dtmin/dvl1b(i)*(F3(i)-F3(i-1))
	end do

	CALL BCa(s)

	do i=2, N-1
		v(i)=2.*s(i)/(d(i)+d(i-1))
	end do
	!CALL BCa(v) !aggiunto
	
	if (o==1) then
		do i=1,N
			lambda_wo(i)=cool(temp(i))
		end do
		CALL luminosity(d, lambda_wo, temp, xa, dxa, lumx)
	else if (o==2) then
		CALL luminosity(d, lambda, temp, xa, dxa, lumx)
	end if
	
	CALL energy(d, v, e, xa, Ecin, Eter)
	Eter=(Eter-EterIN)/1.d51
	Ecin=Ecin/1.d51
	Etot=Eter+Ecin
    write(30, 1004) t/yrsec, lumx, Ecin, Eter, Etot
    
enddo       !! here the "do while" ends !!

if (tmax==1.d4*yrsec) then
    PRINT*, '1'
    if (o==1) then
		open(20,file='results1wo.dat')
	else
		open(20,file='results1.dat')
	end if

	do i=3,N  
		write (20,1000) xa(i)/cmpc,xb(i)/cmpc,d(i),v(i),temp(i),p(i)
	end do
	close(20)
	CALL maximum(d, xa, rho_s(1), r_s(1))
	tempo(1)=tmax
	tmax=2.d4*yrsec
	step=step+1
	go to 1010
else if (tmax==2.d4*yrsec) then
    PRINT*, '2'
    if (o==1) then
		open(20,file='results2wo.dat')
	else
		open(20,file='results2.dat')
	end if

	do i=3,N
		write (20,1000) xa(i)/cmpc,xb(i)/cmpc,d(i),v(i),temp(i),p(i)
	end do
	close(20)
	CALL maximum(d, xa, rho_s(2), r_s(2))
	tempo(2)=tmax
	tmax=6.d4*yrsec
	step=step+1
	go to 1010
else if (tmax==6.d4*yrsec) then
    PRINT*, '3'
    if (o==1) then
		open(20,file='results3wo.dat')
	else
		open(20,file='results3.dat')
	end if

	do i=3,N 
		write (20,1000) xa(i)/cmpc,xb(i)/cmpc,d(i),v(i),temp(i),p(i)
	end do
	close(20)
	CALL maximum(d, xa, rho_s(3), r_s(3))
	tempo(3)=tmax
	tmax=8.d4*yrsec
	step=step+1
	go to 1010
else if (tmax==8.d4*yrsec) then
    PRINT*, '4'
    if (o==1) then
		open(20,file='results4wo.dat')
	else
		open(20,file='results4.dat')
	end if

	do i=3,N
		write (20,1000) xa(i)/cmpc,xb(i)/cmpc,d(i),v(i),temp(i),p(i)
	end do
	close(20)
	CALL maximum(d, xa, rho_s(4), r_s(4))
	tempo(4)=tmax
	tmax=1.d5*yrsec
	step=step+1
	go to 1010
else if (tmax==1.d5*yrsec) then
    PRINT*, '5'
    if (o==1) then
		open(20,file='results5wo.dat')
	else
		open(20,file='results5.dat')
	end if

	do i=3,N  
		write (20,1000) xa(i)/cmpc,xb(i)/cmpc,d(i),v(i),temp(i),p(i)
	end do
	close(20)
	CALL maximum(d, xa, rho_s(5), r_s(5))
	tempo(5)=tmax
	tmax=2.d5*yrsec
	step=step+1
	go to 1010
else if (tmax==2.d5*yrsec) then
    PRINT*, '6'
    if (o==1) then
		open(20,file='results6wo.dat')
	else
		open(20,file='results6.dat')
	end if

	do i=3,N 
		write (20,1000) xa(i)/cmpc,xb(i)/cmpc,d(i),v(i),temp(i),p(i)
	end do
	close(20)
	CALL maximum(d, xa, rho_s(6), r_s(6))
	tempo(6)=tmax
	tmax=3.d5*yrsec
	step=step+1
	go to 1010
else if (tmax==3.d5*yrsec) then
    PRINT*, '7'
    if (o==1) then
		open(20,file='results7wo.dat')
	else
		open(20,file='results7.dat')
	end if

	do i=3,N 
		write (20,1000) xa(i)/cmpc,xb(i)/cmpc,d(i),v(i),temp(i),p(i)
	end do
	close(20)
	CALL maximum(d, xa, rho_s(7), r_s(7))
	tempo(7)=tmax
	tmax=4.d5*yrsec
	step=step+1
	go to 1010
else if (tmax==4.d5*yrsec) then
    PRINT*, '8'
    if (o==1) then
		open(20,file='results8wo.dat')
	else
		open(20,file='results8.dat')
	end if

	do i=3,N 
		write (20,1000) xa(i)/cmpc,xb(i)/cmpc,d(i),v(i),temp(i),p(i)
	end do
	close(20)
	CALL maximum(d, xa, rho_s(8), r_s(8))
	tempo(8)=tmax
	tmax=5.d5*yrsec
	step=step+1
	go to 1010
else if (tmax==5.d5*yrsec) then
    PRINT*, '9'
    if (o==1) then
		open(20,file='results9wo.dat')
	else
		open(20,file='results9.dat')
	end if

	do i=3,N 
		write (20,1000) xa(i)/cmpc,xb(i)/cmpc,d(i),v(i),temp(i),p(i)
	end do
	close(20)
	CALL maximum(d, xa, rho_s(9), r_s(9))
	tempo(9)=tmax
end if
1004 format(5(1pe12.4))
close(30)

1000 format(6(1pe12.4))

if (o==1) then
	open(99, file='Rshock.dat')
else
	open(99, file='Rshock_cool.dat')
end if

do i=1,9
	sed(i)=(2.*e0/rho_s(i))**(1./5.)*tempo(i)**(2./5.)
	write(99,1002) r_s(i)/cmpc, rho_s(i), tempo(i)/yrsec, sed(i)/cmpc
	r_s2(i)=log10(r_s(i)/cmpc)
	tempo2(i)=log10(tempo(i)/yrsec)
end do
close(99)

time(1)=1.d4*yrsec   
dt0=(5.d5-1.d4)*yrsec/(N-1)
rsed(1)=(2.*e0/(2.d-24))**(1./5.)*time(1)**(2./5.)
do i=2,N
	time(i)=time(i-1)+dt0
	rsed(i)= (2.*e0/2.d-24)**(1./5.)*time(i)**(2./5.)
end do

1002 format(4(1pe12.4))
open(99, file='Sedov.dat')

do i=1,N
	write(99, 1003) time(i)/yrsec, rsed(i)/cmpc
end do
close(99)


CALL fit_lineare(tempo2, r_s2, 9, a, b)
print*, a
print*, b

!cooling function expression

open(99, file='cooling.dat')
do i=1,100*N
	tempt(i)=1.d4+(1.d8-1.d4)*(i-1.)/(N*100-1.)
	lambda2(i)=Cool(tempt(i))
	write(99, 1003) tempt(i), lambda2(i)
end do
close(99)

1003 format(2(1pe12.4))

END PROGRAM ZEUS


SUBROUTINE BCa(z1) !corrette BC per velocità e momento (riflessione)
USE DATA
IMPLICIT NONE
real*8, dimension (N) :: z1

z1(2)=0.
z1(1)=-z1(3)
z1(N)=z1(N-1)
z1(1)=z1(2)       !! ouflow !!
z1(N)=z1(N-1)

END SUBROUTINE BCa

SUBROUTINE BCb(z2) ! BC di outflow tradizionali
USE DATA
IMPLICIT NONE
real*8, dimension (N) :: z2
z2(1)=z2(2)
z2(N)=z2(N-1)
END SUBROUTINE BCb

SUBROUTINE maximum(z3, z4, dm, rm)
USE DATA
IMPLICIT NONE
REAL*8, dimension(N) :: z3, z4
REAL*8 :: maxi, rmax, dm, rm

i=1
maxi=z3(1)
rmax=z4(1)
DO WHILE (i<=N)
	if (z3(i)>maxi) then
		maxi=z3(i)
		rmax=z4(i)
	end if
	i=i+1
END DO
 dm=maxi
 rm=rmax
END SUBROUTINE maximum

SUBROUTINE fit_lineare(x, y, nd, c1, c2)
	USE DATA
    IMPLICIT NONE
    INTEGER::nd
    REAL*8::x(nd), y(nd), c1, c2, delta, q(nd), xy(nd)

    DO i=1,nd
        q(i)=x(i)**2
        xy(i)=x(i)*y(i)
    END DO

    delta=nd*SUM(q)-(SUM(x))**2
    c2=(SUM(q)*SUM(y)-SUM(x)*SUM(xy))/delta
    c1=(nd*SUM(xy)-SUM(x)*SUM(y))/delta

END SUBROUTINE

SUBROUTINE luminosity(rho, lamb, tmp, r, dr, lum)
	USE DATA
    IMPLICIT NONE
    REAL*8::rho(N), lamb(N), tmp(N), r(N), dr(N), lum, l
    l=0
	do i=1,N
		if (tmp(i)>1.d6) then
			l=l+lamb(i)*4*pi*r(i)**2*dr(i)*(rho(i)/2.17d-24)**2
		end if
	end do
	lum=l
END SUBROUTINE

SUBROUTINE energy(rho, vel, en,r, k, th)
	USE DATA
	IMPLICIT NONE
	
	REAL*8::rho(N),vel(N), en(N), r(N), k, th
	k=0.
	th=0.
	do i=2,N
		k=k+vel(i)**2*0.5*4./3.*pi*(r(i)**3-r(i-1)**3)*(rho(i))
		th=th+en(i)*4./3.*pi*(r(i)**3-r(i-1)**3)
	end do
END SUBROUTINE

Real*8 FUNCTION Cool(Temp1)
USE DATA
IMPLICIT NONE
Real*8:: Temp1
		if ((Temp1/kkev)>0.02) then
			cool=1.d-22*(8.6*1.d-3*(temp1/kkev)**(-1.7)+0.058*(temp1/kkev)**0.5+0.063)
		else if (((Temp1/kkev)<=0.02).and.((Temp1/kkev)>=0.0017235)) then
			cool=6.72*1.d-22*(temp1/kkev/0.02)**0.6			
		else if ((Temp1/kkev)<=0.0017235) then
			cool=1.544*1.d-22*(temp1/kkev/0.0017235)**6
		end if
 !erg s-1 cm-3
END FUNCTION Cool
