c------------------------------------------------------------------------------
c
c	SUBROUTINE NFDIF
c
c	PURPOSE
c	   Differentiation of a periodacally tabulated function.
c
c	USAGE
c	   call nfdif(fx,dx,n,fc,nfp,dfx,ier)
c
c	DESCRIPTION OF PARAMETERS
c	   fx ... Vector of tabulated function values of lenght n.
c	   dx ... Spacing of x.
c	   n .... Size of fx.
c	   fc ... Critical frequency for LF filtration
c	          suppress frequency of fc and higher
c		  omega = 2 * Pi * f
c	          fc can be from 0.0 to 0.5;
c	          one step from fx(i) to fx(i+1) with frequency f=0.5
c	          changes the period in amount of Pi.
c	   nfp .. LF filter lenght. Defines number of Fourier series
c	          coefficients for aproximation of tansfer function
c	          such that 2*nfp+1 coeficients are taken.
c	   dfx .. Resultant vector of df(x)/dx
c	          dfx(i)=0 i=1,...,nfp
c			   i=n-nfp+1,...,n
c	   ier .. Resultant error code where
c	          ier=0 no error
c	          ier=1 nfp is greater then 50
c	          ier=2 2*nfp+1 is greater then n
c
c	REMARKS
c	   * n must be geater or equal then nfp
c	   * nfp must be smaller then 50
c
c	SUBROUTINES AND FUNCTION SUBPROGRAM REQUIRED
c	   ortpol
c
c	METHOD
c	   uses digital filter described in R.W.HAMMING
c	   'DIGITAL FILTERS',PRENTICE-HALL,INC.
c	   ENGLEWOOD CLIFFS,NEW JERSEY,CHAPTER 6.4,PAGE 109.
c
c------------------------------------------------------------------------------
	subroutine NFDIF(fx,dx,n,fc,nfp,dfx,ier)
	implicit real*8 (a-h,o-z)

	parameter (PI=3.141592654)
	real*8 filter(0:100),fx(*),dfx(*)
	real*8 xh(11),yh(11),ysmot(11),coef(0:3)

	if(nfp.gt.50)then
	  ier=1
	  return
	else if(2*nfp+1.gt.n)then
	  ier=2
	  return
	endif

c filter computation
	omega=2.0*PI*fc
	f3=float(nfp+1)
	do 1 i=1,nfp
	  f1=float(i)
	  f2=f1*omega
	  smarg=PI*f1/f3
	  filter(nfp+i)=(sin(f2)/(f1*f1)-omega*cos(f2)/f1)/PI
     +			 *sin(smarg)/smarg
	  filter(nfp-i)=-filter(nfp+i)
1	continue
	filter(nfp)=0.0

c differentiation
	do 3 i=1,n-2*nfp
	  sdfx=0.0
	  do 4 j=0,2*nfp
	    sdfx=sdfx+filter(j)*fx(i+j)
4	  continue
	  dfx(nfp+i)=sdfx/dx
3	continue

c extrapolation
	nort=5
	ig=2
	do 20 i=1,nort
	  xh(i)=nfp+i
	  yh(i)=dfx(nfp+i)
20	continue
	x0=xh(nort/2+1)
	call ortpol(nort,ig,x0,xh,yh,ysmot,coef)
	do 21 i=1,nfp
	  deltx=i-(nfp+nort/2+1)
	  dfx(i)=coef(0)+coef(1)*deltx+coef(2)*deltx*deltx
21	continue

	do 22 i=1,nort
	  xh(i)=i
	  yh(i)=dfx(n-nfp-nort+i)
22	continue
	x0=xh(nort/2)
	call ortpol(nort,ig,x0,xh,yh,ysmot,coef)
	do 23 i=1,nfp
	  deltx=i+(nort/2+1)
	  dfx(n-nfp+i)=coef(0)+coef(1)*deltx+coef(2)*deltx*deltx
23	continue

	ier=0
	return

	end
