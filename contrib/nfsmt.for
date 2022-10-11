c
c NOTE: MS FORTRAN Ver.4.0 BUG
c	must be compiled with /Od
c------------------------------------------------------------------------------
c
c	SUBROUTINE NFSMT
c
c	PURPOSE
c	   Smoothing of a periodacally tabulated function.
c
c	USAGE
c	   call nfsmt(fx,n,fc,nfp,ier)
c
c	DESCRIPTION OF PARAMETERS
c	   fx ... Vector of tabulated function values of lenght n.
c	   n  ... Size of fx.
c	   fc ... critical frequency for nf filtration
c	          suppress frequency of fc and higher
c		  omega = 2 * Pi * f
c	          fc can be from 0.0 to 0.5;
c	          one step from fx(i) to fx(i+1) with frequency f=0.5
c	          changes the period in amount of Pi.
c	   nfp .. LF filter lenght.defines number of Fourier series
c	          coefficients for aproximation of tansfer function
c	          such that 2*nfp+1 coeficients are taken.
c	   ier .. Resultant error code where
c	          ier=0 no error
c	          ier=1 nfp is greater then 50
c	          ier=2 2*nfp+1 is greater then n
c
c	REMARKS
c	   * n must be geater or equal then 2*nfp+1
c	   * nfp must be smaller then 50
c
c	SUBROUTINES AND FUNCTION SUBPROGRAM REQUIRED
c	   none
c
c	METHOD
c	   uses digital filter described in R.W.HAMMING
c	   'DIGITAL FILTERS',PRENTICE-HALL,INC.
c	   ENGLEWOOD CLIFFS,NEW JERSEY,CHAPTER 6.3,PAGE 106.
c
c------------------------------------------------------------------------------
	subroutine NFSMT(fx,n,fc,nfp,ier)
	implicit real*8 (a-h,o-z)

	parameter (PI=3.141592654)
	real*8 filter(0:100),fx(*)

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
	  filter(nfp+i)=(sin(f2)/(PI*f1))*(sin(smarg)/smarg)
	  filter(nfp-i)=filter(nfp+i)
1	continue
	filter(nfp)=2.0*fc

c smoothing
	iend=n-2*nfp
	do 3 i=1,iend
	  sfx=0.0
	  do 4 j=0,2*nfp
	    sfx=sfx+filter(j)*fx(i+j)
4	  continue
	  fx(nfp+i)=sfx
3	continue

	ier=0
	return

	end
