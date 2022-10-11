	subroutine ORTPOL(n,ig,x0,xh,y,ysmooth,c)
	implicit real*8 (a-h,o-z)
c******************************************************************************
c smoothing procedure using ortogonal polynom aproximation
c n:in-integer*2 # of points
c ig:in-integer*2 polynoms grade
c x0:in-real*4 interval center
c xh:in-real*4 array of x values
c y:in-real*4 array of y values
c ysmooth:out-real*4 array of smoothed y values
c c:coefficients of polynom y=Co+C1x+...+Cigx^ig
c******************************************************************************

	parameter (grade=3,size=100)
	real*8 a(0:grade),c1(0:grade,0:grade),d(0:grade,grade+1),
     *	z(0:grade,size),x(size),xh(*),y(*),ysmooth(*),c(0:*)

	if(x0.ne.0) then
	  do 1 i=1,n
	    x(i)=xh(i)-x0
1	  continue
	end if

c 1st grade
	s1=0.0
	s2=0.0
	a(0)=0.0
	a(1)=0.0
	do 2 k=1,n
	  s1=s1+x(k)
	  s2=s2+x(k)*x(k)
2	continue
	c1(0,0)=1.0
	rnenn=sqrt(n*s2-s1*s1)
	c1(0,1)=-s1/rnenn
	c1(1,1)=n/rnenn
	do 3 k=1,n
	  z(0,k)=c1(0,0)
	  z(1,k)=c1(0,1)+c1(1,1)*x(k)
	  a(0)=a(0)+y(k)
	  a(1)=a(1)+y(k)*z(1,k)
3	continue
	a(0)=a(0)/n
	a(1)=a(1)/n
	do 4 k=1,n
	  ysmooth(k)=a(0)+a(1)*z(1,k)
4	continue
	c(0)=a(0)*c1(0,0)+a(1)*c1(0,1)
	c(1)=a(1)*c1(1,1)

c 2nd > IGth grade
	if(ig.le.1) return
	do 5 i=2,ig
	  a(i)=0.0
	  alfa1=0.0
	  alfa2=0.0
	  alfa3=0.0
	  c(i)=0.0
	  do 6 k=1,n
	    alfa1=alfa1+(x(k)*x(k))*(z(i-1,k))**2
	    alfa2=alfa2+x(k)*(z(i-1,k))**2
	    alfa3=alfa3+x(k)*z(i-1,k)*z(i-2,k)
6	  continue
	  rnenn=sqrt(abs(n*alfa1-alfa2**2-alfa3**2))
	  alfa1=n/rnenn
	  alfa2=-alfa2/rnenn
	  alfa3=-alfa3/rnenn
c z(i,k) ith grade polynom in kth point computing recursiv 
c from i-1st and i-2nd grade polynom
c ysmooth=AoZo+A1Z1+.....+AiZi+.....+AmaxgradeZmaxgrade
	  do 7 k=1,n
	    z(i,k)=(alfa1*x(k)+alfa2)*z(i-1,k)+alfa3*z(i-2,k)
	    a(i)=a(i)+y(k)*z(i,k)
7	  continue
	  a(i)=a(i)/n
c ysmooth(i)=ysmooth(i-1)+A(i)Z(i)
	  do 8 k=1,n
	    ysmooth(k)=ysmooth(k)+a(i)*z(i,k)
8	  continue
	  d(0,1)=0.0
	  d(0,2)=c1(0,i-1)
	  d(0,3)=c1(0,i-2)
	  do 9 k=1,i
	    d(k,1)=c1(k-1,i-1)
	    d(k,2)=c1(k,i-1)
	    d(k,3)=c1(k,i-2)
9	  continue
	  d(i-1,3)=0.0
	  d(i,2)=0.0
	  d(i,3)=0.0
	  do 10 k=0,i
	    c1(k,i)=alfa1*d(k,1)+alfa2*d(k,2)+alfa3*d(k,3)
10	  continue
c 
	  do 11 k=0,i
	    c(k)=c(k)+a(i)*c1(k,i)
11	  continue
5	continue

	return
	end
