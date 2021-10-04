c***********************************************************************
      subroutine  linear( x,h,ram,ee,k,ig,funct )
c
c     this subroutine performs the linear search along the direction spe
c     by the vector h
c       ----------------------------------------------------------------
c       the following subroutine is directly called by this subroutine:
c             funct
c       ----------------------------------------------------------------
c
c     inputs:
c        x:       vector of position
c        h:       search direction
c        k:       dimension of vector x
c
c     outputs:
c        ram:     optimal step width
c        e2:      minimum function value
c        ig:      error code
c
      implicit  real  *8 ( a-h,o-z )
      external funct
      integer  return,sub
      dimension  x(1) , h(1) , x1(50)
      dimension  g(50)
      common     / ccc /  isw , ipr
      common /mpi/nprocs, myrank, ista, iend, ista2, iend2,ised
c
      isw = 1
      if( ram .le. 1.0d-30 )  ram = 0.1d0
      const2 = 1.0d-16
      hnorm = 0.d0
      do 10  i=1,k
   10 hnorm = hnorm + h(i)**2
      hnorm = dsqrt( hnorm )

     
      if(hnorm.gt.1)ram = ram/hnorm

      ram2=ram
      e1 =ee
      ram1 = 0.d0
c
      do 20  i=1,k
   20 x1(i) = x(i) + ram2*h(i)
      call  funct( k,x1,e2,g,ig )
      if(ipr.ge.7.and.myrank.eq.0) write(6,2)  ram2,e2
c
       if( ig .eq. 1 )  go to  50
       if( e2 .gt. e1 ) go to 50
   30 ram3 = ram2*2.d0
      do 40  i=1,k
   40 x1(i) = x(i) + ram3*h(i)
      call  funct( k,x1,e3,g,ig )
      if( ig.eq.1 )  go to  500
      if( ipr.ge.7.and.myrank.eq.0 ) write(6,3)  ram3,e3
      if( e3 .gt. e2 )  go to 70
      ram1 = ram2
      ram2 = ram3
      e1 = e2
      e2 = e3
      go to 30
c
   50 ram3 = ram2
      e3 = e2
      ram2 = ram3*0.1d0
      if( ram2*hnorm .lt. const2 )  go to  400
      do 60  i=1,k
   60 x1(i) = x(i) + ram2*h(i)
      call  funct( k,x1,e2,g,ig )
      if(ipr.ge.7.and.myrank.eq.0) write(6,4)  ram2,e2
      if( e2.gt.e1 )  go to 50
c
   70 ret = 80
      go to 200
c
   80 do 90  i=1,k
   90 x1(i) = x(i) + ram*h(i)
      call  funct( k,x1,ee,g,ig )
      if(ipr.ge.7.and.myrank.eq.0)  write(6,5)  ram,ee
c
      ifg = 0
      sub = 300
      sub =200

   95 ret = 130
      if( ram .gt. ram2 )  go to 110
      if( ee .ge. e2 )  go to 100
      ram3 = ram2
      ram2 = ram
      e3 =e2
      e2 =ee
      if(sub.eq.200) goto 200
      if(sub.eq.300) goto 300
c
  100 ram1 = ram
      e1 = ee
      if(sub.eq.200) goto 200
      if(sub.eq.300) goto 300
c
  110 if( ee .le. e2 )  go to 120
      ram3 = ram
      e3 = ee
      if(sub.eq.200) goto 200
      if(sub.eq.300) goto 300
 
c
  120 ram1 = ram2
      ram2 = ram
      e1 = e2
      e2 = ee
      if(sub.eq.200) goto 200
      if(sub.eq.300) goto 300

c
  130 do 140  i=1,k
  140 x1(i) = x(i) + ram*h(i)
      call  funct( k,x1,ee,g,ig )
      if( ipr.ge.7.and.myrank.eq.0 )  write(6,6)  ram,ee
       sub=200
      ifg = ifg+1
      ifg = 0
      if( ifg .eq. 1 )  go to 95
c
      if( e2 .lt. ee )  ram = ram2
      return
c
c      -------  internal subroutine sub1  -------
  200 a1 = (ram3-ram2)*e1
      a2 = (ram1-ram3)*e2
      a3 = (ram2-ram1)*e3
      b2 = (a1+a2+a3)*2.d0
      b1 = a1*(ram3+ram2) + a2*(ram1+ram3) + a3*(ram2+ram1)
      if( b2 .eq. 0.d0 )  go to 210
      ram = b1 /b2
      if(ret.eq.80) goto 80
      if(ret.eq.130) goto 130
  210 ig = 1
      ram = ram2
      return
c
c      -------  internal subroutine sub2  -------
c
  300 if( ram3-ram2 .gt. ram2-ram1 )  go to 310
      ram = (ram1+ram2)*0.5d0
      if(ret.eq.80) goto 80
      if(ret.eq.130) goto 130
c
  310 ram = (ram2+ram3)*0.5d0
      if(ret.eq.80) goto 80
      if(ret.eq.130) goto 130
c ------------------------------------------------------------
c
  400 ram = 0.d0
      return
c ------------------------------------------------------------
c
  500 ram = (ram2+ram3)*0.5d0
  510 do 520  i=1,k
  520 x1(i) = x(i) + ram*h(i)
      call  funct( k,x1,e3,g,ig )
      if( ipr.ge.7.and.myrank.eq.0 )  write(6,7)  ram,e3
      if( ig.eq.1 )  go to 540
      if( e3.gt.e2 )  go to 530
      ram1 = ram2
      ram2 = ram
      e1 = e2
      e2 = e3
      go to 500
c
  530 ram3 = ram
      go to 70
c
  540 ram = (ram2+ram)*0.5d0
      go to 510
c
c ------------------------------------------------------------
    1 format( 1h ,'lambda =',d18.10, 10x,'e1 =',d25.17 )
    2 format( 1h ,'lambda =',d18.10, 10x,'e2 =',d25.17 )
    3 format( 1h ,'lambda =',d18.10, 10x,'e3 =',d25.17 )
    4 format( 1h ,'lambda =',d18.10, 10x,'e4 =',d25.17 )
    5 format( 1h ,'lambda =',d18.10, 10x,'e5 =',d25.17 )
    6 format( 1h ,'lambda =',d18.10, 10x,'e6 =',d25.17 )
    7 format( 1h ,'lambda =',d18.10, 10x,'e7 =',d25.17 )
      end

c************************************************************************
C************************************************************************

      subroutine  davidn( x,n,ihes,funct )
c
c          minimization by davidon-fletcher-powell procedure
c
c       ----------------------------------------------------------------
c       the following subroutines are directly called by this subroutine
c             funct
c             hesian
c             linear
c       ----------------------------------------------------------------
c          inputs:
c             x:       vector of initial values
c             k:       dimension of the vector x
c             ihes:    =0   inverse of hessian matrix is not available
c                      =1   inverse of hessian matrix is available
c
c          output:
c             x:       vector of minimizing solution
c
      implicit  real * 8  ( a-h , o-z )
      external funct
      dimension  x(50) , dx(50) , g(50) , g0(50) , y(50)
      dimension  h(50,50) , wrk(50) , s(50)
      dimension  r(31,31)
      common     / ccc /  isw
      common     / ddd /  r , f , aic , sd
      common /mpi/nprocs, myrank, ista, iend, ista2, iend2,ised

      data  tau1 , tau2  /  1.0d-6 , 1.0d-6  /
      data  eps1 , eps2  / 1.0d-6 , 1.0d-6  /
      iter=0
      ramda = 0.05d0
      const1 = 1.0d-17
c
c          initial estimate of inverse of hessian
c
      do  20   i=1,n
      do  10   j=1,n
   10 h(i,j) = 0.0d00
      s(i) = 0.0d00
      dx(i) = 0.0d00
   20 h(i,i) = 1.0d00
      isw = 0
c
      call  funct( n,x,xm,g,ig )
c
      if(myrank.eq.0)write( 6,340 )     xm , 2*(xm+n)
c
c          inverse of hessian computation (if available)
c
c      if( ihes .eq. 1 )   call  hesian( x,n,h )
c
      icc = 0
c      iteration
11110 continue
      icc = icc + 1
      do  11111   ic=1,n
      if( ic .eq. 1 .and. icc .eq. 1 )     go to 120
c
      do  40   i=1,n
   40 y(i) = g(i) - g0(i)
      do  60   i=1,n
      sum = 0.0d00
      do  50   j=1,n
   50 sum = sum + y(j) * h(i,j)
   60 wrk(i) = sum
      s1 = 0.0d00
      s2 = 0.0d00
      do  70   i=1,n
      s1 = s1 + wrk(i) * y(i)
   70 s2 = s2 + dx(i) * y(i)
      if( s1.le.const1 .or. s2.le.const1 )  go to 900
      if( s1 .le. s2 )     go to 100
c
c          update the inverse of hessian matrix
c
c               ---  davidon-fletcher-powell type correction  ---
c
      do  90   i=1,n
      do  90   j=i,n
      h(i,j) = h(i,j) + dx(i)*dx(j)/s2 - wrk(i)*wrk(j)/s1
   90 h(j,i) = h(i,j)
      go to  120
c
c               ---  fletcher type correction  ---
c
  100 continue
      stem = s1 / s2 + 1.0d00
      do  110   i=1,n
      do  110   j=i,n
      h(i,j) = h(i,j)- (dx(i)*wrk(j)+wrk(i)*dx(j)-dx(i)*dx(j)*stem)/s2
  110 h(j,i) = h(i,j)
c
c
c
  120 continue
      ss = 0.0d00
      do  150   i=1,n
      sum = 0.0d00
      do  140   j=1,n
  140 sum = sum + h(i,j)*g(j)
      ss = ss + sum * sum
  150 s(i) = -sum
c
c
      s1 = 0.0d00
      s2 = 0.0d00
      do  170   i=1,n
      s1 = s1 + s(i)*g(i)
  170 s2 = s2 + g(i)*g(i)
      ds2 = dsqrt(s2)
      gtem = dabs(s1) / ds2
c     write(6,610)gtem,ds2
      if( gtem .le. tau1  .and.  ds2 .le. tau2 )     go to  900
      if( s1 .lt. 0.0d00 )     go to  200
      do  190   i=1,n
      do  180   j=1,n
  180 h(i,j) = 0.0d00
      h(i,i) = 1.0d00
  190 s(i) = -s(i)
  200 continue
c
      ed = xm
c
c          linear  search
c
      call  linear( x,s,ramda,ed,n,ig,funct )
c
      if(myrank.eq.0)write( 6,330 )     ramda , ed , s1 , s2
c
      iter=iter+1
c     if(mod(iter,10).eq.0) call hes4(n,x,h)
c     if(mod(iter,10).eq.0) call invdet(h,tdet,n,50)
c
      s1 = 0.0d00
      do  210   i=1,n
      dx(i) = s(i) * ramda
      s1 = s1 + dx(i) * dx(i)
      g0(i) = g(i)
  210 x(i) = x(i) + dx(i)
      xmb = xm
      isw = 0
c
      call  funct( n,x,xm,g,ig )
c
      s2 = 0.d0
      do  220     i=1,n
  220 s2 = s2 + g(i)*g(i)
      if( dsqrt(s2) .gt. tau2 )   go to  11111
      if( xmb/xm-1.d0 .lt. eps1  .and.  dsqrt(s1) .lt. eps2 )  go to 900
11111 continue
      if( icc .ge. 10 )     go to 900
      go to 11110
  900 continue
c
      if(myrank.eq.0)then
         write( 6,340 )     xm , 2*(xm+n)
c
         write( 6,600 )
         write( 6,610 )     (x(i),i=1,n)
         write( 6,601 )
         write( 6,610 )     (g(i),i=1,n)
       endif

c@@@
c     call hes4(n,x,h)
c     do 653 i=1,n
c 653 write(6,604) (h(i,j),j=1,n)
c     call invdet(h,hdet,n,50)
c     write(6,602)
c     do 652 i=1,n
c 652 write(6,604) (h(i,j),j=1,n)

  602 format(1h ,'*** estimated inverse hessian ***')
  604 format(1h ,8d15.5/(1h ,8d15.5))
c
      return
  330 format( 1h ,' lmbd = ',d13.7,2x,'-ll = ',d21.15,2x,d9.2,1x,d9.2)
  340 format( 1h ,'- log likelihood =',d23.15,5x,'aic =',f10.1)
  600 format( 1h ,'-----  x  -----' )
  601 format( 1h ,'***  gradient  ***' )
  610 format( 1h ,10d13.5 )
      end


     
    


